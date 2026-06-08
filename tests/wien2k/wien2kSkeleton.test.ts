import assert from "node:assert/strict";
import test from "node:test";
import {
  DEFAULT_WIEN2K_INITIALIZATION_SETTINGS,
  DEFAULT_WIEN2K_SCF_RUN_SETTINGS,
  DEFAULT_WIEN2K_STRUCTURE_CONTROLS,
  buildWien2kInitCommandPlan,
  buildWien2kScfCommandPlan,
  initializedWien2kCaseArtifacts,
  isTerminalWien2kCasePhase,
  listWien2kStructureSources,
  normalizeWien2kCaseName,
  validateWien2kInitializationSettings,
  validateWien2kScfRunSettings,
  validateWien2kStructureControls,
} from "../../src/lib/engines/wien2k";
import {
  applyWien2kBandsTotalKPoints,
  getWien2kBandProjectionOptions,
  getWien2kBandProjectionOptionsFromSites,
  transformWien2kKPathForKlistBand,
} from "../../src/lib/wien2kBandsWizard";
import type { Wien2kCaseReference } from "../../src/lib/engines/wien2k";

const caseRef: Wien2kCaseReference = {
  caseName: "Si",
  remoteCaseDir: "/scratch/qcortado/projects/project-1/Si",
  remoteScratchDir: "/scratch/qcortado/work/Si",
  projectId: "project-1",
  cifId: "cif-1",
};

function crystalDataForSpaceGroup(spaceGroupNumber: number, overrides: Record<string, unknown> = {}) {
  return {
    cell_length_a: { value: 4 },
    cell_length_b: { value: 4 },
    cell_length_c: { value: 4 },
    cell_angle_alpha: { value: 90 },
    cell_angle_beta: { value: 90 },
    cell_angle_gamma: { value: 90 },
    space_group_IT_number: spaceGroupNumber,
    atom_sites: [],
    ...overrides,
  } as any;
}

test("Wien2k case names are filename prefixes, not paths", () => {
  assert.equal(normalizeWien2kCaseName(" Si_case-1.0 "), "Si_case-1.0");
  assert.equal(normalizeWien2kCaseName("../Si"), null);
  assert.equal(normalizeWien2kCaseName("Si case"), null);
  assert.equal(normalizeWien2kCaseName("Si/other"), null);
});

test("Wien2k initialized case artifacts preserve remote case-directory state", () => {
  const basenames = initializedWien2kCaseArtifacts("Si").map((artifact) => artifact.basename);

  assert.ok(basenames.includes("Si.struct"));
  assert.ok(basenames.includes("Si.clmsum"));
  assert.ok(!basenames.includes("Si.scf"));
  assert.ok(basenames.includes("Si.klist"));
  assert.ok(basenames.includes("Si.inc"));
  assert.ok(basenames.includes("Si.inm"));
});

test("Wien2k command plans are data-only remote case plans", () => {
  const initPlan = buildWien2kInitCommandPlan(caseRef, DEFAULT_WIEN2K_INITIALIZATION_SETTINGS);
  const ldaInitPlan = buildWien2kInitCommandPlan(caseRef, {
    ...DEFAULT_WIEN2K_INITIALIZATION_SETTINGS,
    exchangeCorrelation: 5,
  });
  const scfPlan = buildWien2kScfCommandPlan(caseRef, {
    ...DEFAULT_WIEN2K_SCF_RUN_SETTINGS,
    spinMode: "spin_polarized",
  }, true);

  assert.equal(initPlan.program, "init_lapw");
  assert.equal(initPlan.workingDirectory, caseRef.remoteCaseDir);
  assert.ok(initPlan.argv.includes("-rkmax"));
  assert.ok(!initPlan.argv.includes("-red"));
  assert.deepEqual(ldaInitPlan.argv.slice(ldaInitPlan.argv.indexOf("-vxc"), ldaInitPlan.argv.indexOf("-vxc") + 2), ["-vxc", "5"]);
  assert.deepEqual(initPlan.environment, [["SCRATCH", "/scratch/qcortado/work/Si"]]);

  assert.equal(scfPlan.program, "runsp_lapw");
  assert.ok(!scfPlan.argv.includes("-p"));
  assert.ok(scfPlan.argv.includes("-NI"));
  assert.equal(scfPlan.phase, "scf_running");
});

test("WIEN2k SCF v1 validates reviewed serial-native settings", () => {
  assert.equal(validateWien2kInitializationSettings(DEFAULT_WIEN2K_INITIALIZATION_SETTINGS), null);
  for (const exchangeCorrelation of [5, 11, 13, 19]) {
    assert.equal(
      validateWien2kInitializationSettings({ ...DEFAULT_WIEN2K_INITIALIZATION_SETTINGS, exchangeCorrelation }),
      null,
    );
  }
  assert.equal(validateWien2kScfRunSettings(DEFAULT_WIEN2K_SCF_RUN_SETTINGS), null);
  assert.match(
    validateWien2kInitializationSettings({ ...DEFAULT_WIEN2K_INITIALIZATION_SETTINGS, exchangeCorrelation: 999 }) ?? "",
    /native initialization options/,
  );
  assert.match(
    validateWien2kScfRunSettings({ ...DEFAULT_WIEN2K_SCF_RUN_SETTINGS, maxIterations: 0 }) ?? "",
    /positive integer/,
  );
});

test("WIEN2k SCF command plans expose DFT+U and advanced run flags", () => {
  const scfPlan = buildWien2kScfCommandPlan(caseRef, {
    ...DEFAULT_WIEN2K_SCF_RUN_SETTINGS,
    spinMode: "spin_polarized",
    iterativeDiagonalization: true,
    forceMinimization: true,
    dispersionCorrection: "dftd3",
    dftU: {
      enabled: true,
      doubleCounting: "sic",
      targets: [{
        siteIndex: 1,
        element: "Ni",
        manifold: "3d",
        orbitalL: 2,
        uEv: 6,
        jEv: 0,
      }],
    },
  });

  assert.equal(scfPlan.program, "runsp_lapw");
  assert.ok(scfPlan.argv.includes("-orb"));
  assert.ok(scfPlan.argv.includes("-dftd3"));
  assert.ok(scfPlan.argv.includes("-it"));
  assert.ok(scfPlan.argv.includes("-min"));
  assert.equal(validateWien2kScfRunSettings({
    ...DEFAULT_WIEN2K_SCF_RUN_SETTINGS,
    dftU: {
      enabled: true,
      doubleCounting: "sic",
      targets: [{
        siteIndex: 1,
        element: "Ni",
        manifold: "3d",
        orbitalL: 2,
        uEv: 6,
        jEv: 0,
      }],
    },
  })?.includes("requires spin-polarized"), true);
});

test("Wien2k terminal phases are explicit", () => {
  assert.equal(isTerminalWien2kCasePhase("initialized"), false);
  assert.equal(isTerminalWien2kCasePhase("scf_complete"), true);
  assert.equal(isTerminalWien2kCasePhase("failed"), true);
});

test("WIEN2k structure controls validate native refinement limits", () => {
  assert.equal(validateWien2kStructureControls(DEFAULT_WIEN2K_STRUCTURE_CONTROLS), null);
  assert.match(
    validateWien2kStructureControls({
      ...DEFAULT_WIEN2K_STRUCTURE_CONTROLS,
      nnBondlengthFactor: 0,
    }) ?? "",
    /NN bond-length factor/,
  );
  assert.match(
    validateWien2kStructureControls({
      ...DEFAULT_WIEN2K_STRUCTURE_CONTROLS,
      sgroupTolerance: 1e-8,
    }) ?? "",
    /SGROUP tolerance/,
  );
  assert.match(
    validateWien2kStructureControls({
      ...DEFAULT_WIEN2K_STRUCTURE_CONTROLS,
      siteOverrides: [{ siteIndex: 1, npt: 780 }],
    }) ?? "",
    /positive odd integer/,
  );
});

test("WIEN2k structure source adapter selects saved setup artifacts only", () => {
  const sources = listWien2kStructureSources([
    { id: "source", engine_id: "wien2k", calc_type: "engine_setup", parameters: { setup_kind: "structure" } },
    { id: "qe", engine_id: "qe", calc_type: "engine_setup", parameters: { setup_kind: "structure" } },
    { id: "future-scf", engine_id: "wien2k", calc_type: "scf", parameters: {} },
  ]);

  assert.deepEqual(sources.map((source) => source.id), ["source"]);
});

test("WIEN2k bands helper redistributes total k-points along the path", () => {
  const path = applyWien2kBandsTotalKPoints([
    { label: "G", coords: [0, 0, 0], npoints: 1 },
    { label: "X", coords: [0.5, 0, 0], npoints: 0 },
    { label: "L", coords: [0.5, 0.5, 0], npoints: 0 },
  ], 12);

  assert.equal(path.length, 3);
  assert.equal(path[0].npoints + path[1].npoints, 12);
  assert.equal(path[2].npoints, 0);
});

test("WIEN2k bands helper writes FCC path coordinates in WIEN2k conventional reciprocal basis", () => {
  const transformed = transformWien2kKPathForKlistBand([
    { label: "X", coords: [0.5, 0, 0.5], npoints: 0 },
  ], crystalDataForSpaceGroup(225));

  assert.deepEqual(transformed[0].coords, [0, 1, 0]);
});

test("WIEN2k bands helper writes BCC path coordinates in WIEN2k conventional reciprocal basis", () => {
  const transformed = transformWien2kKPathForKlistBand([
    { label: "H", coords: [0.5, -0.5, 0.5], npoints: 0 },
  ], crystalDataForSpaceGroup(229));

  assert.deepEqual(transformed[0].coords, [0, 1, 0]);
});

test("WIEN2k bands helper writes base-centered path coordinates in WIEN2k conventional reciprocal basis", () => {
  const transformed = transformWien2kKPathForKlistBand([
    { label: "A", coords: [0.25, 0.25, 0.5], npoints: 0 },
  ], crystalDataForSpaceGroup(65));

  assert.deepEqual(transformed[0].coords, [0.5, 0, 0.5]);
});

test("WIEN2k bands helper preserves rhombohedral primitive reciprocal coordinates", () => {
  const transformed = transformWien2kKPathForKlistBand([
    { label: "Z", coords: [0.5, 0.5, 0.5], npoints: 0 },
  ], crystalDataForSpaceGroup(166, {
    cell_length_a: { value: 4 },
    cell_length_b: { value: 4 },
    cell_length_c: { value: 10 },
    cell_angle_alpha: { value: 90 },
    cell_angle_beta: { value: 90 },
    cell_angle_gamma: { value: 120 },
  }));

  assert.deepEqual(transformed[0].coords, [0.5, 0.5, 0.5]);
});

test("WIEN2k bands projection options expose atom and orbital choices", () => {
  const options = getWien2kBandProjectionOptions({
    atom_sites: [
      { label: "Ni1", type_symbol: "Ni", fract_x: 0, fract_y: 0, fract_z: 0, occupancy: 1 },
      { label: "O1", type_symbol: "O", fract_x: 0.5, fract_y: 0.5, fract_z: 0.5, occupancy: 1 },
    ],
  } as any);

  assert.equal(options[0]?.value, "all");
  assert.ok(options.some((option) => option.value === "atom:1" && option.label === "Ni"));
  assert.ok(options.some((option) => option.value === "atom:1:orbital:2" && option.label === "Ni d"));
  assert.ok(options.some((option) => option.value === "atom:2:orbital:1" && option.label === "O p"));
});

test("WIEN2k bands projection options collapse expanded decorated sites to element groups", () => {
  const options = getWien2kBandProjectionOptions({
    atom_sites: Array.from({ length: 8 }, (_, index) => ({
      label: `Si0+_${index + 1}`,
      type_symbol: "Si0+",
      fract_x: 0,
      fract_y: 0,
      fract_z: 0,
      occupancy: 1,
    })),
  } as any);

  const atomOptions = options.filter((option) => option.kind === "atom");
  const orbitalOptions = options.filter((option) => option.kind === "orbital");

  assert.deepEqual(atomOptions.map((option) => option.label), ["Si"]);
  assert.deepEqual(orbitalOptions.map((option) => option.label), ["Si s", "Si p", "Si d", "Si f"]);
});

test("WIEN2k bands projection options prefer WIEN2k inequivalent site ordering", () => {
  const options = getWien2kBandProjectionOptionsFromSites([
    { symbol: "Ho", positions: [[0, 0, 0]] },
    { symbol: "P", positions: [[0.5, 0.5, 0.5]] },
  ] as any);

  const atomOptions = options.filter((option) => option.kind === "atom");

  assert.deepEqual(atomOptions.map((option) => [option.atomIndex, option.label]), [
    [1, "Ho"],
    [2, "P"],
  ]);
  assert.ok(options.some((option) => option.value === "jatom:1:orbital:3" && option.label === "Ho f"));
  assert.ok(options.some((option) => option.value === "jatom:2:orbital:1" && option.label === "P p"));
});

test("WIEN2K bands projection options keep repeated inequivalent element sites separate", () => {
  const options = getWien2kBandProjectionOptionsFromSites([
    { symbol: "Si", positions: [[0, 0, 0]] },
    { symbol: "Si", positions: [[0.25, 0.25, 0.25]] },
  ] as any);

  const atomOptions = options.filter((option) => option.kind === "atom");

  assert.deepEqual(atomOptions.map((option) => [option.atomIndex, option.label]), [
    [1, "Si"],
    [2, "Si"],
  ]);
});
