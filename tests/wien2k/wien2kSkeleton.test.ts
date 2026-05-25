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
  validateWien2kStructureControls,
} from "../../src/lib/engines/wien2k";
import type { Wien2kCaseReference } from "../../src/lib/engines/wien2k";

const caseRef: Wien2kCaseReference = {
  caseName: "Si",
  remoteCaseDir: "/scratch/qcortado/projects/project-1/Si",
  remoteScratchDir: "/scratch/qcortado/work/Si",
  projectId: "project-1",
  cifId: "cif-1",
};

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
  assert.ok(basenames.includes("Si.scf"));
  assert.ok(basenames.includes("Si.klist"));
});

test("Wien2k command plans are data-only remote case plans", () => {
  const initPlan = buildWien2kInitCommandPlan(caseRef, DEFAULT_WIEN2K_INITIALIZATION_SETTINGS);
  const scfPlan = buildWien2kScfCommandPlan(caseRef, {
    ...DEFAULT_WIEN2K_SCF_RUN_SETTINGS,
    spinMode: "spin_polarized",
    parallel: true,
  });

  assert.equal(initPlan.program, "init_lapw");
  assert.equal(initPlan.workingDirectory, caseRef.remoteCaseDir);
  assert.ok(initPlan.argv.includes("-rkmax"));
  assert.deepEqual(initPlan.environment, [["SCRATCH", "/scratch/qcortado/work/Si"]]);

  assert.equal(scfPlan.program, "runsp_lapw");
  assert.ok(scfPlan.argv.includes("-p"));
  assert.equal(scfPlan.phase, "scf_running");
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
