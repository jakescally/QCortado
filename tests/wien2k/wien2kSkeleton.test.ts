import assert from "node:assert/strict";
import test from "node:test";
import {
  DEFAULT_WIEN2K_INITIALIZATION_SETTINGS,
  DEFAULT_WIEN2K_SCF_RUN_SETTINGS,
  buildWien2kInitCommandPlan,
  buildWien2kScfCommandPlan,
  initializedWien2kCaseArtifacts,
  isTerminalWien2kCasePhase,
  normalizeWien2kCaseName,
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
