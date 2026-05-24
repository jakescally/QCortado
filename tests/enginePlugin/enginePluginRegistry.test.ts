import assert from "node:assert/strict";
import test from "node:test";
import {
  DEFAULT_ENGINE_ID,
  ENGINE_WORKFLOW_VIEWS,
  QE_ENGINE_PLUGIN_MANIFEST,
  FALLBACK_ENGINE_PLUGIN_MANIFESTS,
  getEngineWorkflowView,
  getFrontendEnginePlugin,
  isEngineWorkflowView,
  resolveEngineWorkflowHostRoute,
} from "../../src/lib/engines";
import { WIEN2K_RESERVED_ENGINE_PLUGIN_MANIFEST } from "../../src/lib/engines/wien2k";

test("QE frontend plugin is the default engine and exposes workflow routes", () => {
  const plugin = getFrontendEnginePlugin(DEFAULT_ENGINE_ID);

  assert.ok(plugin);
  assert.equal(plugin.id, "qe");
  assert.equal(plugin.getWorkflowView("scf"), "scf-wizard");
  assert.equal(plugin.getWorkflowView("bands"), "bands-wizard");
  assert.equal(getEngineWorkflowView("qe", "transport"), "transport-wizard");
});

test("QE plugin manifest separates shared panels from QE-specific panels", () => {
  const scf = QE_ENGINE_PLUGIN_MANIFEST.workflows.find((workflow) => workflow.kind === "scf");
  const bands = QE_ENGINE_PLUGIN_MANIFEST.workflows.find((workflow) => workflow.kind === "bands");

  assert.ok(scf);
  assert.ok(bands);
  assert.ok(scf.sharedPanels.some((panel) => panel.kind === "structure_viewer"));
  assert.ok(scf.sharedPanels.some((panel) => panel.kind === "hpc_run_settings"));
  assert.ok(scf.enginePanels.some((panel) => panel.componentKey === "qe.scf.pseudopotentials"));
  assert.ok(bands.produces.includes("band_dataset"));
});

test("engine workflow host routes every QE workflow through a supported wizard view", () => {
  for (const workflow of QE_ENGINE_PLUGIN_MANIFEST.workflows) {
    const route = resolveEngineWorkflowHostRoute(DEFAULT_ENGINE_ID, workflow.kind);

    assert.ok(route, `missing route for ${workflow.kind}`);
    assert.equal(route.engineId, DEFAULT_ENGINE_ID);
    assert.equal(route.kind, workflow.kind);
    assert.ok(isEngineWorkflowView(route.view));
    assert.ok(ENGINE_WORKFLOW_VIEWS.includes(route.view));
  }
});

test("Wien2k skeleton is reserved and hidden from implemented frontend routes", () => {
  assert.equal(getFrontendEnginePlugin("wien2k"), null);
  assert.equal(getEngineWorkflowView("wien2k", "scf"), null);
  assert.equal(resolveEngineWorkflowHostRoute("wien2k", "scf"), null);
  assert.deepEqual(
    FALLBACK_ENGINE_PLUGIN_MANIFESTS.map((manifest) => manifest.descriptor.id),
    ["qe"],
  );

  assert.equal(WIEN2K_RESERVED_ENGINE_PLUGIN_MANIFEST.descriptor.status, "reserved");
  assert.deepEqual(WIEN2K_RESERVED_ENGINE_PLUGIN_MANIFEST.descriptor.executionModes, ["remote"]);
  assert.equal(WIEN2K_RESERVED_ENGINE_PLUGIN_MANIFEST.hpc.supportsRemoteOnly, true);
  assert.equal(WIEN2K_RESERVED_ENGINE_PLUGIN_MANIFEST.hpc.supportsLocal, false);
});

test("Wien2k reserved manifest models case state without QE pseudopotentials", () => {
  const scf = WIEN2K_RESERVED_ENGINE_PLUGIN_MANIFEST.workflows.find((workflow) => workflow.kind === "scf");

  assert.ok(scf);
  assert.ok(scf.inputRequirements.some((requirement) => requirement.kind === "engine_case_state"));
  assert.ok(scf.enginePanels.every((panel) => !panel.componentKey.includes("pseudo")));
  assert.ok(scf.produces.includes("scf_summary"));
});
