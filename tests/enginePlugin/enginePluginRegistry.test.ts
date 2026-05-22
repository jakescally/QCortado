import assert from "node:assert/strict";
import test from "node:test";
import {
  DEFAULT_ENGINE_ID,
  QE_ENGINE_PLUGIN_MANIFEST,
  getEngineWorkflowView,
  getFrontendEnginePlugin,
} from "../../src/lib/engines";

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
