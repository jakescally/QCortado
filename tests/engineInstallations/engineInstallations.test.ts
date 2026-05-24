import assert from "node:assert/strict";
import test from "node:test";

import {
  INSTALLABLE_ENGINE_DESCRIPTORS,
  QE_ENGINE_DESCRIPTOR,
  buildDefaultEngineInstallForm,
  isEngineAlreadyAvailable,
  isSelectableEngineStatus,
} from "../../src/lib/engines";
import type { HpcProfile } from "../../src/lib/types";

const profile: HpcProfile = {
  id: "andromeda",
  name: "Andromeda",
  cluster: "andromeda",
  host: "andromeda.example.edu",
  port: 22,
  username: "researcher",
  auth_method: "ssh_key",
  remote_qe_bin_dir: "/opt/qe/bin",
  remote_pseudo_dir: "/opt/qe/pseudo",
  remote_workspace_root: "/scratch/researcher/qcortado",
  remote_project_root: "/project/researcher/qcortado",
  resource_mode: "both",
  launcher: "srun",
  default_cpu_resources: {
    resource_type: "cpu",
  },
  default_gpu_resources: {
    resource_type: "gpu",
  },
  credential_persisted: true,
  warning_acknowledged: true,
  created_at: "2026-01-01T00:00:00Z",
  updated_at: "2026-01-01T00:00:00Z",
};

test("engine install candidates include QE and hidden WIEN2k", () => {
  assert.deepEqual(
    INSTALLABLE_ENGINE_DESCRIPTORS.map((descriptor) => descriptor.id),
    ["qe", "wien2k"],
  );
});

test("already selectable engines are not installable again", () => {
  assert.equal(isEngineAlreadyAvailable("qe", [QE_ENGINE_DESCRIPTOR]), true);
  assert.equal(isEngineAlreadyAvailable("wien2k", [QE_ENGINE_DESCRIPTOR]), false);
  assert.equal(
    isEngineAlreadyAvailable("wien2k", [QE_ENGINE_DESCRIPTOR], [
      {
        engineId: "wien2k",
        hpcProfileId: "andromeda",
        remoteInstallRoot: "/opt/WIEN2k",
        remoteWorkspaceRoot: "/scratch/researcher/qcortado",
        remoteProjectRoot: "/project/researcher/qcortado",
        verifiedExecutables: ["x", "init_lapw", "run_lapw", "runsp_lapw"],
        verifiedAt: "2026-01-01T00:00:00Z",
      },
    ]),
    true,
  );
});

test("configured engine status is selectable but reserved is not", () => {
  assert.equal(isSelectableEngineStatus("implemented"), true);
  assert.equal(isSelectableEngineStatus("configured"), true);
  assert.equal(isSelectableEngineStatus("reserved"), false);
});

test("engine install form defaults only choose the selected HPC profile", () => {
  const defaults = buildDefaultEngineInstallForm(profile);

  assert.equal(defaults.hpcProfileId, "andromeda");
  assert.equal(defaults.remoteInstallRoot, "");
  assert.equal(Object.hasOwn(defaults, "remoteWorkspaceRoot"), false);
  assert.equal(Object.hasOwn(defaults, "remoteProjectRoot"), false);
});
