import assert from "node:assert/strict";
import test from "node:test";
import {
  buildHpcQeInputCommandLine,
  buildHpcQeRuntimeSetupLines,
  hpcProfileToEngineRuntimeProfile,
  qeHpcProfileToRuntimeProfile,
} from "../../src/lib/engines/qe";
import type { HpcProfile } from "../../src/lib/types";

const profile: HpcProfile = {
  id: "profile-1",
  name: "Test Cluster",
  cluster: "test",
  host: "hpc.example.edu",
  port: 22,
  username: "user",
  auth_method: "ssh_key",
  ssh_key_path: "~/.ssh/id_rsa",
  remote_qe_bin_dir: "/opt/qe/bin",
  remote_qe_cpu_bin_dir: "/opt/qe/cpu/bin",
  remote_qe_gpu_bin_dir: "/opt/qe/gpu/bin",
  remote_epw_path: " /opt/qe/gpu/bin/epw.x ",
  remote_wannier90_path: " wannier90.x ",
  remote_postw90_path: null,
  remote_pseudo_dir: "/opt/qe/pseudo",
  remote_cpu_pseudo_dir: "/opt/qe/cpu/pseudo",
  remote_gpu_pseudo_dir: "/opt/qe/gpu/pseudo",
  remote_workspace_root: "/scratch/qcortado/work",
  remote_project_root: "/scratch/qcortado/projects",
  resource_mode: "both",
  launcher: "srun",
  launcher_cpu_extra_args: "--cpu-bind=cores",
  launcher_gpu_extra_args: "--gpu-bind=closest",
  launcher_extra_args: null,
  default_cpu_resources: {
    resource_type: "cpu",
    partition: "cpu",
    walltime: "02:00:00",
    nodes: 1,
    ntasks: 4,
    cpus_per_task: 1,
    memory_gb: 16,
    gpus: 0,
    additional_sbatch: ["--exclusive"],
  },
  default_gpu_resources: {
    resource_type: "gpu",
    partition: "gpu",
    walltime: "01:00:00",
    nodes: 1,
    ntasks: 1,
    cpus_per_task: 8,
    memory_gb: 32,
    gpus: 1,
    additional_sbatch: ["--gres=gpu:1"],
  },
  credential_persisted: false,
  warning_acknowledged: false,
  created_at: "2026-05-21T00:00:00.000Z",
  updated_at: "2026-05-21T00:00:00.000Z",
};

test("qeHpcProfileToRuntimeProfile resolves QE paths without changing the source profile", () => {
  const runtimeProfile = qeHpcProfileToRuntimeProfile(profile, {
    resourceType: "gpu",
  });

  assert.equal(runtimeProfile.engineId, "qe");
  assert.equal(runtimeProfile.supported, true);
  assert.equal(runtimeProfile.supportsResourceType, true);
  assert.equal(runtimeProfile.resourceType, "gpu");
  assert.equal(runtimeProfile.pathMode, "path");
  assert.equal(runtimeProfile.moduleLoad, null);
  assert.equal(runtimeProfile.launcherCommand, "srun --gpu-bind=closest");
  assert.equal(runtimeProfile.paths.remoteBinDir, "/opt/qe/gpu/bin");
  assert.equal(runtimeProfile.paths.remotePseudoDir, "/opt/qe/gpu/pseudo");
  assert.equal(runtimeProfile.paths.remoteEpwPath, "/opt/qe/gpu/bin/epw.x");
  assert.equal(runtimeProfile.paths.remoteWannier90Path, "wannier90.x");
  assert.equal(runtimeProfile.paths.remotePostw90Path, null);
  assert.equal(runtimeProfile.paths.usesPseudopotentials, true);
  assert.deepEqual(runtimeProfile.defaultResources.additional_sbatch, ["--gres=gpu:1"]);

  runtimeProfile.defaultResources.additional_sbatch?.push("--mutated");
  assert.deepEqual(profile.default_gpu_resources.additional_sbatch, ["--gres=gpu:1"]);
});

test("QE module mode loads its environment and invokes tools from PATH", () => {
  const moduleProfile: HpcProfile = {
    ...profile,
    qe_path_mode: "module",
    qe_module_use: "/cluster/modulefiles",
    qe_module_load: "quantum-espresso/7.5",
  };

  assert.deepEqual(
    buildHpcQeRuntimeSetupLines(moduleProfile, "cpu"),
    ["module use '/cluster/modulefiles'", "module load 'quantum-espresso/7.5'"],
  );
  assert.equal(
    buildHpcQeInputCommandLine(moduleProfile, "pw.x", "pw.in", "pw.out", undefined, "cpu"),
    "srun --cpu-bind=cores pw.x -pd .true. -in pw.in > pw.out 2>&1",
  );
  const runtimeProfile = qeHpcProfileToRuntimeProfile(moduleProfile);
  assert.equal(runtimeProfile.pathMode, "module");
  assert.equal(runtimeProfile.moduleLoad, "quantum-espresso/7.5");
});

test("hpcProfileToEngineRuntimeProfile does not synthesize future engine paths", () => {
  const runtimeProfile = hpcProfileToEngineRuntimeProfile(profile, "wien2k");

  assert.equal(runtimeProfile.engineId, "wien2k");
  assert.equal(runtimeProfile.supported, false);
  assert.equal(runtimeProfile.profileId, "profile-1");
  assert.match(runtimeProfile.reason, /does not have an HPC runtime profile adapter/i);
});
