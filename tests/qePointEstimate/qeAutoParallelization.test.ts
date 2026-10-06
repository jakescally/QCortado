import assert from "node:assert/strict";
import test from "node:test";
import { applyQeAutoPlan, planQeAutoAllocation, sourceScfLayout } from "../../src/lib/qeAutoParallelization";
import { sourceScfWorkload } from "../../src/lib/qePointEstimate";

test("phonon Auto favors images and one electronic pool when each image is small", () => {
  const plan = planQeAutoAllocation({ nodes: 3, totalCpus: 24, phonon: true, qPoints: 12, kPoints: 20 });
  assert.equal(plan.images, 3);
  assert.equal(plan.pools, 1);
  assert.equal(plan.ranksPerImage, 8);
  assert.equal(plan.threads, 1);
});

test("SCF Auto prefers node-local pools without imposing an eight-rank cap", () => {
  const plan = planQeAutoAllocation({ nodes: 2, totalCpus: 64, phonon: false, kPoints: 16 });
  assert.equal(plan.images, 1);
  assert.equal(plan.pools, 2);
  assert.equal(plan.ranksPerPool, 32);
});

test("known compatible source SCF layout is preserved", () => {
  const plan = planQeAutoAllocation({ nodes: 2, totalCpus: 64, phonon: true, qPoints: 10, kPoints: 20, sourceLayout: { ranks: 8, pools: 2 } });
  assert.equal(plan.ranksPerImage, 8);
  assert.ok(plan.images <= 10);
  assert.equal(plan.pools, 2);
  assert.equal(plan.ranksPerPool, 4);
});

test("small or unknown meshes never create pools or images without useful points", () => {
  for (const phonon of [false, true]) {
    const plan = planQeAutoAllocation({ nodes: 1, totalCpus: 32, phonon, qPoints: 1, kPoints: 1 });
    assert.equal(plan.images, 1);
    assert.equal(plan.pools, 1);
  }
  const unknown = planQeAutoAllocation({ nodes: 2, totalCpus: 16, phonon: true });
  assert.equal(unknown.images, 1);
  assert.equal(unknown.pools, 1);
  assert.ok(unknown.notes.some((note) => note.includes("unavailable")));
});

test("rank caps use integer CPU factors without losing requested CPUs", () => {
  const plan = planQeAutoAllocation({ nodes: 2, totalCpus: 64, phonon: false, kPoints: 16, maxRanks: 20 });
  assert.equal(plan.ranks, 16);
  assert.equal(plan.threads, 4);
  assert.equal(plan.ranks * plan.threads, 64);
  assert.throws(() => planQeAutoAllocation({ nodes: 4, totalCpus: 17, phonon: false, maxRanks: 12 }));
});

test("uneven node allocations let Slurm distribute ranks", () => {
  const plan = planQeAutoAllocation({ nodes: 3, totalCpus: 25, phonon: true, qPoints: 12, kPoints: 20 });
  assert.equal(plan.ranks, 25);
  assert.equal(plan.tasksPerNode, null);
  assert.equal(plan.ranks % (plan.images * plan.pools), 0);
});

test("a range of allocation sizes always produces valid decompositions", () => {
  for (const nodes of [1, 2, 3, 4]) for (const cpus of [4, 8, 12, 24, 32, 64, 96, 128]) {
    for (const kPoints of [1, 3, 10, 40]) for (const qPoints of [1, 5, 12, 30]) {
      const plan = planQeAutoAllocation({ nodes, totalCpus: cpus, phonon: true, kPoints, qPoints });
      assert.equal(plan.ranks * plan.threads, cpus);
      assert.equal(plan.ranks % (plan.images * plan.pools), 0);
      assert.ok(plan.images <= qPoints);
      assert.ok(plan.pools <= kPoints);
      assert.ok(plan.ranks >= nodes);
    }
  }
});

test("concrete Auto resources preserve unrelated cluster settings", () => {
  const plan = planQeAutoAllocation({ nodes: 3, totalCpus: 24, phonon: true, qPoints: 12, kPoints: 20 });
  const request = applyQeAutoPlan({ resource_type: "cpu", account: "project", walltime: "01:00:00", memory_gb: 64, ntasks_per_node: 99, pencil_decomposition: true }, plan);
  assert.equal(request.ntasks, 24);
  assert.equal(request.ntasks_per_node, 8);
  assert.equal(request.qe_total_cpus, 24);
  assert.equal(request.qe_qpoint_images, 3);
  assert.equal(request.account, "project");
  assert.equal(request.memory_gb, 64);
  assert.equal(request.pencil_decomposition, false);
});

test("legacy calculations without layout metadata remain supported", () => {
  assert.equal(sourceScfLayout({}), undefined);
  assert.equal(sourceScfLayout({ qe_point_parallelization: { mpi_tasks: 10, kpoint_pools: 3 } }), undefined);
  assert.deepEqual(sourceScfLayout({ qe_point_parallelization: { mpi_tasks: 8, kpoint_pools: 2 } }), { ranks: 8, pools: 2 });
});

test("800 CPUs on 20 nodes use hybrid parallelism and respect point capacity", () => {
  for (const phonon of [false, true]) {
    const plan = planQeAutoAllocation({ nodes: 20, totalCpus: 800, phonon, kPoints: 20, qPoints: 12 });
    assert.ok(plan.ranks < 800);
    assert.ok(plan.threads >= 2 && plan.threads <= 8);
    assert.equal(plan.ranks * plan.threads, 800);
    assert.equal(plan.ranks % 20, 0);
    assert.ok(plan.ranksPerPool * plan.threads <= 40);
    assert.equal(plan.ranks % (plan.images * plan.pools), 0);
  }
});

test("a large allocation with too little independent work warns instead of creating empty groups", () => {
  const plan = planQeAutoAllocation({ nodes: 20, totalCpus: 800, phonon: true, kPoints: 1, qPoints: 1 });
  assert.equal(plan.images, 1);
  assert.equal(plan.pools, 1);
  assert.ok(plan.ranks < 800);
  assert.equal(plan.ranks * plan.threads, 800);
  assert.ok(plan.notes.some((note) => note.includes("exceed useful parallelism")));
});

test("hybrid candidates preserve source SCF MPI layout when pure MPI cannot", () => {
  const plan = planQeAutoAllocation({ nodes: 20, totalCpus: 800, phonon: true, kPoints: 20, qPoints: 12,
    sourceLayout: { ranks: 8, pools: 2 } });
  assert.equal(plan.ranksPerImage, 8);
  assert.equal(plan.pools, 2);
  assert.ok(plan.images <= 12);
  assert.equal(plan.ranks * plan.threads, 800);
  assert.ok(plan.notes.some((note) => note.includes("Preserving the source")));
});

test("hybrid layouts preserve full budgets and rank caps across uneven node allocations", () => {
  for (const nodes of [1, 3, 20]) for (const totalCpus of [128, 255, 800, 1024]) {
    const plan = planQeAutoAllocation({ nodes, totalCpus, phonon: true, kPoints: 7, qPoints: 11, maxRanks: 128 });
    assert.equal(plan.ranks * plan.threads, totalCpus);
    assert.ok(plan.ranks >= nodes && plan.ranks <= 128);
    assert.equal(plan.ranks % (plan.images * plan.pools), 0);
    assert.ok(plan.images <= 11 && plan.pools <= 7);
    if (plan.tasksPerNode != null) assert.equal(plan.tasksPerNode * nodes, plan.ranks);
  }
});

test("known irreps permit images beyond q-points without exceeding actual work", () => {
  const input = { nodes: 20, totalCpus: 800, phonon: true, kPoints: 20, qPoints: 1 };
  assert.equal(planQeAutoAllocation(input).images, 1);
  const plan = planQeAutoAllocation({ ...input, irreps: 24 });
  assert.ok(plan.images > 1 && plan.images <= 24);
  assert.equal(plan.ranks % (plan.images * plan.pools), 0);
  assert.throws(() => planQeAutoAllocation({ ...input, irreps: 0 }));
});

test("FFT capacity is honored without an arbitrary eight-rank ceiling", () => {
  const plan = planQeAutoAllocation({ nodes: 2, totalCpus: 64, phonon: false, kPoints: 1, fftPlanes: 16 });
  assert.ok(plan.ranksPerPool > 8 && plan.ranksPerPool <= 16);
  assert.ok([1, 2, 4].includes(plan.threads));
  assert.equal(plan.ranks * plan.threads, 64);
  assert.throws(() => planQeAutoAllocation({ nodes: 20, totalCpus: 800, phonon: false, kPoints: 1, fftPlanes: 16 }));
});

test("source restart pools are preserved even before k estimates arrive", () => {
  const plan = planQeAutoAllocation({ nodes: 2, totalCpus: 64, phonon: true, qPoints: 1, irreps: 24,
    sourceLayout: { ranks: 8, pools: 2 } });
  assert.equal(plan.ranksPerImage, 8);
  assert.equal(plan.pools, 2);
  assert.ok(plan.images <= 24);
});

test("large regular allocations prefer benchmarked thread candidates and node-local groups", () => {
  for (const totalCpus of [256, 512, 800, 1024]) {
    const nodes = totalCpus === 800 ? 20 : 8;
    const plan = planQeAutoAllocation({ nodes, totalCpus, phonon: false, kPoints: 64 });
    assert.ok([1, 2, 4].includes(plan.threads));
    assert.ok(plan.ranksPerPool * plan.threads <= totalCpus / nodes);
    assert.equal((plan.ranks / nodes) % plan.ranksPerPool, 0);
  }
});

test("SCF output supplies the smaller FFT plane count and actual band count", () => {
  assert.deepEqual(sourceScfWorkload("Dense grid: FFT dimensions: (60, 60, 180)\nSmooth grid: FFT dimensions: ( 40, 40, 120 )\nnumber of Kohn-Sham states= 96"),
    { fftPlanes: 120, bands: 96 });
  assert.deepEqual(sourceScfWorkload(null), { fftPlanes: undefined, bands: undefined });
});
