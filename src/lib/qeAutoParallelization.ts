import type { SlurmResourceRequest } from "./types";

export interface QeSourceLayout { ranks: number; pools: number; }
export interface QeAutoPlan {
  nodes: number;
  totalCpus: number;
  ranks: number;
  threads: number;
  tasksPerNode: number | null;
  images: number;
  pools: number;
  ranksPerImage: number;
  ranksPerPool: number;
  notes: string[];
}

function positiveInteger(value: number, name: string): number {
  if (!Number.isSafeInteger(value) || value < 1 || value > 0xffffffff) {
    throw new Error(`${name} must be a positive integer no larger than 4,294,967,295.`);
  }
  return value;
}

function divisors(value: number): number[] {
  const values: number[] = [];
  for (let d = 1; d * d <= value; d++) {
    if (value % d) continue;
    values.push(d);
    if (d * d !== value) values.push(value / d);
  }
  return values.sort((a, b) => a - b);
}

export function sourceScfLayout(parameters: Record<string, unknown>): QeSourceLayout | undefined {
  const raw = parameters.qe_point_parallelization as Record<string, unknown> | undefined;
  const ranks = Number(raw?.mpi_tasks);
  const pools = Number(raw?.kpoint_pools);
  return Number.isSafeInteger(ranks) && Number.isSafeInteger(pools) && ranks > 0 && pools > 0 && ranks % pools === 0
    ? { ranks, pools } : undefined;
}

/** A transparent starting heuristic, not a machine-specific performance model. */
export function planQeAutoAllocation(input: {
  nodes: number;
  totalCpus: number;
  phonon: boolean;
  kPoints?: number | null;
  qPoints?: number | null;
  maxRanks?: number | null;
  sourceLayout?: QeSourceLayout;
  irreps?: number | null;
  fftPlanes?: number | null;
  bands?: number | null;
}): QeAutoPlan {
  const nodes = positiveInteger(input.nodes, "Nodes");
  const totalCpus = positiveInteger(input.totalCpus, "Total allocated CPUs");
  const rankLimit = input.maxRanks == null ? totalCpus : Math.min(totalCpus, positiveInteger(input.maxRanks, "MPI rank limit"));
  if (nodes > rankLimit) throw new Error("Use at least one MPI rank per requested node; reduce Nodes or increase the MPI rank limit.");
  const rankCandidates = divisors(totalCpus).filter((n) => n <= rankLimit && n >= nodes);
  if (!rankCandidates.length) throw new Error("This CPU allocation cannot satisfy the workflow's MPI rank limit. Adjust CPUs or Nodes.");
  const notes: string[] = [];
  const reliableCount = (value: number | null | undefined) => Number.isSafeInteger(value) && (value ?? 0) > 0 ? value! : null;
  const kPoints = reliableCount(input.kPoints);
  const qPoints = reliableCount(input.qPoints);
  if (kPoints == null || input.phonon && qPoints == null) notes.push("Point counts unavailable; using a conservative decomposition.");
  const source = input.sourceLayout;
  const validSource = source && Number.isSafeInteger(source.ranks) && Number.isSafeInteger(source.pools)
    && source.ranks > 0 && source.pools > 0 && source.ranks % source.pools === 0;
  const cpusPerNode = totalCpus / nodes;
  const irreps = input.irreps == null ? null : positiveInteger(input.irreps, "Known phonon irreps");
  const fftPlanes = reliableCount(input.fftPlanes);
  const bands = reliableCount(input.bands);
  const imageCapacity = input.phonon ? (irreps ?? qPoints ?? 1) : 1;
  const preferredThreads = nodes > 1 && totalCpus >= 128 ? 2 : 1;
  const imbalance = (points: number | null, count: number) => points == null ? 0 : 1 - points / (Math.ceil(points / count) * count);
  const candidates: { ranks: number; threads: number; images: number; pools: number;
    ranksPerImage: number; ranksPerPool: number; compatible: boolean; score: number }[] = [];
  for (const ranks of rankCandidates) {
    const threads = totalCpus / ranks;
    const standardThreads = [1, 2, 4].includes(threads);
    for (const images of divisors(ranks).filter((count) => count <= imageCapacity)) {
      const ranksPerImage = ranks / images;
      for (const pools of divisors(ranksPerImage).filter((count) => count <= (kPoints ?? 1) || Boolean(input.phonon && validSource && ranksPerImage === source!.ranks && count === source!.pools))) {
        const ranksPerPool = ranksPerImage / pools;
        // With pencil decomposition disabled, do not exceed the known FFT
        // plane capacity. No task-group assumption is made for ph.x.
        if (fftPlanes != null && ranksPerPool > fftPlanes) continue;
        const poolCpus = ranksPerPool * threads;
        const compatible = Boolean(input.phonon && validSource && ranksPerImage === source!.ranks && pools === source!.pools);
        const local = poolCpus <= cpusPerNode;
        const aligned = ranks % nodes === 0 && (local ? (ranks / nodes) % ranksPerPool === 0 : ranksPerPool % (ranks / nodes) === 0);
        // Benchmarks favor exploiting independent work while keeping FFT
        // communication local; none establishes a universal ranks/pool value.
        const score = Math.max(0, Math.log2(poolCpus / cpusPerNode)) * 2
          + Math.max(0, Math.log2(cpusPerNode / (2 * poolCpus))) * 0.7
          + (aligned ? 0 : 0.5)
          + (totalCpus % nodes === 0 && ranks % nodes !== 0 ? 2 : 0)
          + imbalance(kPoints, pools) * 0.8 + (input.phonon ? imbalance(irreps ?? qPoints, images) * 0.8 : 0)
          + Math.log2(pools) * 0.08 // modest replication cost; not a memory estimate
          + Math.log2(images) * 0.03 // image initialization and collection overhead
          + (standardThreads ? Math.abs(Math.log2(threads / preferredThreads)) * 0.12 : 4 + Math.abs(Math.log2(threads / 2)))
          + (fftPlanes == null ? 0 : Math.max(0, Math.log2(ranksPerPool / (fftPlanes / 4))) * 0.6)
          + (bands == null ? 0 : Math.max(0, Math.log2(ranksPerPool / bands)) * 0.3);
        candidates.push({ ranks, threads, images, pools, ranksPerImage, ranksPerPool, compatible, score });
      }
    }
  }
  if (!candidates.length) throw new Error("This allocation exceeds the known FFT plane capacity. Reduce CPUs or Nodes, or use Manual with a verified QE layout.");
  // Distributed restart compatibility outranks guesses about performance.
  const compatibleCandidates = candidates.filter((candidate) => candidate.compatible);
  const eligible = compatibleCandidates.length ? compatibleCandidates : candidates;
  const best = eligible.reduce((best, candidate) => candidate.score < best.score ? candidate : best);
  const { ranks, threads, images, pools, ranksPerImage, ranksPerPool } = best;
  notes.push("Compared MPI/OpenMP layouts using node locality, point balance and group replication. No universal ranks-per-pool optimum is assumed; benchmark representative work to tune further.");
  if (input.phonon) {
    notes.push(irreps == null
      ? "Actual phonon irrep count unavailable: images are conservatively limited to estimated q-points. Supply the selected irrep count in Advanced settings to allow more images."
      : `Images are limited to the supplied ${irreps} selected phonon irreps. Irrep costs can differ, so equal counts do not guarantee equal work.`);
  }
  if (best.compatible) notes.push("Preserving the source SCF rank and electronic-pool layout per image.");
  else if (input.phonon && validSource) notes.push(`Source SCF used ${source!.ranks} ranks and ${source!.pools} pools; distributed-format restart files need a matching layout. Use Manual to match the source.`);
  if (threads > 1) notes.push(`Hybrid MPI/OpenMP: ${ranks} MPI ranks × ${threads} threads per rank.`);
  if (![1, 2, 4].includes(threads)) notes.push("An unusual thread count is required to fit the CPU budget, restart layout or rank limit; test its efficiency before a large run.");
  if (input.maxRanks != null && rankLimit < totalCpus) notes.push(`Respecting the workflow limit of ${rankLimit} MPI ranks.`);
  if (fftPlanes != null) notes.push(`Using ${fftPlanes} FFT planes from the source SCF as guidance; the phonon FFT grid can differ.`);
  else notes.push("FFT dimensions unavailable; node locality is the fallback. Memory requirements are not predicted.");
  if (input.phonon && kPoints != null) notes.push("Source SCF k-point count guides pools; electronic symmetry and workload change with q.");
  if (ranksPerPool * threads > cpusPerNode || threads > 4) notes.push("This allocation may exceed useful parallelism for the available points. Auto preserves the requested CPUs, but a smaller allocation may be more efficient.");
  return { nodes, totalCpus, ranks, threads, tasksPerNode: ranks % nodes === 0 ? ranks / nodes : null,
    images, pools, ranksPerImage, ranksPerPool, notes };
}

export function applyQeAutoPlan(resources: SlurmResourceRequest, plan: QeAutoPlan): SlurmResourceRequest {
  return { ...resources, qe_parallel_mode: "auto", qe_total_cpus: plan.totalCpus,
    pencil_decomposition: false,
    nodes: plan.nodes, ntasks: plan.ranks, cpus_per_task: plan.threads, ntasks_per_node: plan.tasksPerNode,
    qe_qpoint_images: plan.images, qe_kpoint_pools: plan.pools };
}
