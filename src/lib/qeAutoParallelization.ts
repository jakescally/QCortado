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
  // Prefer modest hybrid groups on large allocations. Penalize oversized FFT
  // groups more strongly than deviation from this preferred thread count.
  const preferredThreads = totalCpus >= 128 ? (cpusPerNode >= 32 ? 4 : cpusPerNode >= 16 ? 2 : 1) : 1;
  const imbalance = (points: number | null, count: number) => points == null ? 0 : 1 - points / (Math.ceil(points / count) * count);
  const candidates = rankCandidates.map((ranks) => {
    const threads = totalCpus / ranks;
    let images = 1;
    let pools = 1;
    const compatible = Boolean(input.phonon && validSource && ranks % source!.ranks === 0
      && ranks / source!.ranks <= (qPoints ?? 1));
    if (compatible) {
      images = ranks / source!.ranks;
      pools = source!.pools;
    } else {
      if (input.phonon && qPoints != null) {
        const goal = Math.min(qPoints, Math.max(1, Math.ceil(ranks / 8)));
        const imageCounts = divisors(ranks).filter((count) => count <= qPoints);
        const score = (count: number) => Math.abs(Math.log2(count / goal)) + imbalance(qPoints, count) * 0.5
          + (count % nodes === 0 || nodes % count === 0 ? 0 : 0.1);
        images = imageCounts.reduce((best, count) => score(count) < score(best) ? count : best, 1);
      }
      const perImage = ranks / images;
      const target = Math.min(8, Math.ceil(ranks / nodes));
      if (kPoints != null && perImage > target) {
        const poolCounts = divisors(perImage).filter((count) => count <= kPoints);
        const score = (count: number) => Math.abs(Math.log2((perImage / count) / target)) + imbalance(kPoints, count) * 0.5;
        pools = poolCounts.reduce((best, count) => score(count) < score(best) ? count : best, 1);
      }
    }
    const ranksPerImage = ranks / images;
    const ranksPerPool = ranksPerImage / pools;
    const score = Math.max(0, Math.log2(ranksPerPool / 8)) * 3
      + Math.abs(Math.log2(threads / preferredThreads)) * 0.8
      + imbalance(kPoints, pools) * 0.5 + (input.phonon ? imbalance(qPoints, images) * 0.5 : 0)
      + (totalCpus % nodes === 0 && ranks % nodes !== 0 ? 1 : 0)
      + (threads > 8 ? Math.log2(threads / 8) * 0.5 : 0);
    return { ranks, threads, images, pools, ranksPerImage, ranksPerPool, compatible, score };
  });
  // A saved distributed restart's MPI layout takes priority over performance
  // guesses. Hybrid threads may change without changing that MPI layout.
  const compatibleCandidates = candidates.filter((candidate) => candidate.compatible);
  const eligible = compatibleCandidates.length ? compatibleCandidates : candidates;
  const best = eligible.reduce((best, candidate) => candidate.score < best.score ? candidate : best);
  const { ranks, threads, images, pools, ranksPerImage, ranksPerPool } = best;
  if (best.compatible) notes.push("Preserving the source SCF rank and electronic-pool layout per image.");
  else if (input.phonon && validSource) {
    notes.push(`Source SCF used ${source!.ranks} ranks and ${source!.pools} pools; distributed-format restart files need a matching layout. Use Manual to match the source.`);
  }
  if (threads > 1) notes.push(`Hybrid MPI/OpenMP: ${ranks} MPI ranks × ${threads} threads per rank; threads share each rank's allocated CPUs within a node.`);
  if (input.maxRanks != null && rankLimit < totalCpus) notes.push(`Respecting the workflow limit of ${rankLimit} MPI ranks.`);
  if (ranksPerPool > 8) notes.push(`${ranksPerPool} ranks share each electronic pool; consider a smaller allocation if scaling is poor.`);
  if (threads > 8 || ranksPerPool * threads > 64) {
    notes.push("This allocation may exceed useful parallelism for the available points. Auto preserves the requested CPUs, but a smaller allocation may be faster; compare timings.");
  }
  return { nodes, totalCpus, ranks, threads, tasksPerNode: ranks % nodes === 0 ? ranks / nodes : null,
    images, pools, ranksPerImage, ranksPerPool, notes };
}

export function applyQeAutoPlan(resources: SlurmResourceRequest, plan: QeAutoPlan): SlurmResourceRequest {
  return { ...resources, qe_parallel_mode: "auto", qe_total_cpus: plan.totalCpus,
    pencil_decomposition: false,
    nodes: plan.nodes, ntasks: plan.ranks, cpus_per_task: plan.threads, ntasks_per_node: plan.tasksPerNode,
    qe_qpoint_images: plan.images, qe_kpoint_pools: plan.pools };
}
