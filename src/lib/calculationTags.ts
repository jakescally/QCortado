import type { EngineId } from "./engines/types";
import { isPhononReadyScf } from "./engines/qe/phononReady";

export type CalculationTagType = "info" | "feature" | "special" | "geometry" | "engine" | "muted";

export interface CalculationTagBadge {
  label: string;
  type: CalculationTagType;
}

export interface TaggableCalculation {
  engine_id?: EngineId | null;
  calc_type: string;
  parameters?: any;
  result?: any;
  scf_summary?: any;
  tags?: string[] | null;
}

const HIDDEN_TAGS = new Set(["pinned"]);

const TAG_ALIASES: Record<string, CalculationTagBadge | null> = {
  pinned: null,
  "phonon-ready": { label: "Phonon-Ready", type: "special" },
  "structure-optimized": { label: "Optimized", type: "special" },
  geometry: { label: "Geometry", type: "geometry" },
  hpc: { label: "HPC", type: "feature" },
  downloaded: { label: "Downloaded", type: "feature" },
  failed: { label: "Failed", type: "feature" },
  "wannier-ready": { label: "Wannier-Ready", type: "special" },
  "wien2k-native": { label: "WIEN2k", type: "special" },
  wien2k: { label: "WIEN2k", type: "special" },
  "scalar-relativistic": { label: "Scalar relativistic", type: "special" },
  soc: { label: "SOC", type: "feature" },
  "non-collinear": { label: "Non-collinear", type: "feature" },
  magnetic: { label: "Magnetic", type: "feature" },
  "dft+u": { label: "DFT+U", type: "feature" },
  vdw: { label: "vdW", type: "feature" },
  proj: { label: "Proj", type: "feature" },
  orb: { label: "Proj", type: "feature" },
  "hubbard-lrt": { label: "Hubbard-LRT", type: "special" },
  converged: { label: "Converged", type: "special" },
  qe: { label: "QE", type: "engine" },
  w2k: { label: "W2k", type: "engine" },
  dos: { label: "DOS", type: "feature" },
  dispersion: { label: "Dispersion", type: "feature" },
  "limited results": { label: "Limited results", type: "feature" },
  "needs review": { label: "Needs Review", type: "feature" },
  warning: { label: "Warning", type: "feature" },
};

function tagKey(label: string): string {
  return label.trim().toLowerCase();
}

function pushRawTag(tags: string[], label: string | null | undefined): void {
  const trimmed = String(label ?? "").trim();
  if (!trimmed) return;
  const normalized = tagKey(trimmed);
  if (HIDDEN_TAGS.has(normalized)) return;
  const displayKey = tagKey(TAG_ALIASES[normalized]?.label ?? trimmed);
  if (!tags.some((tag) => tagKey(TAG_ALIASES[tagKey(tag)]?.label ?? tag) === displayKey)) {
    tags.push(trimmed);
  }
}

function pushBadge(tags: CalculationTagBadge[], label: string, type: CalculationTagType): void {
  const trimmed = label.trim();
  if (!trimmed) return;
  const normalized = tagKey(trimmed);
  if (!tags.some((tag) => tagKey(tag.label) === normalized)) {
    tags.push({ label: trimmed, type });
  }
}

function inferTagType(label: string): CalculationTagType {
  const normalized = tagKey(label);
  if (normalized === "qe" || normalized === "w2k") return "engine";
  if (
    /^\d+\s*[x×]\s*\d+\s*[x×]\s*\d+/.test(label)
    || /^\d+(\.\d+)?\s*k-pts$/i.test(label)
    || /^\d+(\.\d+)?\s*(pts|wf|bands|frmsf|bxsf|modes)$/i.test(label)
    || /^rkmax\b/i.test(label)
    || /^[fe]\s+\d/i.test(label)
    || /^[-+]?\d+(\.\d+)?(e[-+]?\d+)?$/i.test(label)
    || /(\bfine k\b|\bq\b|\bk\b|\bmu\b|μ|\bt\b|\bfs\b)$/i.test(label)
    || normalized === "relax"
    || normalized === "vc-relax"
  ) {
    return "info";
  }
  if (normalized === "geometry") return "geometry";
  if (normalized.includes("ready") || normalized === "optimized" || normalized === "wien2k") {
    return "special";
  }
  return "feature";
}

export function normalizeCalculationTagBadges(rawTags?: string[] | null): CalculationTagBadge[] {
  const badges: CalculationTagBadge[] = [];
  for (const rawTag of rawTags ?? []) {
    const trimmed = String(rawTag ?? "").trim();
    if (!trimmed) continue;
    const normalized = tagKey(trimmed);
    if (HIDDEN_TAGS.has(normalized)) continue;
    const alias = TAG_ALIASES[normalized];
    if (alias === null) continue;
    if (alias) {
      pushBadge(badges, alias.label, alias.type);
    } else {
      pushBadge(badges, trimmed, inferTagType(trimmed));
    }
  }
  return badges;
}

export function getCalcTagClass(tag: CalculationTagBadge): string {
  const normalizedLabel = tag.label.trim().toUpperCase();
  return `calc-tag calc-tag-${tag.type}`
    + `${normalizedLabel === "HPC" ? " calc-tag-hpc" : ""}`
    + `${normalizedLabel === "DOWNLOADED" ? " calc-tag-downloaded" : ""}`
    + `${normalizedLabel === "FAILED" ? " calc-tag-failed" : ""}`;
}

export function isHpcCalculation(calc: Pick<TaggableCalculation, "parameters" | "result">): boolean {
  const params = calc.parameters || {};
  const backend = String(params.execution_backend || "").trim().toLowerCase();
  if (backend === "hpc") return true;
  if (params.remote_job_id || params.remote_workdir || params.remote_project_path) return true;
  const rawOutput = typeof calc.result?.raw_output === "string" ? calc.result.raw_output : "";
  return rawOutput.includes("HPC_STAGE|") || rawOutput.includes("HPC_CMD|");
}

export function isHpcArtifactsDownloaded(
  calc: Pick<TaggableCalculation, "parameters" | "tags">,
  downloadedIds?: Set<string>,
  calcId?: string,
): boolean {
  if (calcId && downloadedIds?.has(calcId)) return true;
  if (Array.isArray(calc.tags) && calc.tags.some((tag) => tagKey(String(tag)) === "downloaded")) return true;
  const params = calc.parameters || {};
  if (params.artifacts_downloaded_full === true) return true;
  return String(params.artifact_sync_mode || "").trim().toLowerCase() === "full";
}

export function isFailedCalculation(calc: Pick<TaggableCalculation, "parameters" | "tags">): boolean {
  if (Array.isArray(calc.tags) && calc.tags.some((tag) => tagKey(String(tag)) === "failed")) return true;
  const params = calc.parameters || {};
  const status = String(params.run_status ?? params.status ?? "").trim().toLowerCase();
  return status === "failed";
}

function formatThreshold(value: unknown): string | null {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return null;
  return numeric < 0.001 ? numeric.toExponential(0) : numeric.toString();
}

function formatGrid(value: unknown, suffix = ""): string | null {
  if (!Array.isArray(value) || value.length !== 3) return null;
  const [a, b, c] = value.map((entry) => Number(entry));
  if (![a, b, c].every((entry) => Number.isFinite(entry) && entry > 0)) return null;
  return `${a}×${b}×${c}${suffix}`;
}

function isWannierReadyScf(calc: TaggableCalculation): boolean {
  if (calc.calc_type !== "scf" || !calc.result?.converged) return false;
  const params = calc.parameters || {};
  const nspin = Number(params.nspin ?? 1);
  if (nspin !== 1) return false;
  if (Boolean(params.noncolin) || Boolean(params.lspinorb) || Boolean(params.lda_plus_u)) return false;
  const vdw = String(params.vdw_corr || "").trim().toLowerCase();
  return !vdw || vdw === "none";
}

function getWien2kScfParameters(calc: TaggableCalculation): Record<string, any> {
  if (calc.engine_id !== "wien2k" || calc.calc_type !== "scf") return {};
  return {
    ...(calc.parameters?.initialization ?? {}),
    ...(calc.parameters?.run ?? {}),
    ...(calc.parameters ?? {}),
  };
}

function isWien2kSocCalculation(calc: TaggableCalculation): boolean {
  return calc.engine_id === "wien2k"
    && calc.calc_type === "scf"
    && Boolean(
      calc.parameters?.spin_orbit
      || calc.parameters?.spinOrbit
      || calc.parameters?.soc
      || calc.tags?.some((tag) => tagKey(String(tag)) === "soc"),
    );
}

export function buildCalculationTagList(
  calc: TaggableCalculation,
  options: { downloadedIds?: Set<string>; calcId?: string } = {},
): string[] {
  const tags: string[] = [];
  const params = calc.parameters || {};
  for (const tag of calc.tags ?? []) {
    pushRawTag(tags, tag);
  }

  if (isFailedCalculation(calc)) pushRawTag(tags, "failed");

  if ((calc.engine_id ?? "qe") === "wien2k") {
    pushRawTag(tags, "wien2k-native");
  } else {
    pushRawTag(tags, "QE");
  }

  if (calc.calc_type === "optimization" || calc.calc_type === "relax" || calc.calc_type === "vcrelax") {
    pushRawTag(tags, "geometry");
    pushRawTag(tags, params.optimization_mode === "relax" || calc.calc_type === "relax" ? "Relax" : "VC-Relax");
    pushRawTag(tags, formatGrid(params.kgrid));
    const forceConv = formatThreshold(params.forc_conv_thr);
    if (forceConv) pushRawTag(tags, `F ${forceConv}`);
    const energyConv = formatThreshold(params.etot_conv_thr);
    if (energyConv) pushRawTag(tags, `E ${energyConv}`);
  } else if (calc.calc_type === "scf") {
    if (params.structure_source?.type === "optimization") pushRawTag(tags, "geometry");
    if (isPhononReadyScf(params, calc.tags)) pushRawTag(tags, "phonon-ready");
    if (isWannierReadyScf(calc)) pushRawTag(tags, "wannier-ready");
    pushRawTag(tags, formatGrid(params.kgrid));
    pushRawTag(tags, formatThreshold(params.conv_thr));
    if (params.lspinorb || isWien2kSocCalculation(calc)) pushRawTag(tags, "SOC");
    if (params.nspin === 4) pushRawTag(tags, "Non-collinear");
    if (params.nspin === 2) pushRawTag(tags, "Magnetic");
    if (params.lda_plus_u) pushRawTag(tags, "DFT+U");
    if (params.vdw_corr && params.vdw_corr !== "none") pushRawTag(tags, "vdW");

    const wien2kParams = getWien2kScfParameters(calc);
    pushRawTag(tags, formatGrid(wien2kParams.k_mesh ?? wien2kParams.kMesh));
    const rkmax = Number(wien2kParams.rkmax ?? wien2kParams.rkMax);
    if (Number.isFinite(rkmax) && rkmax > 0) {
      pushRawTag(tags, `RKMax ${rkmax >= 10 || Number.isInteger(rkmax) ? rkmax.toFixed(0) : rkmax.toFixed(1)}`);
    }
    const spinMode = String(wien2kParams.spin_mode ?? wien2kParams.spinMode ?? "").trim().toLowerCase();
    if (spinMode === "spin_polarized") pushRawTag(tags, "Spin-polarized");
  } else if (calc.calc_type === "bands") {
    if (params.total_k_points) pushRawTag(tags, `${params.total_k_points} k-pts`);
    if (params.lspinorb) pushRawTag(tags, "SOC");
    if (params.nspin === 2) pushRawTag(tags, "Magnetic");
    if (params.lda_plus_u) pushRawTag(tags, "DFT+U");
    if (params.vdw_corr && params.vdw_corr !== "none") pushRawTag(tags, "vdW");
    if (params.fat_bands_requested) pushRawTag(tags, "Proj");
  } else if (calc.calc_type === "dos") {
    pushRawTag(tags, formatGrid(params.dos_k_grid, " K"));
    if (params.n_points) pushRawTag(tags, `${params.n_points} pts`);
    if (params.lspinorb) pushRawTag(tags, "SOC");
    if (params.nspin === 2) pushRawTag(tags, "Magnetic");
    if (params.lda_plus_u) pushRawTag(tags, "DFT+U");
    if (params.vdw_corr && params.vdw_corr !== "none") pushRawTag(tags, "vdW");
  } else if (calc.calc_type === "wannier") {
    pushRawTag(tags, formatGrid(params.k_grid, " K"));
    if (params.n_wann) pushRawTag(tags, `${params.n_wann} WF`);
    if (params.n_bands) pushRawTag(tags, `${params.n_bands} bands`);
    if (Array.isArray(params.exclude_bands) && params.exclude_bands.length > 0) pushRawTag(tags, "Disentangled");
    if (params.total_spread != null && Number.isFinite(Number(params.total_spread))) {
      pushRawTag(tags, `Omega ${Number(params.total_spread).toFixed(3)}`);
    }
  } else if (calc.calc_type === "transport") {
    pushRawTag(tags, formatGrid(params.boltz_kmesh, " K"));
    if (params.mu_points) pushRawTag(tags, `${params.mu_points} μ`);
    if (params.temperature_points) pushRawTag(tags, `${params.temperature_points} T`);
    if (Number.isFinite(Number(params.relaxation_time_fs))) {
      pushRawTag(tags, `tau ${Number(params.relaxation_time_fs).toFixed(1)} fs`);
    }
    if (params.engine) pushRawTag(tags, String(params.engine));
    if (params.is_2d) pushRawTag(tags, "2D");
  } else if (calc.calc_type === "fermi_surface") {
    pushRawTag(tags, formatGrid(params.fermi_k_grid, " K"));
    if (params.n_frmsf_files) pushRawTag(tags, `${params.n_frmsf_files} FRMSF`);
    if (params.n_bxsf_files) pushRawTag(tags, `${params.n_bxsf_files} BXSF`);
    if (params.lspinorb) pushRawTag(tags, "SOC");
    if (params.nspin === 2) pushRawTag(tags, "Magnetic");
    if (params.lda_plus_u) pushRawTag(tags, "DFT+U");
    if (params.vdw_corr && params.vdw_corr !== "none") pushRawTag(tags, "vdW");
  } else if (calc.calc_type === "phonon") {
    pushRawTag(tags, formatGrid(params.q_grid, " Q"));
    if (params.n_modes) pushRawTag(tags, `${params.n_modes} modes`);
    if (params.calculate_dos) pushRawTag(tags, "DOS");
    if (params.calculate_dispersion) pushRawTag(tags, "Dispersion");
  } else if (calc.calc_type === "epw") {
    pushRawTag(tags, formatGrid(Array.isArray(params.fine_k_grid) ? params.fine_k_grid : params.k_mesh, " fine K"));
    pushRawTag(tags, formatGrid(Array.isArray(params.coarse_q_grid) ? params.coarse_q_grid : params.q_mesh, " Q"));
    if (params.parse_partial === true) pushRawTag(tags, "Limited results");
  }

  if (isHpcCalculation(calc)) {
    pushRawTag(tags, "HPC");
    if (isHpcArtifactsDownloaded(calc, options.downloadedIds, options.calcId)) {
      pushRawTag(tags, "Downloaded");
    }
  }

  return tags;
}

export function getCalculationTagBadges(
  calc: TaggableCalculation,
  options: { legacyFallback?: boolean; downloadedIds?: Set<string>; calcId?: string } = {},
): CalculationTagBadge[] {
  const stored = normalizeCalculationTagBadges(calc.tags);
  if (!options.legacyFallback) return stored;
  const derived = normalizeCalculationTagBadges(
    buildCalculationTagList(calc, { downloadedIds: options.downloadedIds, calcId: options.calcId }),
  );
  return derived.length > stored.length ? derived : stored;
}
