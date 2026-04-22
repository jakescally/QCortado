export interface EpwSourceRef {
  calc_id?: string;
  calc_type?: string;
}

export interface EpwArtifactManifestEntry {
  source_calc_id?: string;
  source_calc_type?: string;
  rel_path?: string;
  size_bytes?: number;
}

export interface EpwNamedArtifact {
  file_name?: string;
  size_bytes?: number;
}

export interface EpwResultSummary {
  completed?: boolean;
  elapsed_seconds?: number | null;
  core_metrics?: Record<string, number>;
  generated_outputs?: EpwNamedArtifact[];
  unknown_metrics?: Record<string, unknown>;
  parse_partial?: boolean;
  notes?: string[];
  warnings?: string[];
  parse_coverage?: string[];
}

export interface EpwErrorRecord {
  code?: string;
  message?: string;
  hint?: string | null;
}

export interface EpwParsedTable {
  file_name: string;
  family: string;
  title: string;
  column_labels: string[];
  rows: Array<Array<number | null>>;
  skipped?: boolean;
  skip_reason?: string | null;
}

export interface EpwMobilityDataset {
  carrier_type: string;
  method: string;
  iteration?: number | null;
  converged?: boolean | null;
  max_error?: number | null;
  component_labels: string[];
  temperature_values_k: number[];
  fermi_values_ev: number[];
  density_values_cm3: number[];
  population_values: Array<number | null>;
  mobility_values: Array<Array<number | null>>;
}

export interface EpwTransportData {
  mobility?: EpwMobilityDataset[];
  scattering_file_notices?: number;
  notes?: string[];
  warnings?: string[];
}

export interface EpwEliashbergIteration {
  temperature_k?: number | null;
  iteration: number;
  ethr: number;
  znorm: number;
  delta_mev: number;
}

export interface EpwGapSummary {
  temperature_k?: number | null;
  free_energy_mev?: number | null;
  gap_min_mev?: number | null;
  gap_max_mev?: number | null;
}

export interface EpwSuperconductivityData {
  lambda?: number | null;
  lambda_tr?: number | null;
  electron_phonon_coupling?: number | null;
  tc_mcmillan_k?: number | null;
  tc_allen_dynes_k?: number | null;
  tc_sisso_k?: number | null;
  w_log_mev?: number | null;
  bcs_gap_mev?: number | null;
  muc?: number | null;
  frequency_cutoff_ev?: number | null;
  frequency_points?: number | null;
  eliashberg_converged?: boolean | null;
  eliashberg_nsiter?: number | null;
  temperatures_k?: number[];
  eliashberg_iterations?: EpwEliashbergIteration[];
  gap_summaries?: EpwGapSummary[];
  spectral_tables?: EpwParsedTable[];
  notes?: string[];
  warnings?: string[];
}

export interface EpwSelfEnergyMode {
  mode_label: string;
  lambda: number;
  lambda_tr?: number | null;
  gamma_mev?: number | null;
  gamma_tr_mev?: number | null;
  omega_mev?: number | null;
}

export interface EpwSpectralData {
  self_energy_modes?: EpwSelfEnergyMode[];
  tables?: EpwParsedTable[];
  notes?: string[];
  warnings?: string[];
}

export interface EpwViewerData {
  schema_version?: number;
  sources?: {
    phonon?: EpwSourceRef | null;
    wannier?: EpwSourceRef | null;
    scf?: EpwSourceRef | null;
    manifests?: EpwArtifactManifestEntry[];
  };
  input?: {
    prefix?: string;
    outdir?: string;
    dvscf_dir?: string;
    wannier_dir?: string;
    k_mesh?: [number, number, number];
    q_mesh?: [number, number, number];
    coarse_k_mesh?: [number, number, number] | null;
    fine_k_mesh?: [number, number, number] | null;
    coarse_q_mesh?: [number, number, number] | null;
    fine_q_mesh?: [number, number, number] | null;
    epbwrite?: boolean;
    epbread?: boolean;
    epwwrite?: boolean;
    epwread?: boolean;
    wannierize?: boolean;
    fsthick_ev?: number | null;
    degaussw_ev?: number | null;
    nbndsub?: number | null;
  };
  goals?: {
    coupling?: boolean;
    phonon_linewidth_a2f?: boolean;
    electron_self_energy?: boolean;
    transport_mobility?: boolean;
    superconductivity?: boolean;
  };
  extensions?: {
    superconductivity?: unknown;
  };
  runtime?: {
    pools?: number | null;
    max_seconds?: number | null;
    artifact_sync_mode?: string | null;
  };
  artifacts?: EpwNamedArtifact[];
  result_summary?: EpwResultSummary;
  transport?: EpwTransportData | null;
  superconductivity?: EpwSuperconductivityData | null;
  spectral?: EpwSpectralData | null;
  parsed_tables?: EpwParsedTable[];
  errors?: EpwErrorRecord[];
}

export interface EpwViewerPayload {
  data: EpwViewerData;
  rawOutput?: string | null;
}

export type EpwViewerTab = "overview" | "transport" | "superconductivity" | "spectral" | "files";
export type EpwMobilityMetric = "mobility" | "fermi" | "density" | "population";

export interface EpwSeriesPoint {
  index: number;
  x: number;
  y: number | null;
  label?: string;
}

export interface EpwComponentOption {
  value: string;
  label: string;
  tooltip: string;
}

const COMPONENT_ORDER = ["avg", "xx", "yy", "zz", "xy", "xz", "yx", "yz", "zx", "zy"];

function isFiniteNumber(value: unknown): value is number {
  return typeof value === "number" && Number.isFinite(value);
}

export function formatEpwNumber(value: unknown, decimals = 4): string {
  if (value === null || value === undefined) return "N/A";
  if (typeof value !== "number" || !Number.isFinite(value)) return String(value);
  const abs = Math.abs(value);
  if (abs > 0 && (abs < 0.001 || abs >= 100000)) {
    return value.toExponential(Math.min(4, Math.max(1, decimals)));
  }
  if (Number.isInteger(value)) return String(value);
  return value.toFixed(decimals).replace(/0+$/, "").replace(/\.$/, "");
}

export function formatEpwBytes(bytes: number): string {
  if (!Number.isFinite(bytes) || bytes <= 0) return "0 B";
  const units = ["B", "KB", "MB", "GB", "TB"];
  let value = bytes;
  let unitIndex = 0;
  while (value >= 1024 && unitIndex < units.length - 1) {
    value /= 1024;
    unitIndex += 1;
  }
  return `${value.toFixed(value >= 10 || unitIndex === 0 ? 0 : 1)} ${units[unitIndex]}`;
}

export function hasEpwTransportData(data: EpwViewerData | null | undefined): boolean {
  return Array.isArray(data?.transport?.mobility) && data!.transport!.mobility!.length > 0;
}

export function hasEpwSuperconductivityData(data: EpwViewerData | null | undefined): boolean {
  const superconductivity = data?.superconductivity;
  if (!superconductivity) return false;
  return [
    superconductivity.lambda,
    superconductivity.lambda_tr,
    superconductivity.electron_phonon_coupling,
    superconductivity.tc_mcmillan_k,
    superconductivity.tc_allen_dynes_k,
    superconductivity.tc_sisso_k,
    superconductivity.w_log_mev,
    superconductivity.bcs_gap_mev,
  ].some(isFiniteNumber)
    || Boolean(superconductivity.eliashberg_iterations?.length)
    || Boolean(superconductivity.gap_summaries?.length)
    || Boolean(superconductivity.spectral_tables?.length);
}

export function getDefaultEpwTab(data: EpwViewerData | null | undefined): EpwViewerTab {
  if (hasEpwTransportData(data)) return "transport";
  if (hasEpwSuperconductivityData(data)) return "superconductivity";
  return "overview";
}

export function getEpwComponentOptions(labels: string[] | null | undefined): EpwComponentOption[] {
  const normalized = Array.from(new Set((labels || []).map((label) => String(label).trim().toLowerCase()).filter(Boolean)));
  const ordered = [
    ...COMPONENT_ORDER.filter((label) => normalized.includes(label)),
    ...normalized.filter((label) => !COMPONENT_ORDER.includes(label)).sort(),
  ];
  return ordered.map((value) => ({
    value,
    label: value === "avg" ? "Avg" : value,
    tooltip: value === "avg"
      ? "Average of the available diagonal tensor components."
      : `Tensor component ${value}. EPW may sort or rotate components depending on the calculation settings.`,
  }));
}

export function getEpwMobilityValue(dataset: EpwMobilityDataset, component: string, index: number): number | null {
  const labels = (dataset.component_labels || []).map((label) => String(label).trim().toLowerCase());
  if (component === "avg" && !labels.includes("avg")) {
    const diagonalValues = ["xx", "yy", "zz"]
      .map((label) => labels.indexOf(label))
      .filter((componentIndex) => componentIndex >= 0)
      .map((componentIndex) => dataset.mobility_values?.[componentIndex]?.[index])
      .filter(isFiniteNumber);
    if (diagonalValues.length === 0) return null;
    return diagonalValues.reduce((sum, value) => sum + value, 0) / diagonalValues.length;
  }
  const componentIndex = labels.indexOf(component);
  if (componentIndex < 0) return null;
  const value = dataset.mobility_values?.[componentIndex]?.[index];
  return isFiniteNumber(value) ? value : null;
}

export function buildEpwMobilitySeries(
  dataset: EpwMobilityDataset,
  metric: EpwMobilityMetric,
  component: string,
): EpwSeriesPoint[] {
  return (dataset.temperature_values_k || []).map((temperature, index) => {
    let y: number | null = null;
    if (metric === "mobility") {
      y = getEpwMobilityValue(dataset, component, index);
    } else if (metric === "fermi") {
      const value = dataset.fermi_values_ev?.[index];
      y = isFiniteNumber(value) ? value : null;
    } else if (metric === "density") {
      const value = dataset.density_values_cm3?.[index];
      y = isFiniteNumber(value) ? value : null;
    } else {
      const value = dataset.population_values?.[index];
      y = isFiniteNumber(value) ? value : null;
    }
    return {
      index,
      x: temperature,
      y,
      label: `${formatEpwNumber(temperature, 2)} K`,
    };
  });
}

export function buildEpwTableSeries(table: EpwParsedTable, xIndex = 0, yIndex = 1): EpwSeriesPoint[] {
  return (table.rows || []).map((row, index) => ({
    index,
    x: isFiniteNumber(row[xIndex]) ? Number(row[xIndex]) : index,
    y: isFiniteNumber(row[yIndex]) ? Number(row[yIndex]) : null,
  }));
}

export function collectEpwWarnings(data: EpwViewerData): string[] {
  const warnings = [
    ...(data.result_summary?.warnings || []),
    ...(data.transport?.warnings || []),
    ...(data.superconductivity?.warnings || []),
    ...(data.spectral?.warnings || []),
    ...(data.errors || []).map((error) => `[${error.code || "error"}] ${error.message || "No message"}${error.hint ? ` Hint: ${error.hint}` : ""}`),
  ];
  return Array.from(new Set(warnings.filter(Boolean)));
}

export function mergeEpwArtifacts(data: EpwViewerData): EpwNamedArtifact[] {
  const merged = new Map<string, EpwNamedArtifact>();
  const push = (entry: EpwNamedArtifact | null | undefined) => {
    const name = String(entry?.file_name || "").trim();
    if (!name) return;
    merged.set(name, { file_name: name, size_bytes: Number(entry?.size_bytes) || 0 });
  };
  for (const entry of data.result_summary?.generated_outputs || []) push(entry);
  for (const entry of data.artifacts || []) push(entry);
  return Array.from(merged.values()).sort((left, right) => String(left.file_name).localeCompare(String(right.file_name)));
}
