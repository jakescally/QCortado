export interface TransportArtifact {
  file_name: string;
  size_bytes: number;
}

export interface TransportDataset {
  file_name: string;
  component_labels: string[];
  values: Array<Array<Array<number | null>>>;
}

export interface TransportTdfData {
  file_name: string;
  energy_values_ev: number[];
  component_labels: string[];
  values: Array<Array<number | null>>;
}

export interface TransportResult {
  engine: string;
  seedname: string;
  source_wannier_calc_id: string;
  reference_fermi_energy_ev: number;
  mu_values_ev: number[];
  mu_offsets_ev: number[];
  temperature_values_k: number[];
  relaxation_time_fs: number;
  is_2d: boolean;
  boltz_2d_dir?: string | null;
  conductivity: TransportDataset;
  sigma_s: TransportDataset;
  seebeck: TransportDataset;
  kappa: TransportDataset;
  tdf?: TransportTdfData | null;
  notes: string[];
  warnings: string[];
  artifact_manifest: TransportArtifact[];
}

export const TRANSPORT_METRIC_KEYS = [
  "conductivity",
  "sigma_s",
  "seebeck",
  "kappa",
] as const;

export type TransportMetricKey = typeof TRANSPORT_METRIC_KEYS[number];
export type TransportDefinitionKey = TransportMetricKey | "tdf";
export type TransportAxisMode = "mu" | "temperature";
export type TransportTickDensity = "sparse" | "normal" | "dense";
export type TransportPrecisionMode = "auto" | "manual";

export interface TransportMetricDefinition {
  key: TransportDefinitionKey;
  symbol: string;
  shortLabel: string;
  longLabel: string;
  defaultUnit: string;
  tooltip: string;
  tauBehavior: "dependent" | "independent" | "n/a";
  accentColor: string;
}

export interface TransportAxisSettings {
  showAbsoluteMu: boolean;
  customYRange: [number, number] | null;
  keepYRangeAcrossSlices: boolean;
  manualYMinInput: string;
  manualYMaxInput: string;
  tickDensity: TransportTickDensity;
  precisionMode: TransportPrecisionMode;
  precision: number;
  xLabelOverride: string;
  yLabelOverride: string;
  unitSuffixOverride: string;
  showGrid: boolean;
  showZeroLine: boolean;
}

export interface TransportViewerState {
  metric: TransportMetricKey;
  component: string;
  axisMode: TransportAxisMode;
  selectedTemperatureIndex: number;
  selectedMuIndex: number;
  heatmapOpen: boolean;
  settingsOpen: boolean;
  runContextOpen: boolean;
  tdfOpen: boolean;
  tdfComponent: string;
  axisSettings: TransportAxisSettings;
}

export interface TransportComponentOption {
  value: string;
  label: string;
  tooltip: string;
}

export interface TransportSeriesPoint {
  index: number;
  x: number;
  y: number | null;
  muValueEv: number;
  muOffsetEv: number;
  temperatureK: number;
}

export interface TransportHeatmapCell {
  muIndex: number;
  temperatureIndex: number;
  muValueEv: number;
  muOffsetEv: number;
  temperatureK: number;
  value: number | null;
}

export interface TransportSliceSelection {
  selectedTemperatureIndex: number;
  selectedMuIndex: number;
}

const TRANSPORT_DEFINITIONS: Record<TransportDefinitionKey, TransportMetricDefinition> = {
  conductivity: {
    key: "conductivity",
    symbol: "σ",
    shortLabel: "Sigma",
    longLabel: "Electrical conductivity",
    defaultUnit: "S/m",
    tooltip:
      "Electrical conductivity tensor from the BoltzWann interpolation. This quantity scales with the relaxation time chosen for the run.",
    tauBehavior: "dependent",
    accentColor: "#0f766e",
  },
  sigma_s: {
    key: "sigma_s",
    symbol: "σS",
    shortLabel: "Sigma S",
    longLabel: "Conductivity-weighted Seebeck",
    defaultUnit: "A/(m*K)",
    tooltip:
      "Conductivity-weighted Seebeck tensor. This is the product σS used internally in transport formulas and it retains the relaxation-time dependence of σ.",
    tauBehavior: "dependent",
    accentColor: "#b45309",
  },
  seebeck: {
    key: "seebeck",
    symbol: "S",
    shortLabel: "Seebeck",
    longLabel: "Seebeck coefficient",
    defaultUnit: "uV/K",
    tooltip:
      "Seebeck coefficient tensor. Under the current constant-relaxation-time treatment it is displayed as tau-independent.",
    tauBehavior: "independent",
    accentColor: "#be185d",
  },
  kappa: {
    key: "kappa",
    symbol: "K",
    shortLabel: "K",
    longLabel: "BoltzWann K",
    defaultUnit: "W/(m*K)",
    tooltip:
      "Upstream BoltzWann K output. QCortado shows the quantity as written by BoltzWann and does not silently relabel it as electronic thermal conductivity.",
    tauBehavior: "dependent",
    accentColor: "#7c3aed",
  },
  tdf: {
    key: "tdf",
    symbol: "TDF",
    shortLabel: "TDF",
    longLabel: "Transport distribution function",
    defaultUnit: "BoltzWann units",
    tooltip:
      "Transport distribution function versus energy. Use it as a diagnostic view of the interpolated transport kernel rather than as a direct replacement for the mu-T transport grids.",
    tauBehavior: "n/a",
    accentColor: "#2563eb",
  },
};

const DIAGONAL_COMPONENTS = ["xx", "yy", "zz"] as const;

function isFiniteNumber(value: unknown): value is number {
  return Number.isFinite(value);
}

function normalizeComponentLabel(value: unknown): string {
  return String(value ?? "").trim().toLowerCase();
}

function stripTrailingZeros(value: string): string {
  if (!value.includes(".")) {
    return value;
  }
  return value.replace(/(?:\.0+|(\.\d*?[1-9])0+)$/, "$1");
}

function getAutoDecimals(value: number): number {
  const magnitude = Math.abs(value);
  if (magnitude >= 1_000) return 0;
  if (magnitude >= 100) return 1;
  if (magnitude >= 10) return 2;
  if (magnitude >= 1) return 3;
  if (magnitude >= 0.1) return 4;
  return 5;
}

function getPreferredTickCount(density: TransportTickDensity): number {
  switch (density) {
    case "sparse":
      return 4;
    case "dense":
      return 8;
    default:
      return 6;
  }
}

function niceStep(rawStep: number): number {
  if (!Number.isFinite(rawStep) || rawStep <= 0) {
    return 1;
  }
  const exponent = Math.floor(Math.log10(rawStep));
  const fraction = rawStep / 10 ** exponent;
  let niceFraction = 1;
  if (fraction <= 1) {
    niceFraction = 1;
  } else if (fraction <= 2) {
    niceFraction = 2;
  } else if (fraction <= 2.5) {
    niceFraction = 2.5;
  } else if (fraction <= 5) {
    niceFraction = 5;
  } else {
    niceFraction = 10;
  }
  return niceFraction * 10 ** exponent;
}

export function getTransportMetricDefinition(
  key: TransportDefinitionKey,
): TransportMetricDefinition {
  return TRANSPORT_DEFINITIONS[key];
}

export function getDefaultTransportMetric(
  result: Partial<Record<TransportMetricKey, TransportDataset | null | undefined>>,
): TransportMetricKey {
  return TRANSPORT_METRIC_KEYS.find((key) => result[key] != null) ?? "conductivity";
}

export function createDefaultTransportAxisSettings(): TransportAxisSettings {
  return {
    showAbsoluteMu: false,
    customYRange: null,
    keepYRangeAcrossSlices: false,
    manualYMinInput: "",
    manualYMaxInput: "",
    tickDensity: "normal",
    precisionMode: "auto",
    precision: 4,
    xLabelOverride: "",
    yLabelOverride: "",
    unitSuffixOverride: "",
    showGrid: true,
    showZeroLine: true,
  };
}

export function findComponentIndex(labels: string[], target: string): number {
  return labels.findIndex((label) => normalizeComponentLabel(label) === target);
}

export function getTransportComponentOptions(
  dataset: Pick<TransportDataset, "component_labels"> | Pick<TransportTdfData, "component_labels"> | null | undefined,
): TransportComponentOption[] {
  const labels = Array.isArray(dataset?.component_labels)
    ? dataset.component_labels.map((label) => normalizeComponentLabel(label)).filter(Boolean)
    : [];
  const options: TransportComponentOption[] = [];
  const diagonalCount = DIAGONAL_COMPONENTS.filter((label) => findComponentIndex(labels, label) >= 0).length;

  if (diagonalCount > 0) {
    options.push({
      value: "avg",
      label: "avg",
      tooltip:
        "Arithmetic average of the available diagonal tensor components. QCortado averages the diagonal entries that are present in the dataset.",
    });
  }

  for (const preferred of DIAGONAL_COMPONENTS) {
    if (findComponentIndex(labels, preferred) >= 0) {
      options.push({
        value: preferred,
        label: preferred,
        tooltip: `Tensor component ${preferred}.`,
      });
    }
  }

  for (const label of labels) {
    if (options.some((option) => option.value === label)) {
      continue;
    }
    options.push({
      value: label,
      label,
      tooltip: `Tensor component ${label}.`,
    });
  }

  if (options.length === 0) {
    options.push({
      value: "value",
      label: "value",
      tooltip: "Single-component dataset value.",
    });
  }

  return options;
}

export function getDatasetComponentValue(
  dataset: TransportDataset | null | undefined,
  component: string,
  temperatureIndex: number,
  muIndex: number,
): number | null {
  if (!dataset) {
    return null;
  }
  const labels = Array.isArray(dataset.component_labels)
    ? dataset.component_labels.map((label) => normalizeComponentLabel(label))
    : [];
  const values = Array.isArray(dataset.values) ? dataset.values : [];

  const readValue = (componentIndex: number): number | null => {
    const raw = values?.[componentIndex]?.[temperatureIndex]?.[muIndex];
    return isFiniteNumber(raw) ? raw : null;
  };

  if (component === "avg") {
    const diagonalIndices = DIAGONAL_COMPONENTS
      .map((label) => findComponentIndex(labels, label))
      .filter((index) => index >= 0);
    const diagonalValues = diagonalIndices
      .map(readValue)
      .filter((value): value is number => value != null);
    if (diagonalValues.length === 0) {
      return null;
    }
    return diagonalValues.reduce((sum, value) => sum + value, 0) / diagonalValues.length;
  }

  const index = findComponentIndex(labels, component);
  if (index < 0) {
    return null;
  }
  return readValue(index);
}

export function getTdfComponentValue(
  dataset: TransportTdfData | null | undefined,
  component: string,
  energyIndex: number,
): number | null {
  if (!dataset) {
    return null;
  }
  const labels = Array.isArray(dataset.component_labels)
    ? dataset.component_labels.map((label) => normalizeComponentLabel(label))
    : [];
  const values = Array.isArray(dataset.values) ? dataset.values : [];

  const readValue = (componentIndex: number): number | null => {
    const raw = values?.[componentIndex]?.[energyIndex];
    return isFiniteNumber(raw) ? raw : null;
  };

  if (component === "avg") {
    const diagonalIndices = DIAGONAL_COMPONENTS
      .map((label) => findComponentIndex(labels, label))
      .filter((index) => index >= 0);
    const diagonalValues = diagonalIndices
      .map(readValue)
      .filter((value): value is number => value != null);
    if (diagonalValues.length === 0) {
      return null;
    }
    return diagonalValues.reduce((sum, value) => sum + value, 0) / diagonalValues.length;
  }

  const index = findComponentIndex(labels, component);
  if (index < 0) {
    return null;
  }
  return readValue(index);
}

export function clampTransportIndex(index: number, length: number): number {
  if (length <= 0) {
    return 0;
  }
  return Math.min(Math.max(0, index), length - 1);
}

export function findClosestValueIndex(values: number[], target: number): number {
  if (values.length === 0) {
    return 0;
  }
  let closestIndex = 0;
  let smallestDelta = Number.POSITIVE_INFINITY;
  values.forEach((value, index) => {
    const delta = Math.abs(value - target);
    if (delta < smallestDelta) {
      smallestDelta = delta;
      closestIndex = index;
    }
  });
  return closestIndex;
}

export function getDefaultTransportTemperatureIndex(values: number[]): number {
  return findClosestValueIndex(values, 300);
}

export function getDefaultTransportMuIndex(values: number[]): number {
  return findClosestValueIndex(values, 0);
}

export function buildTransportSeries(args: {
  result: TransportResult;
  dataset: TransportDataset | null | undefined;
  component: string;
  axisMode: TransportAxisMode;
  selectedTemperatureIndex: number;
  selectedMuIndex: number;
  showAbsoluteMu: boolean;
}): TransportSeriesPoint[] {
  const {
    result,
    dataset,
    component,
    axisMode,
    selectedTemperatureIndex,
    selectedMuIndex,
    showAbsoluteMu,
  } = args;

  const safeTemperatureIndex = clampTransportIndex(
    selectedTemperatureIndex,
    result.temperature_values_k.length,
  );
  const safeMuIndex = clampTransportIndex(selectedMuIndex, result.mu_offsets_ev.length);

  if (axisMode === "mu") {
    return result.mu_offsets_ev.map((muOffsetEv, index) => ({
      index,
      x: showAbsoluteMu ? (result.mu_values_ev[index] ?? muOffsetEv) : muOffsetEv,
      y: getDatasetComponentValue(dataset, component, safeTemperatureIndex, index),
      muValueEv: result.mu_values_ev[index] ?? muOffsetEv,
      muOffsetEv,
      temperatureK: result.temperature_values_k[safeTemperatureIndex] ?? 0,
    }));
  }

  const muValueEv = result.mu_values_ev[safeMuIndex] ?? 0;
  const muOffsetEv = result.mu_offsets_ev[safeMuIndex] ?? 0;
  return result.temperature_values_k.map((temperatureK, index) => ({
    index,
    x: temperatureK,
    y: getDatasetComponentValue(dataset, component, index, safeMuIndex),
    muValueEv,
    muOffsetEv,
    temperatureK,
  }));
}

export function buildTransportHeatmapGrid(
  result: TransportResult,
  dataset: TransportDataset | null | undefined,
  component: string,
): TransportHeatmapCell[][] {
  return result.temperature_values_k.map((temperatureK, temperatureIndex) =>
    result.mu_offsets_ev.map((muOffsetEv, muIndex) => ({
      muIndex,
      temperatureIndex,
      muValueEv: result.mu_values_ev[muIndex] ?? muOffsetEv,
      muOffsetEv,
      temperatureK,
      value: getDatasetComponentValue(dataset, component, temperatureIndex, muIndex),
    })),
  );
}

export function getNextTransportSliceSelection(
  axisMode: TransportAxisMode,
  clickedTemperatureIndex: number,
  clickedMuIndex: number,
  currentSelection: TransportSliceSelection,
): TransportSliceSelection {
  if (axisMode === "mu") {
    return {
      selectedTemperatureIndex: clickedTemperatureIndex,
      selectedMuIndex: currentSelection.selectedMuIndex,
    };
  }
  return {
    selectedTemperatureIndex: currentSelection.selectedTemperatureIndex,
    selectedMuIndex: clickedMuIndex,
  };
}

export function getTransportDataExtent(values: Array<number | null | undefined>): [number, number] {
  const finiteValues = values.filter((value): value is number => isFiniteNumber(value));
  if (finiteValues.length === 0) {
    return [0, 1];
  }
  let min = Math.min(...finiteValues);
  let max = Math.max(...finiteValues);
  if (Math.abs(max - min) < 1e-12) {
    const pad = Math.max(Math.abs(min) * 0.2, 1);
    min -= pad;
    max += pad;
    return [min, max];
  }
  const span = max - min;
  const pad = span * 0.1;
  min -= pad;
  max += pad;
  return [min, max];
}

export function buildTransportTicks(
  min: number,
  max: number,
  density: TransportTickDensity,
): number[] {
  if (!Number.isFinite(min) || !Number.isFinite(max)) {
    return [0];
  }
  if (max < min) {
    return [min, max];
  }
  if (Math.abs(max - min) < 1e-12) {
    return [min];
  }

  const targetCount = getPreferredTickCount(density);
  const step = niceStep((max - min) / Math.max(1, targetCount - 1));
  const start = Math.ceil(min / step) * step;
  const ticks: number[] = [];
  for (let tick = start; tick <= max + step * 0.25; tick += step) {
    ticks.push(Number(stripTrailingZeros(tick.toFixed(12))));
    if (ticks.length > 32) {
      break;
    }
  }

  if (ticks.length === 0) {
    return [min, max];
  }

  return ticks;
}

export function formatTransportNumber(
  value: number | null | undefined,
  precisionMode: TransportPrecisionMode = "auto",
  precision = 4,
): string {
  if (!isFiniteNumber(value)) {
    return "N/A";
  }
  const absValue = Math.abs(value);
  const decimals = precisionMode === "manual" ? precision : getAutoDecimals(value);
  if ((absValue >= 1e4 || (absValue > 0 && absValue < 1e-3)) && precisionMode === "auto") {
    return value.toExponential(3);
  }
  return stripTrailingZeros(value.toFixed(decimals));
}

export function resolveTransportYAxisLabel(
  definition: TransportMetricDefinition,
  settings: TransportAxisSettings,
): string {
  const base = settings.yLabelOverride.trim() || definition.symbol;
  const unit = resolveTransportUnitSuffix(definition, settings);
  return unit ? `${base} (${unit})` : base;
}

export function resolveTransportXAxisLabel(
  axisMode: TransportAxisMode,
  settings: TransportAxisSettings,
): string {
  const override = settings.xLabelOverride.trim();
  if (override) {
    return override;
  }
  if (axisMode === "temperature") {
    return "Temperature (K)";
  }
  return settings.showAbsoluteMu ? "μ (eV)" : "Δμ (eV)";
}

export function resolveTransportUnitSuffix(
  definition: TransportMetricDefinition,
  settings: Pick<TransportAxisSettings, "unitSuffixOverride">,
): string {
  const override = settings.unitSuffixOverride.trim();
  return override || definition.defaultUnit;
}

export function formatTransportComponentLabel(component: string): string {
  return component.trim().toLowerCase() || "value";
}

export function getTransportTauBadge(definition: TransportMetricDefinition): string {
  switch (definition.tauBehavior) {
    case "independent":
      return "tau-free";
    case "dependent":
      return "tau";
    default:
      return "derived";
  }
}
