import { useState, useRef, useMemo, useCallback, useEffect, useId, useLayoutEffect } from "react";
import { useTheme } from "../lib/ThemeContext";
import type {
  BandDataset,
  BandProjectionDataset as NormalizedBandProjectionDataset,
  BandProjectionGroup as NormalizedBandProjectionGroup,
} from "../lib/viewers/bands/types";

interface HighSymmetryMarker {
  k_distance: number;
  label: string;
}

interface BandGap {
  value: number;
  is_direct: boolean;
  vbm_k: number;
  cbm_k: number;
  vbm_energy: number;
  cbm_energy: number;
}

export interface BandProjectionGroup {
  id: string;
  label: string;
  kind: "atom" | "orbital" | string;
  weights: number[][];
}

export interface BandProjectionData {
  source: string;
  atom_groups: BandProjectionGroup[];
  orbital_groups: BandProjectionGroup[];
  element_orbital_groups?: BandProjectionGroup[];
}

export interface BandData {
  k_points: number[];
  energies: number[][];
  fermi_energy: number;
  high_symmetry_points: HighSymmetryMarker[];
  n_bands: number;
  n_kpoints: number;
  band_gap: BandGap | null;
  energy_range: [number, number];
  projections?: BandProjectionData | null;
}

export type BandPlotData = BandData | BandDataset;

interface BandPlotProps {
  data: BandPlotData;
  width?: number;
  height?: number;
  energyRange?: [number, number];
  showFermiLevel?: boolean;
  /** Actual Fermi energy from SCF calculation (bands always reports 0) */
  scfFermiEnergy?: number;
  yAxisLabel?: string;
  pointLabel?: string;
  valueLabel?: string;
  valueUnit?: string;
  valueDecimals?: number;
  primaryCountLabel?: string;
  secondaryCountLabel?: string;
  scrollHint?: string;
  yClampRange?: [number, number] | null;
  viewerType?: "electronic" | "phonon";
  sharedSettings?: BandPlotSharedSettings | null;
  showSidebar?: boolean;
  projectionSelection?: string | null;
  enableWheelRangeControl?: boolean;
  enableHoverScrollLock?: boolean;
  comparisonOptions?: BandPlotComparisonOption[];
  comparisonTitle?: string;
  comparisonNoneLabel?: string;
  windowOverlays?: BandPlotWindowOverlay[];
  onWindowOverlayChange?: (update: BandPlotWindowOverlayUpdate) => void;
  windowOverlayHint?: string;
  showBandGapOverlayOverride?: boolean;
  calculationParameters?: Record<string, unknown> | null;
  onPersistSelectedValenceBandIndex?: (bandIndex: number | null) => Promise<void> | void;
}

interface HoveredPoint {
  band: number;
  k: number;
  energy: number;
  x: number;
  y: number;
  projectionWeight?: number;
  projectionWeightNormalized?: number;
}

interface OverlayDragState {
  overlayId: string;
  pointerId: number;
  mode: BandPlotWindowOverlayInteraction;
  startClientY: number;
  startMin: number;
  startMax: number;
}

export type ColorMode = "single" | "rainbow";
export type RainbowPalette = "jet" | "sinebow";
type ProjectionMode = "atom" | "orbital";
type ProjectionNormalizeMode = "global" | "band";
type ProjectionRenderMode = "fat" | "color";
type FatColorMode = "accent" | "band";
export type FermiReferenceMode = "scf" | "bands";
export type ZeroEnergyReferenceMode = "calculated-fermi" | "vbm";

export interface BandPlotSharedSettings {
  fermiReferenceMode: FermiReferenceMode;
  zeroEnergyReferenceMode: ZeroEnergyReferenceMode;
  lineWidth: number;
  lineOpacity: number;
  plotTextScale: number;
  colorMode: ColorMode;
  singleBandColor: string;
  rainbowPalette: RainbowPalette;
  plotBgWhite: boolean;
  showBandGapOverlay: boolean;
}

export interface BandPlotProjectionOption {
  value: string;
  label: string;
}

export interface BandPlotComparisonOption {
  id: string;
  label: string;
  data: BandPlotData;
}

interface NormalizedBandPlotComparisonOption extends Omit<BandPlotComparisonOption, "data"> {
  data: BandData;
}

function isBandDataset(data: BandPlotData): data is BandDataset {
  return Boolean(
    data &&
    typeof data === "object" &&
    (data as BandDataset).schema === "cortado.band_path.v1" &&
    Array.isArray((data as BandDataset).series),
  );
}

function finiteNumberOrNull(value: unknown): number | null {
  const parsed = Number(value);
  return Number.isFinite(parsed) ? parsed : null;
}

function finiteNumberOrFallback(value: unknown, fallback: number): number {
  return finiteNumberOrNull(value) ?? fallback;
}

function numberArrayFromMetadata(value: unknown): [number, number] | null {
  if (!Array.isArray(value) || value.length < 2) {
    return null;
  }
  const min = finiteNumberOrNull(value[0]);
  const max = finiteNumberOrNull(value[1]);
  return min != null && max != null ? [min, max] : null;
}

function calculateEnergyRangeFromBands(energies: number[][]): [number, number] {
  let min = Number.POSITIVE_INFINITY;
  let max = Number.NEGATIVE_INFINITY;

  for (const band of energies) {
    for (const value of band) {
      if (!Number.isFinite(value)) {
        continue;
      }
      min = Math.min(min, value);
      max = Math.max(max, value);
    }
  }

  return Number.isFinite(min) && Number.isFinite(max) ? [min, max] : [0, 0];
}

function mapBandDatasetProjections(
  projections: NormalizedBandProjectionDataset | null | undefined,
): BandProjectionData | null {
  if (!projections) {
    return null;
  }

  const toLegacyGroup = (group: NormalizedBandProjectionGroup): BandProjectionGroup => ({
    id: group.id,
    label: group.label,
    kind: group.kind,
    weights: group.weights,
  });

  const atomGroups: BandProjectionGroup[] = [];
  const orbitalGroups: BandProjectionGroup[] = [];
  const elementOrbitalGroups: BandProjectionGroup[] = [];

  for (const group of projections.groups) {
    if (group.kind === "atom") {
      atomGroups.push(toLegacyGroup(group));
    } else if (group.kind === "element_orbital") {
      elementOrbitalGroups.push(toLegacyGroup(group));
    } else {
      orbitalGroups.push(toLegacyGroup(group));
    }
  }

  return {
    source: projections.source,
    atom_groups: atomGroups,
    orbital_groups: orbitalGroups,
    element_orbital_groups: elementOrbitalGroups,
  };
}

export function bandDatasetToBandData(dataset: BandDataset): BandData {
  const energies = dataset.series.map((series) =>
    series.values.map((value) => (typeof value === "number" && Number.isFinite(value) ? value : Number.NaN)),
  );
  const metadata = dataset.metadata ?? {};
  const energyRange = numberArrayFromMetadata(metadata.energyRangeEv)
    ?? calculateEnergyRangeFromBands(energies);
  const bandGap = dataset.bandGap
    ? {
      value: dataset.bandGap.valueEv,
      is_direct: dataset.bandGap.isDirect,
      vbm_k: finiteNumberOrFallback(dataset.bandGap.vbmX, 0),
      cbm_k: finiteNumberOrFallback(dataset.bandGap.cbmX, 0),
      vbm_energy: finiteNumberOrFallback(dataset.bandGap.vbmEnergyEv, 0),
      cbm_energy: finiteNumberOrFallback(dataset.bandGap.cbmEnergyEv, 0),
    }
    : null;

  return {
    k_points: dataset.x,
    energies,
    fermi_energy: finiteNumberOrFallback(
      dataset.referenceEnergyEv,
      finiteNumberOrFallback(metadata.legacyFermiEnergyEv, 0),
    ),
    high_symmetry_points: dataset.markers.map((marker) => ({
      k_distance: marker.x,
      label: marker.label,
    })),
    n_bands: finiteNumberOrFallback(metadata.nBands, dataset.series.length),
    n_kpoints: finiteNumberOrFallback(metadata.nKpoints, dataset.x.length),
    band_gap: bandGap,
    energy_range: energyRange,
    projections: mapBandDatasetProjections(dataset.projections),
  };
}

export function normalizeBandPlotData(data: BandPlotData): BandData {
  return isBandDataset(data) ? bandDatasetToBandData(data) : data;
}

export type BandPlotWindowOverlaySide = "left" | "right";
export type BandPlotWindowOverlayInteraction =
  | "move"
  | "resize-min"
  | "resize-max"
  | "slider-min"
  | "slider-max";

export interface BandPlotWindowOverlay {
  id: string;
  min: number;
  max: number;
  color: string;
  side: BandPlotWindowOverlaySide;
  label: string;
}

export interface BandPlotWindowOverlayUpdate {
  id: string;
  min: number;
  max: number;
  interaction: BandPlotWindowOverlayInteraction;
}

interface OrbitalElementOption {
  key: string;
  label: string;
  groups: BandProjectionGroup[];
}

interface ProjectionSelectionEntry extends BandPlotProjectionOption {
  mode: ProjectionMode;
  elementKey?: string;
}

const BAND_GAP_TOLERANCE_EV = 0.01;
const DIRECT_GAP_K_TOLERANCE = 0.01;

// Format high-symmetry point labels (handle Greek letters)
function formatLabel(label: string): string {
  const greekMap: Record<string, string> = {
    G: "Γ",
    Gamma: "Γ",
    GAMMA: "Γ",
    "Σ": "Σ",
    Sigma: "Σ",
    Delta: "Δ",
    Lambda: "Λ",
  };
  return greekMap[label] || label;
}

function clamp01(value: number): number {
  if (value <= 0) return 0;
  if (value >= 1) return 1;
  return value;
}

function clampToRange(value: number, min: number, max: number): number {
  return Math.min(Math.max(value, min), max);
}

function distanceSquaredToSegment(
  px: number,
  py: number,
  x1: number,
  y1: number,
  x2: number,
  y2: number,
): number {
  const dx = x2 - x1;
  const dy = y2 - y1;
  const lengthSquared = dx * dx + dy * dy;

  if (lengthSquared <= 1e-12) {
    const pointDx = px - x1;
    const pointDy = py - y1;
    return pointDx * pointDx + pointDy * pointDy;
  }

  const t = clampToRange(((px - x1) * dx + (py - y1) * dy) / lengthSquared, 0, 1);
  const projX = x1 + t * dx;
  const projY = y1 + t * dy;
  const distX = px - projX;
  const distY = py - projY;
  return distX * distX + distY * distY;
}

const MIN_WINDOW_SPAN_EV = 1e-4;

function rgbString(r: number, g: number, b: number): string {
  const rc = Math.round(Math.max(0, Math.min(255, r)));
  const gc = Math.round(Math.max(0, Math.min(255, g)));
  const bc = Math.round(Math.max(0, Math.min(255, b)));
  return `rgb(${rc}, ${gc}, ${bc})`;
}

// A compact "jet-like" map to mimic many legacy fat-band plots.
function jetColor(t: number): string {
  const x = clamp01(t);
  const r = 255 * clamp01(1.5 - Math.abs(4 * x - 3));
  const g = 255 * clamp01(1.5 - Math.abs(4 * x - 2));
  const b = 255 * clamp01(1.5 - Math.abs(4 * x - 1));
  return rgbString(r, g, b);
}

function sinebowColor(t: number): string {
  const x = clamp01(t);
  const a = Math.PI * 2 * (0.5 - x);
  const r = 255 * Math.pow(Math.sin(a), 2);
  const g = 255 * Math.pow(Math.sin(a + (2 * Math.PI) / 3), 2);
  const b = 255 * Math.pow(Math.sin(a + (4 * Math.PI) / 3), 2);
  return rgbString(r, g, b);
}

function projectionScaleColor(t: number): string {
  const x = clamp01(t);
  const blue = { r: 25, g: 70, b: 238 };
  const red = { r: 232, g: 48, b: 38 };
  return rgbString(
    blue.r + (red.r - blue.r) * x,
    blue.g + (red.g - blue.g) * x,
    blue.b + (red.b - blue.b) * x,
  );
}

function projectionWeightColor(value: number, min: number, max: number): string {
  const span = max - min;
  const t = span > 1e-12 ? (value - min) / span : 0;
  return projectionScaleColor(t);
}

function buildProjectionLegendTicks(
  min: number,
  max: number,
  tickCount: number,
): number[] {
  const count = Math.max(2, Math.min(12, Math.round(tickCount)));
  if (!Number.isFinite(min) || !Number.isFinite(max) || max <= min) {
    return [];
  }
  return Array.from({ length: count }, (_, index) => {
    const t = count <= 1 ? 0 : index / (count - 1);
    return max - (max - min) * t;
  });
}

function formatProjectionLegendLabel(label: string): string {
  const trimmed = label.trim() || "Selected";
  if (/projection$/i.test(trimmed)) {
    return trimmed;
  }
  return `${trimmed} orbital projection`;
}

function bandColorForIndex(
  bandIndex: number,
  totalBands: number,
  colorMode: ColorMode,
  singleColor: string,
  rainbowPalette: RainbowPalette,
): string {
  if (colorMode === "single" || totalBands <= 1) {
    return singleColor;
  }

  const t = totalBands <= 1 ? 0.5 : bandIndex / (totalBands - 1);
  return rainbowPalette === "jet" ? jetColor(t) : sinebowColor(t);
}

function getYAxisTickStep(range: number): number {
  if (range > 500) return 100;
  if (range > 200) return 50;
  if (range > 100) return 25;
  if (range > 20) return 5;
  if (range > 10) return 2;
  if (range > 5) return 1;
  return 0.5;
}

function getYAxisTicks(eMin: number, eMax: number, customStep?: number | null): number[] {
  if (!Number.isFinite(eMin) || !Number.isFinite(eMax) || eMax < eMin) {
    return [0];
  }
  const step = customStep != null && Number.isFinite(customStep) && customStep > 0
    ? customStep
    : getYAxisTickStep(eMax - eMin);
  const ticks: number[] = [];
  let tick = Math.ceil(eMin / step) * step;
  while (tick <= eMax + 1e-9) {
    ticks.push(tick);
    tick += step;
  }
  return ticks;
}

function formatAxisInputValue(value: number): string {
  if (!Number.isFinite(value)) return "";
  return Number.parseFloat(value.toFixed(6)).toString();
}

function getTickDecimals(step: number | null): number {
  if (step == null || !Number.isFinite(step) || step <= 0) {
    return 1;
  }
  const normalized = formatAxisInputValue(step);
  const decimalPart = normalized.split(".")[1] ?? "";
  return Math.min(6, decimalPart.length);
}

function getElectronicReferenceLabelSubscript(label: string): string | null {
  const normalized = label.replace(/\s+/g, "").replace(/−/g, "-").toLowerCase();
  if (normalized === "e-e_f(ev)") {
    return "F";
  }
  if (normalized === "e-e_vbm(ev)") {
    return "VBM";
  }
  return null;
}

function normalizeSymmetryLabel(label: string): string {
  const trimmed = String(label ?? "").trim();
  if (!trimmed) return "";
  return formatLabel(trimmed).replace(/\s+/g, "").toUpperCase();
}

function remapComparisonKPoints(
  referenceData: BandData,
  comparisonData: BandData,
): number[] | null {
  const comparisonKPoints = comparisonData.k_points ?? [];
  if (comparisonKPoints.length === 0) {
    return null;
  }

  const referenceMarkers = referenceData.high_symmetry_points ?? [];
  const comparisonMarkers = comparisonData.high_symmetry_points ?? [];
  if (referenceMarkers.length >= 2 && comparisonMarkers.length >= 2) {
    const sameMarkerCount = referenceMarkers.length === comparisonMarkers.length;
    const labelsMatch = sameMarkerCount && referenceMarkers.every((marker, index) => (
      normalizeSymmetryLabel(marker.label) === normalizeSymmetryLabel(comparisonMarkers[index]?.label)
    ));

    if (labelsMatch) {
      const mapped: number[] = [];
      let segmentIndex = 0;
      for (const kPoint of comparisonKPoints) {
        while (
          segmentIndex < comparisonMarkers.length - 2 &&
          kPoint > comparisonMarkers[segmentIndex + 1].k_distance + 1e-9
        ) {
          segmentIndex += 1;
        }

        const leftComparison = comparisonMarkers[segmentIndex].k_distance;
        const rightComparison = comparisonMarkers[segmentIndex + 1].k_distance;
        const leftReference = referenceMarkers[segmentIndex].k_distance;
        const rightReference = referenceMarkers[segmentIndex + 1].k_distance;
        const span = rightComparison - leftComparison;
        const t = span > 1e-12 ? clamp01((kPoint - leftComparison) / span) : 0;
        mapped.push(leftReference + t * (rightReference - leftReference));
      }
      return mapped;
    }
  }

  const referenceStart = referenceData.k_points[0];
  const referenceEnd = referenceData.k_points[referenceData.k_points.length - 1];
  const comparisonStart = comparisonKPoints[0];
  const comparisonEnd = comparisonKPoints[comparisonKPoints.length - 1];
  const comparisonSpan = comparisonEnd - comparisonStart;
  const referenceSpan = referenceEnd - referenceStart;

  if (
    Number.isFinite(referenceStart) &&
    Number.isFinite(referenceEnd) &&
    Number.isFinite(comparisonStart) &&
    Number.isFinite(comparisonEnd) &&
    Math.abs(comparisonSpan) > 1e-12
  ) {
    return comparisonKPoints.map((kPoint) => (
      referenceStart + ((kPoint - comparisonStart) / comparisonSpan) * referenceSpan
    ));
  }

  return null;
}

export function resolveBandPlotFermiContext(
  inputData: BandPlotData,
  scfFermiEnergy?: number | null,
  requestedMode: FermiReferenceMode = "bands",
): {
  fermiEnergy: number;
  mode: FermiReferenceMode;
  hasScfFermi: boolean;
  hasBandsFermi: boolean;
} {
  const data = normalizeBandPlotData(inputData);
  const hasScfFermi = scfFermiEnergy != null && Number.isFinite(scfFermiEnergy);
  const hasBandsFermi = Number.isFinite(data.fermi_energy);

  let mode: FermiReferenceMode;
  if (requestedMode === "scf" && hasScfFermi) {
    mode = "scf";
  } else if (requestedMode === "bands" && hasBandsFermi) {
    mode = "bands";
  } else if (hasBandsFermi) {
    mode = "bands";
  } else if (hasScfFermi) {
    mode = "scf";
  } else {
    mode = "bands";
  }

  const fermiEnergy = mode === "scf" && hasScfFermi
    ? (scfFermiEnergy as number)
    : hasBandsFermi
      ? data.fermi_energy
      : hasScfFermi
        ? (scfFermiEnergy as number)
        : 0;

  return {
    fermiEnergy,
    mode,
    hasScfFermi,
    hasBandsFermi,
  };
}

export function getDefaultBandPlotEnergyRange(
  inputData: BandPlotData,
  scfFermiEnergy?: number | null,
  requestedMode: FermiReferenceMode = "bands",
  zeroEnergyReferenceMode: ZeroEnergyReferenceMode = "calculated-fermi",
  selectedValenceBandMaximum?: number | null,
): [number, number] {
  const data = normalizeBandPlotData(inputData);
  const { zeroEnergy } = resolveBandPlotEnergyReference(
    data,
    scfFermiEnergy,
    requestedMode,
    zeroEnergyReferenceMode,
    selectedValenceBandMaximum,
  );
  let eMin = data.energy_range[0] - zeroEnergy;
  let eMax = data.energy_range[1] - zeroEnergy;
  const span = Math.max(eMax - eMin, 1e-6);
  const padding = span * 0.1;
  eMin -= padding;
  eMax += padding;

  const maxRange = 20;
  if (eMax - eMin > maxRange * 2) {
    eMin = -maxRange;
    eMax = maxRange;
  }

  return [eMin, eMax];
}

function extractSelectedValenceBandIndex(
  parameters?: Record<string, unknown> | null,
  bandCount?: number,
): number | null {
  const rawMetadata = parameters?.band_viewer_metadata;
  if (!rawMetadata || typeof rawMetadata !== "object") {
    return null;
  }

  const rawIndex = (rawMetadata as { selected_valence_band_index?: unknown })
    .selected_valence_band_index;
  const index = Number(rawIndex);
  if (!Number.isInteger(index) || index < 0) {
    return null;
  }
  if (bandCount != null && index >= bandCount) {
    return null;
  }
  return index;
}

function findSelectedValenceBandMaximum(
  energies: number[][],
  kPoints: number[],
  selectedBandIndex: number | null,
): { bandIndex: number; vbmEnergy: number; vbmK: number } | null {
  if (selectedBandIndex == null || selectedBandIndex < 0 || selectedBandIndex >= energies.length) {
    return null;
  }

  const band = energies[selectedBandIndex];
  const pointCount = Math.min(band.length, kPoints.length);
  if (pointCount === 0) {
    return null;
  }

  let vbmEnergy = Number.NEGATIVE_INFINITY;
  let vbmK = 0;
  for (let index = 0; index < pointCount; index += 1) {
    const energy = band[index];
    if (!Number.isFinite(energy)) {
      continue;
    }
    if (energy > vbmEnergy) {
      vbmEnergy = energy;
      vbmK = kPoints[index];
    }
  }

  if (!Number.isFinite(vbmEnergy)) {
    return null;
  }

  return {
    bandIndex: selectedBandIndex,
    vbmEnergy,
    vbmK,
  };
}

export function calculateDisplayedBandGapFromSelectedValenceBand(
  energies: number[][],
  kPoints: number[],
  selectedBandIndex: number | null,
  zeroEnergy = 0,
): BandGap | null {
  const selectedValenceBand = findSelectedValenceBandMaximum(energies, kPoints, selectedBandIndex);
  if (!selectedValenceBand) {
    return null;
  }

  let cbmEnergy = Number.POSITIVE_INFINITY;
  let cbmK = 0;

  for (const band of energies) {
    const pointCount = Math.min(band.length, kPoints.length);
    for (let index = 0; index < pointCount; index += 1) {
      const energy = band[index];
      if (!Number.isFinite(energy) || energy <= selectedValenceBand.vbmEnergy + BAND_GAP_TOLERANCE_EV) {
        continue;
      }
      if (energy < cbmEnergy) {
        cbmEnergy = energy;
        cbmK = kPoints[index];
      }
    }
  }

  if (!Number.isFinite(cbmEnergy)) {
    return null;
  }

  const gapValue = cbmEnergy - selectedValenceBand.vbmEnergy;
  if (!(gapValue > BAND_GAP_TOLERANCE_EV)) {
    return null;
  }

  return {
    value: gapValue,
    is_direct: Math.abs(selectedValenceBand.vbmK - cbmK) < DIRECT_GAP_K_TOLERANCE,
    vbm_k: selectedValenceBand.vbmK,
    cbm_k: cbmK,
    vbm_energy: selectedValenceBand.vbmEnergy - zeroEnergy,
    cbm_energy: cbmEnergy - zeroEnergy,
  };
}

interface BandEdgeSummary {
  vbmEnergy: number;
  vbmK: number;
  cbmEnergy: number;
  cbmK: number;
  hasCrossingBand: boolean;
}

function findBandEdges(
  energies: number[][],
  kPoints: number[],
  referenceEnergy: number,
): BandEdgeSummary | null {
  if (energies.length === 0 || kPoints.length === 0) {
    return null;
  }

  let vbmEnergy = Number.NEGATIVE_INFINITY;
  let vbmK = 0;
  let cbmEnergy = Number.POSITIVE_INFINITY;
  let cbmK = 0;
  let hasCrossingBand = false;

  for (const band of energies) {
    const pointCount = Math.min(band.length, kPoints.length);
    if (pointCount === 0) continue;

    let bandMin = Number.POSITIVE_INFINITY;
    let bandMax = Number.NEGATIVE_INFINITY;

    for (let index = 0; index < pointCount; index += 1) {
      const energy = band[index];
      if (!Number.isFinite(energy)) continue;

      if (energy < bandMin) bandMin = energy;
      if (energy > bandMax) bandMax = energy;

      if (energy <= referenceEnergy + BAND_GAP_TOLERANCE_EV && energy > vbmEnergy) {
        vbmEnergy = energy;
        vbmK = kPoints[index];
      }

      if (energy >= referenceEnergy - BAND_GAP_TOLERANCE_EV && energy < cbmEnergy) {
        cbmEnergy = energy;
        cbmK = kPoints[index];
      }
    }

    if (
      Number.isFinite(bandMin) &&
      Number.isFinite(bandMax) &&
      bandMin < referenceEnergy - BAND_GAP_TOLERANCE_EV &&
      bandMax > referenceEnergy + BAND_GAP_TOLERANCE_EV
    ) {
      hasCrossingBand = true;
    }
  }

  if (!Number.isFinite(vbmEnergy) || !Number.isFinite(cbmEnergy)) {
    return null;
  }

  return {
    vbmEnergy,
    vbmK,
    cbmEnergy,
    cbmK,
    hasCrossingBand,
  };
}

export function resolveBandPlotEnergyReference(
  inputData: BandPlotData,
  scfFermiEnergy?: number | null,
  requestedFermiMode: FermiReferenceMode = "bands",
  requestedZeroMode: ZeroEnergyReferenceMode = "calculated-fermi",
  selectedValenceBandMaximum?: number | null,
): {
  fermiEnergy: number;
  fermiMode: FermiReferenceMode;
  hasScfFermi: boolean;
  hasBandsFermi: boolean;
  zeroEnergy: number;
  zeroEnergyReferenceMode: ZeroEnergyReferenceMode;
  valenceBandMaximum: number | null;
} {
  const data = normalizeBandPlotData(inputData);
  const {
    fermiEnergy,
    mode: fermiMode,
    hasScfFermi,
    hasBandsFermi,
  } = resolveBandPlotFermiContext(data, scfFermiEnergy, requestedFermiMode);
  const explicitValenceBandMaximum = Number.isFinite(selectedValenceBandMaximum)
    ? (selectedValenceBandMaximum as number)
    : null;
  const resolvedValenceBandMaximum = explicitValenceBandMaximum;
  const zeroEnergyReferenceMode = requestedZeroMode === "vbm" && resolvedValenceBandMaximum != null
    ? "vbm"
    : "calculated-fermi";
  const zeroEnergy = zeroEnergyReferenceMode === "vbm"
    ? resolvedValenceBandMaximum!
    : fermiEnergy;

  return {
    fermiEnergy,
    fermiMode,
    hasScfFermi,
    hasBandsFermi,
    zeroEnergy,
    zeroEnergyReferenceMode,
    valenceBandMaximum: resolvedValenceBandMaximum,
  };
}

export function calculateDisplayedBandGap(
  energies: number[][],
  kPoints: number[],
  referenceEnergy: number,
  zeroEnergy: number = referenceEnergy,
): BandGap | null {
  const bandEdges = findBandEdges(energies, kPoints, referenceEnergy);
  if (!bandEdges || bandEdges.hasCrossingBand) {
    return null;
  }

  const gapValue = bandEdges.cbmEnergy - bandEdges.vbmEnergy;
  if (!(gapValue > BAND_GAP_TOLERANCE_EV)) {
    return null;
  }

  return {
    value: gapValue,
    is_direct: Math.abs(bandEdges.vbmK - bandEdges.cbmK) < DIRECT_GAP_K_TOLERANCE,
    vbm_k: bandEdges.vbmK,
    cbm_k: bandEdges.cbmK,
    vbm_energy: bandEdges.vbmEnergy - zeroEnergy,
    cbm_energy: bandEdges.cbmEnergy - zeroEnergy,
  };
}

function createZeroWeightGrid(nBands: number, nKpoints: number): number[][] {
  return Array.from({ length: nBands }, () => Array(nKpoints).fill(0));
}

function addWeightsInPlace(target: number[][], source: number[][]): void {
  for (let bandIndex = 0; bandIndex < target.length; bandIndex += 1) {
    const sourceBand = source[bandIndex];
    if (!sourceBand) continue;
    for (let kIndex = 0; kIndex < target[bandIndex].length; kIndex += 1) {
      const value = sourceBand[kIndex];
      if (Number.isFinite(value)) {
        target[bandIndex][kIndex] += value;
      }
    }
  }
}

function parseElementIdentityFromGroup(group: BandProjectionGroup): {
  key: string;
  display: string;
} {
  const label = group.label.trim();
  const symbolFromParentheses = label.match(/\(([A-Za-z][A-Za-z]?)\)/)?.[1];
  const symbolFromElementId = group.id.startsWith("element-")
    ? group.id.slice("element-".length)
    : "";
  const fallbackToken = label.split(/\s+/)[0] ?? group.id;
  const rawSymbol =
    symbolFromParentheses ||
    symbolFromElementId ||
    fallbackToken.match(/^([A-Za-z][A-Za-z]?)/)?.[1] ||
    "X";

  const normalized =
    rawSymbol.length <= 1
      ? rawSymbol.toUpperCase()
      : rawSymbol.charAt(0).toUpperCase() + rawSymbol.slice(1).toLowerCase();
  return { key: normalized.toLowerCase(), display: normalized };
}

function aggregateElementProjectionGroups(
  groups: BandProjectionGroup[],
  nBands: number,
  nKpoints: number,
): BandProjectionGroup[] {
  const grouped = new Map<string, BandProjectionGroup>();
  for (const group of groups) {
    const element = parseElementIdentityFromGroup(group);
    let aggregate = grouped.get(element.key);
    if (!aggregate) {
      aggregate = {
        id: `element-${element.key}`,
        label: `${element.display} total`,
        kind: "atom",
        weights: createZeroWeightGrid(nBands, nKpoints),
      };
      grouped.set(element.key, aggregate);
    }
    addWeightsInPlace(aggregate.weights, group.weights);
  }

  return Array.from(grouped.values()).sort((a, b) => a.label.localeCompare(b.label));
}

function buildTotalProjectionGroup(
  groups: BandProjectionGroup[],
  id: string,
  label: string,
  nBands: number,
  nKpoints: number,
  kind: BandProjectionGroup["kind"] = "orbital",
): BandProjectionGroup | null {
  if (groups.length === 0) return null;
  const totalGroup: BandProjectionGroup = {
    id,
    label,
    kind,
    weights: createZeroWeightGrid(nBands, nKpoints),
  };
  for (const group of groups) {
    addWeightsInPlace(totalGroup.weights, group.weights);
  }
  return totalGroup;
}

function parseElementOrbitalIdentity(group: BandProjectionGroup): {
  elementKey: string;
  elementLabel: string;
  orbitalKey: string;
} {
  const idMatch = group.id.match(/^element-orbital-([a-z0-9_+]+)-([a-z0-9_+\-]+)$/i);
  if (idMatch) {
    const elementKey = idMatch[1].toLowerCase();
    const elementLabel = elementKey.length <= 1
      ? elementKey.toUpperCase()
      : elementKey.charAt(0).toUpperCase() + elementKey.slice(1).toLowerCase();
    return {
      elementKey,
      elementLabel,
      orbitalKey: idMatch[2].toLowerCase(),
    };
  }

  const labelMatch = group.label.match(/^([A-Za-z][A-Za-z]?)\s+([A-Za-z0-9]+)\s+orbitals$/i);
  if (labelMatch) {
    return {
      elementKey: labelMatch[1].toLowerCase(),
      elementLabel:
        labelMatch[1].length <= 1
          ? labelMatch[1].toUpperCase()
          : labelMatch[1].charAt(0).toUpperCase() + labelMatch[1].slice(1).toLowerCase(),
      orbitalKey: labelMatch[2].toLowerCase(),
    };
  }

  return {
    elementKey: "all",
    elementLabel: "All elements",
    orbitalKey: "other",
  };
}

function orbitalSortRank(group: BandProjectionGroup): number {
  const orbitalKey = parseElementOrbitalIdentity(group).orbitalKey;
  const orbitalOrder: Record<string, number> = {
    s: 0,
    p: 1,
    d: 2,
    f: 3,
    g: 4,
    other: 99,
  };
  if (orbitalKey in orbitalOrder) {
    return orbitalOrder[orbitalKey];
  }
  const lMatch = orbitalKey.match(/^l(\d+)$/i);
  if (lMatch) {
    return 10 + Number.parseInt(lMatch[1], 10);
  }
  return 98;
}

function formatProjectionOrbitalLabel(orbitalKey: string): string {
  return orbitalKey.trim().toLowerCase() || "other";
}

function formatGlobalOrbitalLabel(group: BandProjectionGroup): string {
  const idMatch = group.id.match(/(?:^|-)orbital-([a-z0-9_+\-]+)$/i);
  if (idMatch) {
    return formatProjectionOrbitalLabel(idMatch[1]);
  }

  const labelMatch = group.label.trim().match(/^([A-Za-z0-9_+\-]+)\s+orbitals$/i);
  if (labelMatch) {
    return formatProjectionOrbitalLabel(labelMatch[1]);
  }

  return group.label.trim() || group.id;
}

function buildProjectionSelectionEntries(data: BandData): ProjectionSelectionEntry[] {
  const hasElementResolvedOrbitals = (data.projections?.element_orbital_groups?.length ?? 0) > 0;
  const elementSelections = new Map<string, {
    elementKey: string;
    elementLabel: string;
    fullSelectionValue: string | null;
    fullSelectionMode: ProjectionMode;
    orbitalGroups: BandProjectionGroup[];
  }>();

  for (const group of data.projections?.atom_groups ?? []) {
    const element = parseElementIdentityFromGroup(group);
    const existing = elementSelections.get(element.key);
    if (existing) {
      if (!existing.fullSelectionValue) {
        existing.fullSelectionValue = `element-${element.key}`;
        existing.fullSelectionMode = "atom";
      }
      continue;
    }

    elementSelections.set(element.key, {
      elementKey: element.key,
      elementLabel: element.display,
      fullSelectionValue: `element-${element.key}`,
      fullSelectionMode: "atom",
      orbitalGroups: [],
    });
  }

  for (const group of data.projections?.element_orbital_groups ?? []) {
    const parsed = parseElementOrbitalIdentity(group);
    const existing = elementSelections.get(parsed.elementKey);
    if (existing) {
      existing.orbitalGroups.push(group);
      if (!existing.fullSelectionValue) {
        existing.fullSelectionValue = `element-orbital-total-${parsed.elementKey}`;
        existing.fullSelectionMode = "orbital";
      }
      continue;
    }

    elementSelections.set(parsed.elementKey, {
      elementKey: parsed.elementKey,
      elementLabel: parsed.elementLabel,
      fullSelectionValue: `element-orbital-total-${parsed.elementKey}`,
      fullSelectionMode: "orbital",
      orbitalGroups: [group],
    });
  }

  const entries: ProjectionSelectionEntry[] = [];
  const sortedElements = Array.from(elementSelections.values()).sort((a, b) =>
    a.elementLabel.localeCompare(b.elementLabel),
  );
  for (const element of sortedElements) {
    if (element.fullSelectionValue) {
      entries.push({
        value: element.fullSelectionValue,
        label: element.elementLabel,
        mode: element.fullSelectionMode,
        elementKey: element.fullSelectionMode === "orbital" ? element.elementKey : undefined,
      });
    }

    const sortedOrbitalGroups = [...element.orbitalGroups].sort((a, b) => {
      const rankDiff = orbitalSortRank(a) - orbitalSortRank(b);
      if (rankDiff !== 0) return rankDiff;
      return a.label.localeCompare(b.label);
    });
    for (const group of sortedOrbitalGroups) {
      const orbital = parseElementOrbitalIdentity(group);
      entries.push({
        value: group.id,
        label: `${element.elementLabel}-${formatProjectionOrbitalLabel(orbital.orbitalKey)}`,
        mode: "orbital",
        elementKey: orbital.elementKey,
      });
    }
  }

  if (entries.length > 0 && hasElementResolvedOrbitals) {
    return entries;
  }

  const globalOrbitalGroups = [...(data.projections?.orbital_groups ?? [])].sort((a, b) => {
    const rankDiff = orbitalSortRank(a) - orbitalSortRank(b);
    if (rankDiff !== 0) return rankDiff;
    return a.label.localeCompare(b.label);
  });
  if (globalOrbitalGroups.length === 0) {
    return [];
  }

  return [
    ...entries,
    {
      value: "orbital-total",
      label: "all orbitals",
      mode: "orbital",
      elementKey: "",
    },
    ...globalOrbitalGroups.map((group) => ({
      value: group.id,
      label: formatGlobalOrbitalLabel(group),
      mode: "orbital" as const,
      elementKey: "",
    })),
  ];
}

export function buildBandPlotProjectionOptions(inputData: BandPlotData): BandPlotProjectionOption[] {
  const data = normalizeBandPlotData(inputData);
  return [
    { value: "none", label: "none" },
    ...buildProjectionSelectionEntries(data).map(({ value, label }) => ({ value, label })),
  ];
}

export function BandPlot({
  data: inputData,
  width = 700,
  height = 500,
  energyRange,
  showFermiLevel = true,
  scfFermiEnergy,
  yAxisLabel = "E − E_F (eV)",
  pointLabel = "Band index",
  valueLabel = "E − E_F",
  valueUnit = "eV",
  valueDecimals = 3,
  primaryCountLabel = "bands",
  secondaryCountLabel = "k-points",
  scrollHint = "Scroll: zoom Y | Shift+Scroll: pan energy",
  yClampRange = [-25, 25],
  viewerType = "electronic",
  sharedSettings = null,
  showSidebar = true,
  projectionSelection,
  enableWheelRangeControl = true,
  enableHoverScrollLock = true,
  comparisonOptions,
  comparisonTitle = "Band Comparison",
  comparisonNoneLabel = "None",
  windowOverlays,
  onWindowOverlayChange,
  windowOverlayHint = "Drag window edges or side sliders to tune selected ranges.",
  showBandGapOverlayOverride,
  calculationParameters = null,
  onPersistSelectedValenceBandIndex,
}: BandPlotProps) {
  const { isDark } = useTheme();
  const data = useMemo(() => normalizeBandPlotData(inputData), [inputData]);
  const colors = useMemo(() => isDark
    ? { bg: "#1e1e2e", axis: "#718096", grid: "#4a5568", text: "#e2e8f0", tooltip: "#2d3748", tooltipBorder: "#4a5568", tooltipText: "#e2e8f0" }
    : { bg: "#ffffff", axis: "#333", grid: "#999", text: "#000", tooltip: "#fff", tooltipBorder: "#ccc", tooltipText: "#333" },
  [isDark]);

  const svgRef = useRef<SVGSVGElement>(null);
  const plotCanvasRef = useRef<HTMLDivElement>(null);
  const clipPathId = useId();
  const projectionGradientId = useId();
  const [hoveredPoint, setHoveredPoint] = useState<HoveredPoint | null>(null);
  const [isHoveringPlot, setIsHoveringPlot] = useState(false);

  // Self-sizing via ResizeObserver
  const [containerSize, setContainerSize] = useState<{ width: number; height: number } | null>(null);

  useLayoutEffect(() => {
    const el = plotCanvasRef.current;
    if (!el) return;

    const observer = new ResizeObserver((entries) => {
      const entry = entries[0];
      if (entry) {
        const nextWidth = Math.floor(entry.contentRect.width);
        const nextHeight = Math.floor(entry.contentRect.height);
        setContainerSize((prev) => {
          if (prev && prev.width === nextWidth && prev.height === nextHeight) {
            return prev;
          }
          return {
            width: nextWidth,
            height: nextHeight,
          };
        });
      }
    });
    observer.observe(el);
    return () => observer.disconnect();
  }, []);

  // Y-axis energy window (adjustable via scroll)
  const [yMin, setYMin] = useState<number | null>(null);
  const [yMax, setYMax] = useState<number | null>(null);
  const [manualYMinInput, setManualYMinInput] = useState("");
  const [manualYMaxInput, setManualYMaxInput] = useState("");
  const [manualYTickStepInput, setManualYTickStepInput] = useState("");
  const [manualRangeError, setManualRangeError] = useState<string | null>(null);

  // Appearance controls
  const [lineWidth, setLineWidth] = useState(1.5);
  const [lineOpacity, setLineOpacity] = useState(0.85);
  const [plotTextScale, setPlotTextScale] = useState(1);
  const [colorMode, setColorMode] = useState<ColorMode>("rainbow");
  const [singleBandColor, setSingleBandColor] = useState("#1565c0");
  const [rainbowPalette, setRainbowPalette] = useState<RainbowPalette>("jet");

  // Fat-band controls
  const [fatBandsEnabled, setFatBandsEnabled] = useState(false);
  const [projectionRenderMode, setProjectionRenderMode] = useState<ProjectionRenderMode>("fat");
  const [projectionMode, setProjectionMode] = useState<ProjectionMode>("atom");
  const [selectedOrbitalElementKey, setSelectedOrbitalElementKey] = useState("");
  const [selectedProjectionId, setSelectedProjectionId] = useState("");
  const [projectionNormalizeMode, setProjectionNormalizeMode] =
    useState<ProjectionNormalizeMode>("global");
  const [fatScale, setFatScale] = useState(8);
  const [fatOpacity, setFatOpacity] = useState(0.45);
  const [fatColorMode, setFatColorMode] = useState<FatColorMode>("band");
  const [fatAccentColor, setFatAccentColor] = useState("#f57c00");
  const [showLinesWithFat, setShowLinesWithFat] = useState(true);
  const [projectionColorMin, setProjectionColorMin] = useState(0.2);
  const [projectionColorMax, setProjectionColorMax] = useState(0.4);
  const [projectionLegendTickCount, setProjectionLegendTickCount] = useState(9);

  // Plot background toggle
  const [plotBgWhite, setPlotBgWhite] = useState(true);
  const [showBandGapOverlay, setShowBandGapOverlay] = useState(true);

  // UI section toggles
  const [appearanceExpanded, setAppearanceExpanded] = useState(true);
  const [comparisonExpanded, setComparisonExpanded] = useState(true);
  const [bandGapExpanded, setBandGapExpanded] = useState(true);
  const [projectionExpanded, setProjectionExpanded] = useState(false);
  const [exportNote, setExportNote] = useState("");
  const [selectedComparisonId, setSelectedComparisonId] = useState("");
  const [overlayDragState, setOverlayDragState] = useState<OverlayDragState | null>(null);

  const requestedFermiReferenceMode = sharedSettings?.fermiReferenceMode ?? null;
  const requestedZeroEnergyReferenceMode = sharedSettings?.zeroEnergyReferenceMode ?? null;
  const fallbackFermiReferenceMode = (
    scfFermiEnergy != null && Number.isFinite(scfFermiEnergy)
  ) ? "scf" : "bands";
  const [fermiReferenceMode, setFermiReferenceMode] = useState<FermiReferenceMode>(
    requestedFermiReferenceMode ?? fallbackFermiReferenceMode,
  );
  const [zeroEnergyReferenceMode, setZeroEnergyReferenceMode] = useState<ZeroEnergyReferenceMode>(
    requestedZeroEnergyReferenceMode ?? "calculated-fermi",
  );
  const persistedSelectedValenceBandIndex = useMemo(
    () => extractSelectedValenceBandIndex(calculationParameters, data.energies.length),
    [calculationParameters, data.energies.length],
  );
  const [selectedValenceBandIndex, setSelectedValenceBandIndex] = useState<number | null>(
    persistedSelectedValenceBandIndex,
  );
  const [isSelectingValenceBand, setIsSelectingValenceBand] = useState(false);
  const [valenceBandSelectionStatus, setValenceBandSelectionStatus] = useState<string | null>(null);

  useEffect(() => {
    setSelectedValenceBandIndex(persistedSelectedValenceBandIndex);
  }, [persistedSelectedValenceBandIndex]);

  useEffect(() => {
    if (sharedSettings) {
      return;
    }
    const hasScfFermi = scfFermiEnergy != null && Number.isFinite(scfFermiEnergy);
    const hasBandsFermi = Number.isFinite(data.fermi_energy);
    if (fermiReferenceMode === "scf" && !hasScfFermi) {
      setFermiReferenceMode("bands");
      return;
    }
    if (fermiReferenceMode === "bands" && !hasBandsFermi && hasScfFermi) {
      setFermiReferenceMode("scf");
    }
  }, [data.fermi_energy, fermiReferenceMode, scfFermiEnergy, sharedSettings]);

  useEffect(() => {
    const requestedZeroMode = requestedZeroEnergyReferenceMode ?? zeroEnergyReferenceMode;
    if (viewerType !== "electronic" || requestedZeroMode !== "vbm") {
      setIsSelectingValenceBand(false);
      setValenceBandSelectionStatus(null);
      return;
    }
    if (selectedValenceBandIndex == null) {
      setIsSelectingValenceBand(true);
      if (!valenceBandSelectionStatus) {
        setValenceBandSelectionStatus("Click the valence band on the plot.");
      }
    }
  }, [
    requestedZeroEnergyReferenceMode,
    selectedValenceBandIndex,
    valenceBandSelectionStatus,
    viewerType,
    zeroEnergyReferenceMode,
  ]);

  const {
    fermiEnergy,
    fermiMode: resolvedFermiReferenceMode,
    hasScfFermi,
    hasBandsFermi,
    zeroEnergy,
    zeroEnergyReferenceMode: resolvedZeroEnergyReferenceMode,
    valenceBandMaximum,
  } = useMemo(() => {
    const selectedValenceBand = findSelectedValenceBandMaximum(
      data.energies,
      data.k_points,
      selectedValenceBandIndex,
    );
    return resolveBandPlotEnergyReference(
      data,
      scfFermiEnergy,
      requestedFermiReferenceMode ?? fermiReferenceMode,
      requestedZeroEnergyReferenceMode ?? zeroEnergyReferenceMode,
      selectedValenceBand?.vbmEnergy ?? null,
    );
  },
    [
      data,
      fermiReferenceMode,
      requestedFermiReferenceMode,
      requestedZeroEnergyReferenceMode,
      scfFermiEnergy,
      selectedValenceBandIndex,
      zeroEnergyReferenceMode,
    ],
  );
  const selectedValenceBand = useMemo(
    () => findSelectedValenceBandMaximum(data.energies, data.k_points, selectedValenceBandIndex),
    [data.energies, data.k_points, selectedValenceBandIndex],
  );

  const effectiveLineWidth = sharedSettings?.lineWidth ?? lineWidth;
  const effectiveLineOpacity = sharedSettings?.lineOpacity ?? lineOpacity;
  const effectivePlotTextScale = sharedSettings?.plotTextScale ?? plotTextScale;
  const effectiveColorMode = sharedSettings?.colorMode ?? colorMode;
  const effectiveSingleBandColor = sharedSettings?.singleBandColor ?? singleBandColor;
  const effectiveRainbowPalette = sharedSettings?.rainbowPalette ?? rainbowPalette;
  const effectivePlotBgWhite = sharedSettings?.plotBgWhite ?? plotBgWhite;
  const effectiveShowBandGapOverlay = showBandGapOverlayOverride
    ?? sharedSettings?.showBandGapOverlay
    ?? showBandGapOverlay;

  const activeFermiSourceLabel = resolvedFermiReferenceMode === "scf" ? "SCF" : "Bands run";
  const activeFermiDisplay = Number.isFinite(fermiEnergy) ? `${fermiEnergy.toFixed(3)} eV` : "N/A";
  const selectedValenceBandLabel = selectedValenceBandIndex != null
    ? `Band ${selectedValenceBandIndex + 1}`
    : "Not selected";
  const activeZeroReferenceLabel = resolvedZeroEnergyReferenceMode === "vbm" ? "VBM" : "Calculated E_F";
  const activeZeroReferenceDisplay = resolvedZeroEnergyReferenceMode === "vbm" && valenceBandMaximum != null
    ? `${valenceBandMaximum.toFixed(3)} eV`
    : activeFermiDisplay;
  const zeroLineLabel = resolvedZeroEnergyReferenceMode === "vbm" ? "VBM" : "E_F";
  const resolvedYAxisLabel = viewerType === "electronic" && yAxisLabel === "E − E_F (eV)"
    ? (resolvedZeroEnergyReferenceMode === "vbm" ? "E − E_VBM (eV)" : yAxisLabel)
    : yAxisLabel;
  const resolvedValueLabel = viewerType === "electronic" && valueLabel === "E − E_F"
    ? (resolvedZeroEnergyReferenceMode === "vbm" ? "E − E_VBM" : valueLabel)
    : valueLabel;
  const availableComparisonOptions = useMemo<NormalizedBandPlotComparisonOption[]>(
    () => (comparisonOptions ?? []).map((option) => ({
      ...option,
      data: normalizeBandPlotData(option.data),
    })),
    [comparisonOptions],
  );
  const hasComparisonControls = comparisonOptions !== undefined;
  const editableWindowOverlays = useMemo(
    () =>
      (windowOverlays ?? [])
        .filter((overlay) =>
          Number.isFinite(overlay.min)
          && Number.isFinite(overlay.max)
          && Number.isFinite(overlay.max - overlay.min),
        )
        .map((overlay) => ({
          ...overlay,
          min: Math.min(overlay.min, overlay.max),
          max: Math.max(overlay.min, overlay.max),
        })),
    [windowOverlays],
  );
  const windowOverlayById = useMemo(
    () => new Map(editableWindowOverlays.map((overlay) => [overlay.id, overlay])),
    [editableWindowOverlays],
  );
  const overlayEditorEnabled = editableWindowOverlays.length > 0 && typeof onWindowOverlayChange === "function";

  useEffect(() => {
    if (availableComparisonOptions.length === 0) {
      setSelectedComparisonId("");
      return;
    }
    if (!selectedComparisonId) {
      return;
    }
    const stillValid = availableComparisonOptions.some((option) => option.id === selectedComparisonId);
    if (!stillValid) {
      setSelectedComparisonId("");
    }
  }, [availableComparisonOptions, selectedComparisonId]);

  useEffect(() => {
    if (!overlayDragState) {
      return;
    }
    if (!windowOverlayById.has(overlayDragState.overlayId)) {
      setOverlayDragState(null);
    }
  }, [overlayDragState, windowOverlayById]);

  // Shift all energies relative to the selected zero-energy reference.
  const shiftedEnergies = useMemo(() => {
    return data.energies.map((band) => band.map((e) => e - zeroEnergy));
  }, [data.energies, zeroEnergy]);

  const selectedComparison = useMemo(
    () => availableComparisonOptions.find((option) => option.id === selectedComparisonId) ?? null,
    [availableComparisonOptions, selectedComparisonId],
  );

  const comparisonShiftedEnergies = useMemo(() => {
    if (!selectedComparison) {
      return [];
    }
    return selectedComparison.data.energies.map((band) => band.map((e) => e - zeroEnergy));
  }, [selectedComparison, zeroEnergy]);

  const comparisonKPointsForPlot = useMemo(() => {
    if (!selectedComparison) {
      return [];
    }
    return remapComparisonKPoints(data, selectedComparison.data) ?? selectedComparison.data.k_points;
  }, [data, selectedComparison]);

  const displayedBandGap = useMemo(() => {
    if (viewerType !== "electronic") return null;
    if (resolvedZeroEnergyReferenceMode === "vbm") {
      return calculateDisplayedBandGapFromSelectedValenceBand(
        data.energies,
        data.k_points,
        selectedValenceBandIndex,
        zeroEnergy,
      );
    }
    return calculateDisplayedBandGap(data.energies, data.k_points, fermiEnergy, zeroEnergy);
  }, [
    data.energies,
    data.k_points,
    fermiEnergy,
    resolvedZeroEnergyReferenceMode,
    selectedValenceBandIndex,
    viewerType,
    zeroEnergy,
  ]);

  const yDomain = useMemo((): [number, number] => {
    if (yMin !== null && yMax !== null) {
      return [yMin, yMax];
    }
    if (energyRange) {
      return energyRange;
    }
    return getDefaultBandPlotEnergyRange(
      data,
      scfFermiEnergy,
      requestedFermiReferenceMode ?? fermiReferenceMode,
      requestedZeroEnergyReferenceMode ?? zeroEnergyReferenceMode,
      selectedValenceBand?.vbmEnergy ?? null,
    );
  }, [
    data,
    energyRange,
    fermiReferenceMode,
    requestedFermiReferenceMode,
    requestedZeroEnergyReferenceMode,
    scfFermiEnergy,
    selectedValenceBand,
    zeroEnergyReferenceMode,
    yMax,
    yMin,
  ]);

  useEffect(() => {
    if (!showSidebar) {
      return;
    }
    setManualYMinInput(formatAxisInputValue(yDomain[0]));
    setManualYMaxInput(formatAxisInputValue(yDomain[1]));
  }, [showSidebar, yDomain[0], yDomain[1]]);

  const yTickStep = useMemo(() => {
    const trimmed = manualYTickStepInput.trim();
    if (!trimmed) {
      return null;
    }
    const parsed = Number.parseFloat(trimmed);
    return Number.isFinite(parsed) && parsed > 0 ? parsed : null;
  }, [manualYTickStepInput]);

  const axisTickFontSize = Math.max(8, 11 * effectivePlotTextScale);
  const axisLabelFontSize = Math.max(10, 14 * effectivePlotTextScale);
  const symmetryLabelFontSize = axisLabelFontSize;
  const symmetryLabelYOffset = Math.max(20, symmetryLabelFontSize * 1.35);
  const yTickLabelYOffset = axisTickFontSize * 0.35;
  const yTicks = useMemo(
    () => getYAxisTicks(yDomain[0], yDomain[1], yTickStep),
    [yDomain, yTickStep],
  );
  const yTickDecimals = useMemo(() => getTickDecimals(yTickStep), [yTickStep]);
  const maxYTickLabelChars = useMemo(() => {
    if (yTicks.length === 0) return 3;
    let maxChars = 0;
    for (const tick of yTicks) {
      maxChars = Math.max(maxChars, tick.toFixed(yTickDecimals).length);
    }
    return maxChars;
  }, [yTickDecimals, yTicks]);
  const estimatedYTickLabelWidth = Math.max(
    axisTickFontSize * 3,
    maxYTickLabelChars * axisTickFontSize * 0.62,
  );
  const fermiLabelFontSize = axisTickFontSize;
  const fermiLabelYOffset = fermiLabelFontSize * 0.35;
  const projectionLegendReservedWidth =
    viewerType === "electronic" && fatBandsEnabled && projectionRenderMode === "color"
      ? Math.round(Math.max(84, 92 * effectivePlotTextScale))
      : 0;

  // Margins & dimensions — scale from actual axis text footprint.
  const margin = useMemo(() => ({
    top: 30,
    right: Math.round(Math.max(
      30,
      16 + fermiLabelFontSize * 2.2,
      projectionLegendReservedWidth,
    )),
    bottom: Math.round(Math.max(50, 18 + symmetryLabelFontSize * 2.2)),
    left: Math.round(
      Math.max(70, 20 + estimatedYTickLabelWidth + axisLabelFontSize * 1.6),
    ),
  }), [
    axisLabelFontSize,
    estimatedYTickLabelWidth,
    fermiLabelFontSize,
    projectionLegendReservedWidth,
    symmetryLabelFontSize,
  ]);

  const fallbackWidth = Math.max(240, width);
  const fallbackHeight = Math.max(220, height);
  const resolvedWidth = containerSize ? Math.max(1, containerSize.width) : fallbackWidth;
  const resolvedHeight = containerSize ? Math.max(1, containerSize.height) : fallbackHeight;
  const plotWidth = Math.max(1, resolvedWidth - margin.left - margin.right);
  const plotHeight = Math.max(1, resolvedHeight - margin.top - margin.bottom);
  const yAxisLabelXOffset = Math.max(
    axisLabelFontSize + 14,
    Math.min(
      margin.left - 8,
      estimatedYTickLabelWidth + axisLabelFontSize * 0.95 + 16,
    ),
  );
  const yTickLabelX = -(8 + axisTickFontSize * 0.18);
  const yTickMarkX = -(4 + axisTickFontSize * 0.1);
  const fermiLabelXOffset = Math.max(5, fermiLabelFontSize * 0.45);

  // Calculate scales - X is always fixed, Y is adjustable
  const scales = useMemo(() => {
    const kMin = data.k_points.length > 0 ? data.k_points[0] : 0;
    const kMax = data.k_points[data.k_points.length - 1] ?? 1;
    const kSpan = Math.abs(kMax - kMin) > 1e-12 ? kMax - kMin : 1;
    const [eMin, eMax] = yDomain;

    const eSpan = Math.abs(eMax - eMin) > 1e-12 ? eMax - eMin : 1;

    return {
      kMin,
      kMax,
      eMin,
      eMax,
      xScale: (k: number) => ((k - kMin) / kSpan) * plotWidth,
      yScale: (e: number) => plotHeight - ((e - eMin) / eSpan) * plotHeight,
    };
  }, [data.k_points, plotWidth, plotHeight, yDomain]);

  const beginOverlayDrag = useCallback(
    (
      event: React.PointerEvent<SVGElement>,
      overlay: BandPlotWindowOverlay,
      mode: BandPlotWindowOverlayInteraction,
    ) => {
      if (!overlayEditorEnabled) return;
      event.preventDefault();
      event.stopPropagation();
      setOverlayDragState({
        overlayId: overlay.id,
        pointerId: event.pointerId,
        mode,
        startClientY: event.clientY,
        startMin: overlay.min,
        startMax: overlay.max,
      });
    },
    [overlayEditorEnabled],
  );

  useEffect(() => {
    if (!overlayDragState || !onWindowOverlayChange) {
      return;
    }

    const handlePointerMove = (event: PointerEvent) => {
      if (event.pointerId !== overlayDragState.pointerId) {
        return;
      }
      const activeOverlay = windowOverlayById.get(overlayDragState.overlayId);
      if (!activeOverlay) {
        return;
      }
      event.preventDefault();

      const span = Math.max(scales.eMax - scales.eMin, 1e-12);
      const deltaEnergy = -((event.clientY - overlayDragState.startClientY) / Math.max(plotHeight, 1)) * span;
      const minLimit = Math.min(scales.eMin, scales.eMax);
      const maxLimit = Math.max(scales.eMin, scales.eMax);
      const minSpan = Math.max(MIN_WINDOW_SPAN_EV, (maxLimit - minLimit) * 1e-6);

      let nextMin = overlayDragState.startMin;
      let nextMax = overlayDragState.startMax;

      if (overlayDragState.mode === "move") {
        const windowSpan = Math.max(overlayDragState.startMax - overlayDragState.startMin, minSpan);
        nextMin = overlayDragState.startMin + deltaEnergy;
        nextMax = overlayDragState.startMax + deltaEnergy;
        if (nextMin < minLimit) {
          nextMin = minLimit;
          nextMax = minLimit + windowSpan;
        }
        if (nextMax > maxLimit) {
          nextMax = maxLimit;
          nextMin = maxLimit - windowSpan;
        }
      } else if (overlayDragState.mode === "resize-min" || overlayDragState.mode === "slider-min") {
        nextMin = clampToRange(
          overlayDragState.startMin + deltaEnergy,
          minLimit,
          overlayDragState.startMax - minSpan,
        );
        nextMax = overlayDragState.startMax;
      } else {
        nextMax = clampToRange(
          overlayDragState.startMax + deltaEnergy,
          overlayDragState.startMin + minSpan,
          maxLimit,
        );
        nextMin = overlayDragState.startMin;
      }

      onWindowOverlayChange({
        id: activeOverlay.id,
        min: nextMin,
        max: nextMax,
        interaction: overlayDragState.mode,
      });
    };

    const handlePointerUp = (event: PointerEvent) => {
      if (event.pointerId !== overlayDragState.pointerId) {
        return;
      }
      setOverlayDragState(null);
    };

    window.addEventListener("pointermove", handlePointerMove);
    window.addEventListener("pointerup", handlePointerUp);
    window.addEventListener("pointercancel", handlePointerUp);

    return () => {
      window.removeEventListener("pointermove", handlePointerMove);
      window.removeEventListener("pointerup", handlePointerUp);
      window.removeEventListener("pointercancel", handlePointerUp);
    };
  }, [
    onWindowOverlayChange,
    overlayDragState,
    plotHeight,
    scales.eMax,
    scales.eMin,
    windowOverlayById,
  ]);

  const windowOverlayRenderData = useMemo(
    () =>
      editableWindowOverlays.map((overlay) => {
        const minEnergy = Math.min(overlay.min, overlay.max);
        const maxEnergy = Math.max(overlay.min, overlay.max);
        const topY = scales.yScale(maxEnergy);
        const bottomY = scales.yScale(minEnergy);
        const y = Math.min(topY, bottomY);
        const heightPx = Math.max(2, Math.abs(bottomY - topY));
        const sliderX = overlay.side === "left" ? -12 : plotWidth + 12;
        return {
          overlay,
          y,
          heightPx,
          topY,
          bottomY,
          sliderX,
        };
      }),
    [editableWindowOverlays, plotWidth, scales],
  );

  const bandColors = useMemo(
    () =>
      shiftedEnergies.map((_, idx) =>
        bandColorForIndex(
          idx,
          shiftedEnergies.length,
          effectiveColorMode,
          effectiveSingleBandColor,
          effectiveRainbowPalette,
        ),
      ),
    [
      shiftedEnergies,
      effectiveColorMode,
      effectiveSingleBandColor,
      effectiveRainbowPalette,
    ],
  );

  const comparisonModeActive = viewerType === "electronic" && selectedComparison !== null;
  const primaryBandStrokeColors = useMemo(
    () =>
      comparisonModeActive
        ? shiftedEnergies.map(() => "#1565c0")
        : bandColors,
    [bandColors, comparisonModeActive, shiftedEnergies],
  );
  const comparisonBandStrokeColor = "#d97706";

  const bandGapOverlay = useMemo(() => {
    if (viewerType !== "electronic" || !displayedBandGap) {
      return null;
    }

    const vbmX = scales.xScale(displayedBandGap.vbm_k);
    const cbmX = scales.xScale(displayedBandGap.cbm_k);
    const vbmY = scales.yScale(displayedBandGap.vbm_energy);
    const cbmY = scales.yScale(displayedBandGap.cbm_energy);
    const gapTopY = Math.min(vbmY, cbmY);
    const gapBottomY = Math.max(vbmY, cbmY);
    const gapMidY = (gapTopY + gapBottomY) / 2;
    const gapMidX = (vbmX + cbmX) / 2;
    const primaryLabel = `${displayedBandGap.value.toFixed(3)} eV`;
    const secondaryLabel = displayedBandGap.is_direct ? "Direct gap" : "Indirect gap";
    const gapColor = displayedBandGap.is_direct ? "#2e7d32" : "#00796b";
    const labelFontSize = Math.max(10, 11 * effectivePlotTextScale);
    const secondaryFontSize = Math.max(9, labelFontSize - 1);
    const edgeLabelFontSize = Math.max(9, 10 * effectivePlotTextScale);
    const approxCharWidth = labelFontSize * 0.58;
    const labelWidth = Math.max(
      118,
      Math.max(primaryLabel.length, secondaryLabel.length) * approxCharWidth + 18,
    );
    const labelHeight = labelFontSize + secondaryFontSize + 18;
    const labelX = clampToRange(gapMidX - labelWidth / 2, 8, plotWidth - labelWidth - 8);
    const labelY = clampToRange(gapMidY - labelHeight / 2, 8, plotHeight - labelHeight - 8);

    return {
      vbmX,
      cbmX,
      vbmY,
      cbmY,
      gapTopY,
      gapBottomY,
      gapMidX,
      gapMidY,
      labelX,
      labelY,
      labelWidth,
      labelHeight,
      labelFontSize,
      secondaryFontSize,
      edgeLabelFontSize,
      primaryLabel,
      secondaryLabel,
      gapColor,
    };
  }, [displayedBandGap, effectivePlotTextScale, plotHeight, plotWidth, scales, viewerType]);

  const bandGapOverlayVisible =
    viewerType === "electronic" && displayedBandGap !== null && effectiveShowBandGapOverlay;

  const hasProjectionData = useMemo(() => {
    const atomCount = data.projections?.atom_groups?.length ?? 0;
    const orbitalCount = data.projections?.orbital_groups?.length ?? 0;
    const elementOrbitalCount = data.projections?.element_orbital_groups?.length ?? 0;
    return atomCount > 0 || orbitalCount > 0 || elementOrbitalCount > 0;
  }, [data.projections]);

  const projectionSelectionEntries = useMemo(
    () => buildProjectionSelectionEntries(data),
    [data],
  );

  const orbitalElementOptions = useMemo((): OrbitalElementOption[] => {
    const groups = data.projections?.element_orbital_groups ?? [];
    if (groups.length === 0) return [];

    const grouped = new Map<string, OrbitalElementOption>();
    for (const group of groups) {
      const parsed = parseElementOrbitalIdentity(group);
      const existing = grouped.get(parsed.elementKey);
      if (existing) {
        existing.groups.push(group);
      } else {
        grouped.set(parsed.elementKey, {
          key: parsed.elementKey,
          label: parsed.elementLabel,
          groups: [group],
        });
      }
    }

    const options = Array.from(grouped.values());
    for (const option of options) {
      option.groups.sort((a, b) => {
        const rankDiff = orbitalSortRank(a) - orbitalSortRank(b);
        if (rankDiff !== 0) return rankDiff;
        return a.label.localeCompare(b.label);
      });
    }
    options.sort((a, b) => a.label.localeCompare(b.label));
    return options;
  }, [data.projections]);

  const hasElementResolvedOrbitals = orbitalElementOptions.length > 0;

  const allElementProjectionGroups = useMemo(() => {
    const atomGroups = data.projections?.atom_groups ?? [];
    return aggregateElementProjectionGroups(atomGroups, data.n_bands, data.n_kpoints);
  }, [data.n_bands, data.n_kpoints, data.projections]);

  const allElementOrbitalTotalGroups = useMemo(() => {
    if (!hasElementResolvedOrbitals) return [];

    return orbitalElementOptions
      .map((option) =>
        buildTotalProjectionGroup(
          option.groups,
          `element-orbital-total-${option.key}`,
          `${option.label} total (all orbitals)`,
          data.n_bands,
          data.n_kpoints,
          "element_orbital",
        ),
      )
      .filter((group): group is BandProjectionGroup => group !== null);
  }, [
    data.n_bands,
    data.n_kpoints,
    hasElementResolvedOrbitals,
    orbitalElementOptions,
  ]);

  const orbitalProjectionGroups = useMemo(() => {
    if (projectionMode !== "orbital") return [];
    if (!hasElementResolvedOrbitals) {
      const globalOrbitalGroups = data.projections?.orbital_groups ?? [];
      const totalGroup = buildTotalProjectionGroup(
        globalOrbitalGroups,
        "orbital-total",
        "Total (all orbitals)",
        data.n_bands,
        data.n_kpoints,
      );
      return totalGroup ? [totalGroup, ...globalOrbitalGroups] : globalOrbitalGroups;
    }
    const selectedElement =
      orbitalElementOptions.find((option) => option.key === selectedOrbitalElementKey) ??
      orbitalElementOptions[0];
    const selectedGroups = selectedElement?.groups ?? [];
    if (!selectedElement || selectedGroups.length === 0) {
      return selectedGroups;
    }
    const totalGroup = buildTotalProjectionGroup(
      selectedGroups,
      `element-orbital-total-${selectedElement.key}`,
      `${selectedElement.label} total (all orbitals)`,
      data.n_bands,
      data.n_kpoints,
      "element_orbital",
    );
    return totalGroup ? [totalGroup, ...selectedGroups] : selectedGroups;
  }, [
    data.n_bands,
    data.n_kpoints,
    data.projections,
    hasElementResolvedOrbitals,
    orbitalElementOptions,
    projectionMode,
    selectedOrbitalElementKey,
  ]);

  const projectionGroups = useMemo(() => {
    if (!hasProjectionData) return [];
    if (projectionMode === "atom") {
      return allElementProjectionGroups;
    }
    return orbitalProjectionGroups;
  }, [
    allElementProjectionGroups,
    hasProjectionData,
    orbitalProjectionGroups,
    projectionMode,
  ]);

  useEffect(() => {
    if (projectionMode !== "orbital") return;
    if (!hasElementResolvedOrbitals) {
      setSelectedOrbitalElementKey("");
      return;
    }
    const stillValid = orbitalElementOptions.some(
      (element) => element.key === selectedOrbitalElementKey,
    );
    if (!stillValid) {
      setSelectedOrbitalElementKey(orbitalElementOptions[0].key);
    }
  }, [
    hasElementResolvedOrbitals,
    orbitalElementOptions,
    projectionMode,
    selectedOrbitalElementKey,
  ]);

  useEffect(() => {
    if (projectionGroups.length === 0) {
      setSelectedProjectionId("");
      return;
    }
    const stillValid = projectionGroups.some((group) => group.id === selectedProjectionId);
    if (!stillValid) {
      setSelectedProjectionId(projectionGroups[0].id);
    }
  }, [projectionGroups, selectedProjectionId]);

  useEffect(() => {
    if (!hasProjectionData) {
      setFatBandsEnabled(false);
    }
  }, [hasProjectionData]);

  useEffect(() => {
    if (projectionSelection === undefined) {
      return;
    }

    if (!hasProjectionData || projectionSelection == null || projectionSelection === "none") {
      setFatBandsEnabled(false);
      return;
    }

    const selectedEntry = projectionSelectionEntries.find(
      (entry) => entry.value === projectionSelection,
    );
    if (!selectedEntry) {
      setFatBandsEnabled(false);
      return;
    }

    setFatBandsEnabled(true);
    setProjectionMode(selectedEntry.mode);
    if (selectedEntry.mode === "orbital") {
      setSelectedOrbitalElementKey(selectedEntry.elementKey ?? "");
    } else {
      setSelectedOrbitalElementKey("");
    }
    setSelectedProjectionId(selectedEntry.value);
  }, [hasProjectionData, projectionSelection, projectionSelectionEntries]);

  const selectedProjectionGroup = useMemo(() => {
    if (!selectedProjectionId) return null;
    return projectionGroups.find((group) => group.id === selectedProjectionId) || null;
  }, [projectionGroups, selectedProjectionId]);

  const selectedProjectionEntryValue = useMemo(() => {
    const exactMatch = projectionSelectionEntries.find((entry) =>
      entry.value === selectedProjectionId &&
      entry.mode === projectionMode &&
      (entry.mode !== "orbital" || (entry.elementKey ?? "") === selectedOrbitalElementKey)
    );
    if (exactMatch) {
      return exactMatch.value;
    }
    return projectionSelectionEntries.some((entry) => entry.value === selectedProjectionId)
      ? selectedProjectionId
      : projectionSelectionEntries[0]?.value ?? "";
  }, [
    projectionMode,
    projectionSelectionEntries,
    selectedOrbitalElementKey,
    selectedProjectionId,
  ]);

  const applyProjectionSelectionEntry = useCallback(
    (value: string) => {
      const selectedEntry = projectionSelectionEntries.find((entry) => entry.value === value);
      if (!selectedEntry) {
        return;
      }
      setFatBandsEnabled(true);
      setProjectionMode(selectedEntry.mode);
      if (selectedEntry.mode === "orbital") {
        setSelectedOrbitalElementKey(selectedEntry.elementKey ?? "");
      } else {
        setSelectedOrbitalElementKey("");
      }
      setSelectedProjectionId(selectedEntry.value);
    },
    [projectionSelectionEntries],
  );

  const projectionNormalizationGroups = useMemo(() => {
    if (!selectedProjectionGroup) return [];

    if (projectionMode === "atom") {
      return allElementProjectionGroups;
    }

    if (hasElementResolvedOrbitals) {
      if (selectedProjectionGroup.id.startsWith("element-orbital-total-")) {
        return allElementOrbitalTotalGroups;
      }
      return data.projections?.element_orbital_groups ?? [];
    }

    if (selectedProjectionGroup.id === "orbital-total") {
      return projectionGroups.filter((group) => group.id === "orbital-total");
    }

    return data.projections?.orbital_groups ?? [];
  }, [
    allElementOrbitalTotalGroups,
    allElementProjectionGroups,
    data.projections,
    hasElementResolvedOrbitals,
    projectionGroups,
    projectionMode,
    selectedProjectionGroup,
  ]);

  const normalizedProjectionWeights = useMemo(() => {
    const weights = selectedProjectionGroup?.weights;
    if (!weights || weights.length === 0) return null;

    const normalizationGroups = projectionNormalizationGroups;
    if (normalizationGroups.length === 0) return null;

    if (projectionNormalizeMode === "band") {
      const perBandDenominators = weights.map((_, bandIndex) => {
        let maxBandWeight = 0;
        for (const group of normalizationGroups) {
          const bandWeights = group.weights[bandIndex];
          if (!bandWeights) continue;
          for (const value of bandWeights) {
            if (Number.isFinite(value) && value > maxBandWeight) {
              maxBandWeight = value;
            }
          }
        }
        return maxBandWeight > 0 ? maxBandWeight : 1;
      });
      return weights.map((bandWeights, bandIndex) => {
        const denom = perBandDenominators[bandIndex] ?? 1;
        return bandWeights.map((value) => clamp01((Number.isFinite(value) ? value : 0) / denom));
      });
    }

    let globalMax = 0;
    for (const group of normalizationGroups) {
      for (const bandWeights of group.weights) {
        for (const value of bandWeights) {
          if (Number.isFinite(value) && value > globalMax) {
            globalMax = value;
          }
        }
      }
    }
    const denom = globalMax > 0 ? globalMax : 1;
    return weights.map((bandWeights) =>
      bandWeights.map((value) => clamp01((Number.isFinite(value) ? value : 0) / denom)),
    );
  }, [projectionNormalizationGroups, selectedProjectionGroup, projectionNormalizeMode]);

  const projectionColorRangeValid =
    Number.isFinite(projectionColorMin) &&
    Number.isFinite(projectionColorMax) &&
    projectionColorMax > projectionColorMin;

  const projectionActive =
    viewerType === "electronic" &&
    fatBandsEnabled &&
    selectedProjectionGroup !== null;

  const fatBandsActive =
    projectionActive &&
    projectionRenderMode === "fat" &&
    normalizedProjectionWeights !== null;

  const colorProjectionActive =
    projectionActive &&
    projectionRenderMode === "color" &&
    projectionColorRangeValid;

  // Generate SVG path for a band
  const bandToPath = useCallback(
    (band: number[], kPoints: number[]) => {
      if (band.length === 0 || kPoints.length === 0) return "";

      let path = `M ${scales.xScale(kPoints[0])},${scales.yScale(band[0])}`;
      for (let i = 1; i < band.length && i < kPoints.length; i++) {
        const command = kPoints[i] <= kPoints[i - 1] + 1e-10 ? "M" : "L";
        // Restart the polyline at k-path discontinuities to avoid non-physical spikes.
        path += ` ${command} ${scales.xScale(kPoints[i])},${scales.yScale(band[i])}`;
      }
      return path;
    },
    [scales],
  );

  const fatPoints = useMemo(() => {
    if (!fatBandsActive || !normalizedProjectionWeights) return [];

    const points: {
      key: string;
      bandIdx: number;
      cx: number;
      cy: number;
      r: number;
      fill: string;
      opacity: number;
    }[] = [];

    for (let bandIdx = 0; bandIdx < shiftedEnergies.length; bandIdx++) {
      const band = shiftedEnergies[bandIdx];
      const weights = normalizedProjectionWeights[bandIdx] || [];
      for (let kIdx = 0; kIdx < band.length && kIdx < data.k_points.length; kIdx++) {
        const normalizedWeight = weights[kIdx] || 0;
        if (normalizedWeight < 0.01) continue;

        const energy = band[kIdx];
        if (energy < scales.eMin - 1 || energy > scales.eMax + 1) continue;

        const radius = Math.min(16, 0.35 + Math.sqrt(normalizedWeight) * fatScale);
        const fill =
          fatColorMode === "accent"
            ? fatAccentColor
            : primaryBandStrokeColors[bandIdx] || fatAccentColor;
        const opacity = Math.max(0.06, Math.min(1, fatOpacity * (0.3 + normalizedWeight * 0.8)));

        points.push({
          key: `fat-${bandIdx}-${kIdx}`,
          bandIdx,
          cx: scales.xScale(data.k_points[kIdx]),
          cy: scales.yScale(energy),
          r: radius,
          fill,
          opacity,
        });
      }
    }

    return points;
  }, [
    fatAccentColor,
    fatBandsActive,
    fatColorMode,
    fatOpacity,
    fatScale,
    data.k_points,
    normalizedProjectionWeights,
    primaryBandStrokeColors,
    scales,
    shiftedEnergies,
  ]);

  const colorProjectionSegments = useMemo(() => {
    if (!colorProjectionActive || !selectedProjectionGroup) return [];

    const segments: {
      key: string;
      bandIdx: number;
      x1: number;
      y1: number;
      x2: number;
      y2: number;
      stroke: string;
    }[] = [];

    for (let bandIdx = 0; bandIdx < shiftedEnergies.length; bandIdx += 1) {
      const band = shiftedEnergies[bandIdx];
      const weights = selectedProjectionGroup.weights[bandIdx] || [];
      const pointCount = Math.min(band.length, data.k_points.length, weights.length);
      for (let kIdx = 1; kIdx < pointCount; kIdx += 1) {
        if (data.k_points[kIdx] <= data.k_points[kIdx - 1] + 1e-10) {
          continue;
        }

        const energy1 = band[kIdx - 1];
        const energy2 = band[kIdx];
        if (
          (energy1 < scales.eMin - 1 && energy2 < scales.eMin - 1) ||
          (energy1 > scales.eMax + 1 && energy2 > scales.eMax + 1)
        ) {
          continue;
        }

        const weight1 = Number.isFinite(weights[kIdx - 1]) ? weights[kIdx - 1] : 0;
        const weight2 = Number.isFinite(weights[kIdx]) ? weights[kIdx] : 0;
        const weight = (weight1 + weight2) / 2;

        segments.push({
          key: `color-projection-${bandIdx}-${kIdx}`,
          bandIdx,
          x1: scales.xScale(data.k_points[kIdx - 1]),
          y1: scales.yScale(energy1),
          x2: scales.xScale(data.k_points[kIdx]),
          y2: scales.yScale(energy2),
          stroke: projectionWeightColor(weight, projectionColorMin, projectionColorMax),
        });
      }
    }

    return segments;
  }, [
    colorProjectionActive,
    data.k_points,
    projectionColorMax,
    projectionColorMin,
    scales,
    selectedProjectionGroup,
    shiftedEnergies,
  ]);

  const projectionLegendTicks = useMemo(
    () => buildProjectionLegendTicks(
      projectionColorMin,
      projectionColorMax,
      projectionLegendTickCount,
    ),
    [projectionColorMax, projectionColorMin, projectionLegendTickCount],
  );

  const projectionLegendDecimals = useMemo(() => {
    const span = projectionColorMax - projectionColorMin;
    const tickStep = projectionLegendTicks.length > 1
      ? Math.abs(projectionLegendTicks[0] - projectionLegendTicks[1])
      : span;
    return Math.max(2, getTickDecimals(tickStep));
  }, [projectionColorMax, projectionColorMin, projectionLegendTicks]);

  // Handle mouse move for hover
  const handleMouseMove = useCallback(
    (e: React.MouseEvent<SVGSVGElement>) => {
      if (!svgRef.current) return;

      const rect = svgRef.current.getBoundingClientRect();
      const x = e.clientX - rect.left - margin.left;
      const y = e.clientY - rect.top - margin.top;

      // Find closest point
      if (x < 0 || x > plotWidth || y < 0 || y > plotHeight) {
        setHoveredPoint(null);
        return;
      }

      // Convert to data coordinates
      const k = scales.kMin + (x / plotWidth) * (scales.kMax - scales.kMin);
      const energy = scales.eMin + (1 - y / plotHeight) * (scales.eMax - scales.eMin);

      // Find closest k-point index
      let closestKIdx = 0;
      let minDist = Infinity;
      for (let i = 0; i < data.k_points.length; i++) {
        const dist = Math.abs(data.k_points[i] - k);
        if (dist < minDist) {
          minDist = dist;
          closestKIdx = i;
        }
      }

      // Find closest band at that k-point (using shifted energies)
      let closestBandIdx = 0;
      minDist = Infinity;
      for (let b = 0; b < shiftedEnergies.length; b++) {
        const dist = Math.abs(shiftedEnergies[b][closestKIdx] - energy);
        if (dist < minDist) {
          minDist = dist;
          closestBandIdx = b;
        }
      }

      // Only show tooltip if close enough
      const bandEnergy = shiftedEnergies[closestBandIdx][closestKIdx];
      const pixelDist = Math.abs(scales.yScale(bandEnergy) - y);
      if (pixelDist < 15) {
        const projectionWeight =
          selectedProjectionGroup?.weights?.[closestBandIdx]?.[closestKIdx];
        const projectionWeightNormalized =
          normalizedProjectionWeights?.[closestBandIdx]?.[closestKIdx];

        setHoveredPoint({
          band: closestBandIdx + 1,
          k: data.k_points[closestKIdx],
          energy: bandEnergy,
          x: scales.xScale(data.k_points[closestKIdx]),
          y: scales.yScale(bandEnergy),
          projectionWeight: Number.isFinite(projectionWeight ?? NaN)
            ? projectionWeight
            : undefined,
          projectionWeightNormalized: Number.isFinite(projectionWeightNormalized ?? NaN)
            ? projectionWeightNormalized
            : undefined,
        });
      } else {
        setHoveredPoint(null);
      }
    },
    [
      data.k_points,
      margin.left,
      margin.top,
      normalizedProjectionWeights,
      plotHeight,
      plotWidth,
      scales,
      selectedProjectionGroup,
      shiftedEnergies,
    ],
  );

  // Handle scroll to adjust Y-axis range
  const handleWheel = useCallback(
    (e: React.WheelEvent) => {
      if (!enableWheelRangeControl) {
        return;
      }
      e.preventDefault();

      const currentMin = yMin ?? scales.eMin;
      const currentMax = yMax ?? scales.eMax;
      const range = currentMax - currentMin;

      let newMin: number;
      let newMax: number;

      // Normal scroll zooms the Y range
      // Shift+scroll shifts the energy window (pan up/down)
      if (e.shiftKey) {
        const shift = e.deltaY > 0 ? range * 0.02 : -range * 0.02;
        newMin = currentMin + shift;
        newMax = currentMax + shift;
      } else {
        const factor = e.deltaY > 0 ? 1.03 : 0.97;
        const center = (currentMax + currentMin) / 2;
        const newRange = range * factor;
        newMin = center - newRange / 2;
        newMax = center + newRange / 2;
      }

      // Clamp to limits if configured (preserve window size)
      if (yClampRange) {
        const [minLimit, maxLimit] = yClampRange;
        const windowSize = newMax - newMin;

        // If panning too far down, clamp to lower limit
        if (newMin < minLimit) {
          newMin = minLimit;
          newMax = newMin + windowSize;
        }

        // If panning too far up, clamp to upper limit
        if (newMax > maxLimit) {
          newMax = maxLimit;
          newMin = newMax - windowSize;
        }
      }

      setYMin(newMin);
      setYMax(newMax);
      setManualRangeError(null);
    },
    [enableWheelRangeControl, yMin, yMax, scales, yClampRange],
  );

  const applyManualRange = useCallback(() => {
    const parsedMin = Number.parseFloat(manualYMinInput.trim());
    const parsedMax = Number.parseFloat(manualYMaxInput.trim());
    if (!Number.isFinite(parsedMin) || !Number.isFinite(parsedMax)) {
      setManualRangeError("Enter valid numeric values for Y min and Y max.");
      return;
    }
    if (parsedMax <= parsedMin) {
      setManualRangeError("Y max must be greater than Y min.");
      return;
    }

    setYMin(parsedMin);
    setYMax(parsedMax);
    setManualRangeError(null);
  }, [manualYMaxInput, manualYMinInput]);

  const handleManualRangeKeyDown = useCallback(
    (event: React.KeyboardEvent<HTMLInputElement>) => {
      if (event.key !== "Enter") return;
      event.preventDefault();
      applyManualRange();
    },
    [applyManualRange],
  );

  // Reset view
  const resetView = useCallback(() => {
    setYMin(null);
    setYMax(null);
    setManualRangeError(null);
  }, []);

  const handleExportPlaceholder = useCallback(() => {
    setExportNote("Export UI stub added. Wire export logic when ready.");
  }, []);

  const startValenceBandSelection = useCallback(() => {
    setIsSelectingValenceBand(true);
    setValenceBandSelectionStatus("Click the valence band on the plot.");
  }, []);

  const handleValenceBandSelection = useCallback(
    async (bandIndex: number) => {
      setSelectedValenceBandIndex(bandIndex);
      setIsSelectingValenceBand(false);
      setValenceBandSelectionStatus(`Valence band set to band ${bandIndex + 1}.`);

      if (!onPersistSelectedValenceBandIndex) {
        return;
      }

      try {
        await onPersistSelectedValenceBandIndex(bandIndex);
      } catch (error) {
        const message = error instanceof Error ? error.message : String(error);
        setValenceBandSelectionStatus(
          `Valence band set locally, but saving failed: ${message}`,
        );
      }
    },
    [onPersistSelectedValenceBandIndex],
  );

  const handlePlotClick = useCallback(
    (event: React.MouseEvent<SVGSVGElement>) => {
      if (!isSelectingValenceBand || !svgRef.current) {
        return;
      }

      const rect = svgRef.current.getBoundingClientRect();
      const x = event.clientX - rect.left - margin.left;
      const y = event.clientY - rect.top - margin.top;
      if (x < 0 || x > plotWidth || y < 0 || y > plotHeight) {
        return;
      }

      let closestBandIndex: number | null = null;
      let minDistanceSquared = Number.POSITIVE_INFINITY;

      for (let bandIdx = 0; bandIdx < shiftedEnergies.length; bandIdx += 1) {
        const band = shiftedEnergies[bandIdx];
        const pointCount = Math.min(band.length, data.k_points.length);
        if (pointCount === 0) {
          continue;
        }

        for (let pointIdx = 0; pointIdx < pointCount; pointIdx += 1) {
          const pointX = scales.xScale(data.k_points[pointIdx]);
          const pointY = scales.yScale(band[pointIdx]);

          if (pointIdx === 0 || data.k_points[pointIdx] <= data.k_points[pointIdx - 1] + 1e-10) {
            const pointDistX = x - pointX;
            const pointDistY = y - pointY;
            const pointDistanceSquared = pointDistX * pointDistX + pointDistY * pointDistY;
            if (pointDistanceSquared < minDistanceSquared) {
              minDistanceSquared = pointDistanceSquared;
              closestBandIndex = bandIdx;
            }
            continue;
          }

          const previousX = scales.xScale(data.k_points[pointIdx - 1]);
          const previousY = scales.yScale(band[pointIdx - 1]);
          const segmentDistanceSquared = distanceSquaredToSegment(
            x,
            y,
            previousX,
            previousY,
            pointX,
            pointY,
          );
          if (segmentDistanceSquared < minDistanceSquared) {
            minDistanceSquared = segmentDistanceSquared;
            closestBandIndex = bandIdx;
          }
        }
      }

      const selectionThresholdPx = 20;
      if (
        closestBandIndex == null
        || minDistanceSquared > selectionThresholdPx * selectionThresholdPx
      ) {
        setValenceBandSelectionStatus("Click closer to the target band.");
        return;
      }

      void handleValenceBandSelection(closestBandIndex);
    },
    [
      data.k_points,
      handleValenceBandSelection,
      isSelectingValenceBand,
      margin.left,
      margin.top,
      plotHeight,
      plotWidth,
      scales,
      shiftedEnergies,
    ],
  );

  // Prevent page scroll while interacting with the plot area.
  useEffect(() => {
    if (!enableHoverScrollLock) return;
    if (!isHoveringPlot || typeof document === "undefined") return;

    const { body, documentElement } = document;
    const prevBodyOverflow = body.style.overflow;
    const prevHtmlOverflow = documentElement.style.overflow;

    body.style.overflow = "hidden";
    documentElement.style.overflow = "hidden";

    return () => {
      body.style.overflow = prevBodyOverflow;
      documentElement.style.overflow = prevHtmlOverflow;
    };
  }, [enableHoverScrollLock, isHoveringPlot]);

  const drawBandLines = !colorProjectionActive && (!fatBandsActive || showLinesWithFat);
  const projectionLabel = selectedProjectionGroup?.label || "None";
  const projectionDisplayLabel =
    projectionSelectionEntries.find((entry) => entry.value === selectedProjectionEntryValue)
      ?.label ?? projectionLabel;
  const projectionLegendLabel = formatProjectionLegendLabel(projectionDisplayLabel);
  const showProjectionSummary = projectionActive && selectedProjectionGroup !== null;

  const svgBgFill = effectivePlotBgWhite ? "#ffffff" : colors.bg;

  return (
    <div className="band-plot-layout">
      <div className="band-plot-main">
        <div className="band-plot-canvas" ref={plotCanvasRef}>
          <svg
            ref={svgRef}
            width={resolvedWidth}
            height={resolvedHeight}
            onMouseEnter={() => setIsHoveringPlot(true)}
            onClick={handlePlotClick}
            onMouseMove={handleMouseMove}
            onMouseLeave={() => {
              setIsHoveringPlot(false);
              setHoveredPoint(null);
            }}
            onWheel={handleWheel}
            style={{ cursor: "crosshair", display: "block" }}
          >
            <defs>
              <clipPath id={clipPathId}>
                <rect x={0} y={0} width={plotWidth} height={plotHeight} />
              </clipPath>
              <linearGradient id={projectionGradientId} x1="0" y1="1" x2="0" y2="0">
                <stop offset="0%" stopColor={projectionScaleColor(0)} />
                <stop offset="100%" stopColor={projectionScaleColor(1)} />
              </linearGradient>
            </defs>

            {/* Background */}
            <rect width={resolvedWidth} height={resolvedHeight} fill={svgBgFill} />

            {/* Plot area */}
            <g transform={`translate(${margin.left}, ${margin.top})`}>
              {/* Grid lines */}
              <g className="grid-lines" opacity={0.3}>
                {yTicks.map((tick) => (
                  <line
                    key={tick}
                    x1={0}
                    x2={plotWidth}
                    y1={scales.yScale(tick)}
                    y2={scales.yScale(tick)}
                    stroke={colors.grid}
                    strokeDasharray="2,2"
                  />
                ))}
              </g>

              {overlayEditorEnabled && (
                <>
                  <g className="band-plot-window-overlays" clipPath={`url(#${clipPathId})`}>
                    {windowOverlayRenderData.map((entry) => {
                      const isDragging = overlayDragState?.overlayId === entry.overlay.id;
                      return (
                        <g key={`window-overlay-${entry.overlay.id}`}>
                          <rect
                            x={0}
                            y={entry.y}
                            width={plotWidth}
                            height={entry.heightPx}
                            fill={entry.overlay.color}
                            opacity={0.18}
                          />
                          <line
                            x1={0}
                            x2={plotWidth}
                            y1={entry.topY}
                            y2={entry.topY}
                            stroke={entry.overlay.color}
                            strokeWidth={1.5}
                            opacity={0.9}
                          />
                          <line
                            x1={0}
                            x2={plotWidth}
                            y1={entry.bottomY}
                            y2={entry.bottomY}
                            stroke={entry.overlay.color}
                            strokeWidth={1.5}
                            opacity={0.9}
                          />
                          <rect
                            x={0}
                            y={entry.y}
                            width={plotWidth}
                            height={entry.heightPx}
                            fill="transparent"
                            onPointerDown={(event) => beginOverlayDrag(event, entry.overlay, "move")}
                            style={{ cursor: isDragging ? "grabbing" : "grab", touchAction: "none" }}
                          />
                          <rect
                            x={0}
                            y={entry.topY - 5}
                            width={plotWidth}
                            height={10}
                            fill="transparent"
                            onPointerDown={(event) => beginOverlayDrag(event, entry.overlay, "resize-max")}
                            style={{ cursor: "ns-resize", touchAction: "none" }}
                          />
                          <rect
                            x={0}
                            y={entry.bottomY - 5}
                            width={plotWidth}
                            height={10}
                            fill="transparent"
                            onPointerDown={(event) => beginOverlayDrag(event, entry.overlay, "resize-min")}
                            style={{ cursor: "ns-resize", touchAction: "none" }}
                          />
                          <text
                            x={8}
                            y={Math.max(14, entry.y + 15)}
                            fill={colors.text}
                            fontSize={Math.max(10, 11 * effectivePlotTextScale)}
                            fontWeight={700}
                            pointerEvents="none"
                          >
                            {entry.overlay.label}
                          </text>
                        </g>
                      );
                    })}
                  </g>
                  <g className="band-plot-window-sliders">
                    {windowOverlayRenderData.map((entry) => {
                      const sliderBodyRawHeight = Math.abs(entry.bottomY - entry.topY);
                      const sliderBodyHeight = Math.max(12, sliderBodyRawHeight);
                      const sliderBodyY = Math.min(entry.topY, entry.bottomY) - (sliderBodyHeight - sliderBodyRawHeight) / 2;
                      return (
                        <g key={`window-slider-${entry.overlay.id}`} transform={`translate(${entry.sliderX}, 0)`}>
                        <line
                          x1={0}
                          x2={0}
                          y1={0}
                          y2={plotHeight}
                          stroke={entry.overlay.color}
                          strokeOpacity={0.3}
                          strokeWidth={2}
                        />
                        <line
                          x1={0}
                          x2={0}
                          y1={entry.topY}
                          y2={entry.bottomY}
                          stroke={entry.overlay.color}
                          strokeOpacity={0.85}
                          strokeWidth={5}
                        />
                        <rect
                          x={-8}
                          y={sliderBodyY}
                          width={16}
                          height={sliderBodyHeight}
                          fill="transparent"
                          onPointerDown={(event) => beginOverlayDrag(event, entry.overlay, "move")}
                          style={{
                            cursor: overlayDragState?.overlayId === entry.overlay.id ? "grabbing" : "grab",
                            touchAction: "none",
                          }}
                        />
                        <circle
                          cx={0}
                          cy={entry.topY}
                          r={6}
                          fill={entry.overlay.color}
                          stroke={svgBgFill}
                          strokeWidth={1.5}
                          onPointerDown={(event) => beginOverlayDrag(event, entry.overlay, "slider-max")}
                          style={{ cursor: "ns-resize", touchAction: "none" }}
                        />
                        <circle
                          cx={0}
                          cy={entry.bottomY}
                          r={6}
                          fill={entry.overlay.color}
                          stroke={svgBgFill}
                          strokeWidth={1.5}
                          onPointerDown={(event) => beginOverlayDrag(event, entry.overlay, "slider-min")}
                          style={{ cursor: "ns-resize", touchAction: "none" }}
                        />
                        </g>
                      );
                    })}
                  </g>
                </>
              )}

              {/* High-symmetry point vertical lines */}
              {data.high_symmetry_points.map((point, i) => (
                <g key={i}>
                  <line
                    x1={scales.xScale(point.k_distance)}
                    x2={scales.xScale(point.k_distance)}
                    y1={0}
                    y2={plotHeight}
                    stroke={colors.axis}
                    strokeWidth={0.5}
                  />
                  <text
                    x={scales.xScale(point.k_distance)}
                    y={plotHeight + symmetryLabelYOffset}
                    textAnchor="middle"
                    fill={colors.text}
                    fontSize={symmetryLabelFontSize}
                    fontFamily="serif"
                    fontStyle="italic"
                  >
                    {formatLabel(point.label)}
                  </text>
                </g>
              ))}

              {/* Zero-energy reference line */}
              {showFermiLevel && (
                <g>
                  <line
                    x1={0}
                    x2={plotWidth}
                    y1={scales.yScale(0)}
                    y2={scales.yScale(0)}
                    stroke="#d32f2f"
                    strokeWidth={1}
                    strokeDasharray="4,4"
                  />
                  <text
                    x={plotWidth + fermiLabelXOffset}
                    y={scales.yScale(0) + fermiLabelYOffset}
                    fill="#d32f2f"
                    fontSize={fermiLabelFontSize}
                  >
                    {zeroLineLabel}
                  </text>
                </g>
              )}

              {viewerType === "electronic" && isSelectingValenceBand && (
                <g pointerEvents="none">
                  <rect
                    x={12}
                    y={12}
                    width={230}
                    height={34}
                    rx={8}
                    fill={effectivePlotBgWhite ? "rgba(255,255,255,0.94)" : colors.tooltip}
                    stroke="#ef6c00"
                    strokeOpacity={0.85}
                  />
                  <text
                    x={24}
                    y={34}
                    fill="#ef6c00"
                    fontSize={12}
                    fontWeight={600}
                  >
                    Click the valence band to set VBM zero
                  </text>
                </g>
              )}

              {/* Band lines + fat points */}
              <g clipPath={`url(#${clipPathId})`}>
                {bandGapOverlayVisible && bandGapOverlay && (
                  <g opacity={1} pointerEvents="none">
                    <rect
                      x={0}
                      y={bandGapOverlay.gapTopY}
                      width={plotWidth}
                      height={Math.max(1.5, bandGapOverlay.gapBottomY - bandGapOverlay.gapTopY)}
                      fill={bandGapOverlay.gapColor}
                      opacity={0.1}
                    />
                    <line
                      x1={0}
                      x2={plotWidth}
                      y1={bandGapOverlay.gapTopY}
                      y2={bandGapOverlay.gapTopY}
                      stroke={bandGapOverlay.gapColor}
                      strokeOpacity={0.35}
                      strokeDasharray="5,4"
                    />
                    <line
                      x1={0}
                      x2={plotWidth}
                      y1={bandGapOverlay.gapBottomY}
                      y2={bandGapOverlay.gapBottomY}
                      stroke={bandGapOverlay.gapColor}
                      strokeOpacity={0.35}
                      strokeDasharray="5,4"
                    />
                  </g>
                )}

                {drawBandLines &&
                  shiftedEnergies.map((band, bandIdx) => (
                    <g key={bandIdx}>
                      <path
                        d={bandToPath(band, data.k_points)}
                        fill="none"
                        stroke={primaryBandStrokeColors[bandIdx]}
                        strokeWidth={effectiveLineWidth}
                        opacity={effectiveLineOpacity}
                      />
                      {viewerType === "electronic" && isSelectingValenceBand && (
                        <path
                          d={bandToPath(band, data.k_points)}
                          fill="none"
                          stroke="transparent"
                          strokeWidth={12}
                          style={{ cursor: "crosshair", pointerEvents: "stroke" }}
                          onClick={() => {
                            void handleValenceBandSelection(bandIdx);
                          }}
                        />
                      )}
                    </g>
                  ))}

                {colorProjectionActive &&
                  colorProjectionSegments.map((segment) => (
                    <line
                      key={segment.key}
                      x1={segment.x1}
                      y1={segment.y1}
                      x2={segment.x2}
                      y2={segment.y2}
                      stroke={segment.stroke}
                      strokeWidth={effectiveLineWidth}
                      strokeLinecap="round"
                      opacity={effectiveLineOpacity}
                      style={isSelectingValenceBand ? { cursor: "crosshair" } : undefined}
                      onClick={
                        viewerType === "electronic" && isSelectingValenceBand
                          ? () => {
                            void handleValenceBandSelection(segment.bandIdx);
                          }
                          : undefined
                      }
                    />
                  ))}

                {drawBandLines &&
                  selectedComparison &&
                  comparisonShiftedEnergies.map((band, bandIdx) => (
                    <path
                      key={`comparison-${selectedComparison.id}-${bandIdx}`}
                      d={bandToPath(band, comparisonKPointsForPlot)}
                      fill="none"
                      stroke={comparisonBandStrokeColor}
                      strokeWidth={effectiveLineWidth}
                      opacity={Math.min(1, effectiveLineOpacity + 0.08)}
                    />
                  ))}

                {fatBandsActive &&
                  fatPoints.map((point) => (
                    <circle
                      key={point.key}
                      cx={point.cx}
                      cy={point.cy}
                      r={point.r}
                      fill={point.fill}
                      opacity={point.opacity}
                      style={isSelectingValenceBand ? { cursor: "crosshair" } : undefined}
                      onClick={
                        viewerType === "electronic" && isSelectingValenceBand
                          ? () => {
                            void handleValenceBandSelection(point.bandIdx);
                          }
                          : undefined
                      }
                    />
                  ))}

                {bandGapOverlayVisible && bandGapOverlay && (
                  <g pointerEvents="none">
                    <line
                      x1={bandGapOverlay.vbmX}
                      y1={bandGapOverlay.vbmY}
                      x2={bandGapOverlay.cbmX}
                      y2={bandGapOverlay.cbmY}
                      stroke={bandGapOverlay.gapColor}
                      strokeWidth={displayedBandGap?.is_direct ? 2.5 : 2}
                      strokeDasharray={displayedBandGap?.is_direct ? undefined : "7,5"}
                      strokeLinecap="round"
                      opacity={0.9}
                    />
                    <circle
                      cx={bandGapOverlay.vbmX}
                      cy={bandGapOverlay.vbmY}
                      r={4.5}
                      fill="#2e7d32"
                      stroke={svgBgFill}
                      strokeWidth={1.5}
                    />
                    <circle
                      cx={bandGapOverlay.cbmX}
                      cy={bandGapOverlay.cbmY}
                      r={4.5}
                      fill="#ef6c00"
                      stroke={svgBgFill}
                      strokeWidth={1.5}
                    />
                    <text
                      x={bandGapOverlay.vbmX}
                      y={Math.max(
                        bandGapOverlay.edgeLabelFontSize + 2,
                        bandGapOverlay.vbmY - 9,
                      )}
                      textAnchor="middle"
                      fill="#2e7d32"
                      fontSize={bandGapOverlay.edgeLabelFontSize}
                      fontWeight={600}
                      stroke={svgBgFill}
                      strokeWidth={3}
                      paintOrder="stroke"
                    >
                      VBM
                    </text>
                    <text
                      x={bandGapOverlay.cbmX}
                      y={Math.min(
                        plotHeight - 4,
                        bandGapOverlay.cbmY + bandGapOverlay.edgeLabelFontSize + 10,
                      )}
                      textAnchor="middle"
                      fill="#ef6c00"
                      fontSize={bandGapOverlay.edgeLabelFontSize}
                      fontWeight={600}
                      stroke={svgBgFill}
                      strokeWidth={3}
                      paintOrder="stroke"
                    >
                      CBM
                    </text>
                    <rect
                      x={bandGapOverlay.labelX}
                      y={bandGapOverlay.labelY}
                      width={bandGapOverlay.labelWidth}
                      height={bandGapOverlay.labelHeight}
                      rx={7}
                      fill={effectivePlotBgWhite ? "rgba(255,255,255,0.96)" : colors.tooltip}
                      stroke={bandGapOverlay.gapColor}
                      strokeOpacity={0.8}
                    />
                    <text
                      x={bandGapOverlay.labelX + bandGapOverlay.labelWidth / 2}
                      y={bandGapOverlay.labelY + bandGapOverlay.labelFontSize + 4}
                      textAnchor="middle"
                      fill={colors.text}
                      fontSize={bandGapOverlay.labelFontSize}
                      fontWeight={700}
                    >
                      {bandGapOverlay.primaryLabel}
                    </text>
                    <text
                      x={bandGapOverlay.labelX + bandGapOverlay.labelWidth / 2}
                      y={
                        bandGapOverlay.labelY +
                        bandGapOverlay.labelFontSize +
                        bandGapOverlay.secondaryFontSize +
                        10
                      }
                      textAnchor="middle"
                      fill={bandGapOverlay.gapColor}
                      fontSize={bandGapOverlay.secondaryFontSize}
                      fontWeight={600}
                    >
                      {bandGapOverlay.secondaryLabel}
                    </text>
                  </g>
                )}
              </g>

              {/* Hover point */}
              {hoveredPoint && (
                <g>
                  <circle
                    cx={hoveredPoint.x}
                    cy={hoveredPoint.y}
                    r={5}
                    fill="#ff9800"
                    stroke="#fff"
                    strokeWidth={2}
                  />
                </g>
              )}

              {/* Y-axis */}
              <g>
                <line x1={0} y1={0} x2={0} y2={plotHeight} stroke={colors.axis} />
                {yTicks.map((tick) => (
                  <g key={tick}>
                    <line
                      x1={yTickMarkX}
                      x2={0}
                      y1={scales.yScale(tick)}
                      y2={scales.yScale(tick)}
                      stroke={colors.axis}
                    />
                    <text
                      x={yTickLabelX}
                      y={scales.yScale(tick) + yTickLabelYOffset}
                      textAnchor="end"
                      fill={colors.text}
                      fontSize={axisTickFontSize}
                    >
                      {tick.toFixed(yTickDecimals)}
                    </text>
                  </g>
                ))}
                <text
                  x={-yAxisLabelXOffset}
                  y={plotHeight / 2}
                  transform={`rotate(-90 ${-yAxisLabelXOffset} ${plotHeight / 2})`}
                  textAnchor="middle"
                  dominantBaseline="middle"
                  fill={colors.text}
                  fontSize={axisLabelFontSize}
                  fontFamily="serif"
                  fontStyle="italic"
                >
                  {getElectronicReferenceLabelSubscript(resolvedYAxisLabel) ? (
                    <>
                      <tspan>E − E</tspan>
                      <tspan
                        baselineShift="sub"
                        fontSize={Math.max(8, axisLabelFontSize * 0.72)}
                      >
                        {getElectronicReferenceLabelSubscript(resolvedYAxisLabel)}
                      </tspan>
                      <tspan baselineShift="baseline"> (eV)</tspan>
                    </>
                  ) : (
                    resolvedYAxisLabel
                  )}
                </text>
              </g>

              {/* X-axis line */}
              <line x1={0} y1={plotHeight} x2={plotWidth} y2={plotHeight} stroke={colors.axis} />

              {colorProjectionActive && (
                <g className="projection-color-key" pointerEvents="none">
                  {(() => {
                    const legendX = plotWidth + Math.max(28, axisTickFontSize * 2.2);
                    const legendWidth = Math.max(10, axisTickFontSize * 0.95);
                    const tickLength = Math.max(4, axisTickFontSize * 0.45);
                    const labelX = legendX + legendWidth + Math.max(36, axisTickFontSize * 3.2);
                    return (
                      <>
                        <rect
                          x={legendX}
                          y={0}
                          width={legendWidth}
                          height={plotHeight}
                          fill={`url(#${projectionGradientId})`}
                          stroke={colors.axis}
                          strokeWidth={0.5}
                        />
                        {projectionLegendTicks.map((tick) => {
                          const tickY = plotHeight
                            - ((tick - projectionColorMin) /
                              (projectionColorMax - projectionColorMin)) * plotHeight;
                          return (
                            <g key={tick}>
                              <line
                                x1={legendX + legendWidth}
                                x2={legendX + legendWidth + tickLength}
                                y1={tickY}
                                y2={tickY}
                                stroke={colors.axis}
                                strokeWidth={0.75}
                              />
                              <text
                                x={legendX + legendWidth + tickLength + 3}
                                y={tickY + yTickLabelYOffset}
                                fill={colors.text}
                                fontSize={Math.max(7, axisTickFontSize * 0.78)}
                              >
                                {tick.toFixed(projectionLegendDecimals)}
                              </text>
                            </g>
                          );
                        })}
                        <text
                          x={labelX}
                          y={plotHeight / 2}
                          transform={`rotate(-90 ${labelX} ${plotHeight / 2})`}
                          textAnchor="middle"
                          dominantBaseline="middle"
                          fill={colors.text}
                          fontSize={Math.max(8, axisTickFontSize * 0.95)}
                          fontWeight={600}
                        >
                          {projectionLegendLabel}
                        </text>
                      </>
                    );
                  })()}
                </g>
              )}
            </g>

            {/* Tooltip */}
            {hoveredPoint && (
              <g
                transform={`translate(${margin.left + hoveredPoint.x + 10}, ${
                  margin.top + hoveredPoint.y - 10
                })`}
              >
                <rect
                  x={0}
                  y={-42}
                  width={190}
                  height={
                    hoveredPoint.projectionWeight !== undefined &&
                    hoveredPoint.projectionWeightNormalized !== undefined
                      ? 66
                      : 48
                  }
                  fill={colors.tooltip}
                  stroke={colors.tooltipBorder}
                  rx={4}
                  filter="drop-shadow(0 2px 4px rgba(0,0,0,0.1))"
                />
                <text x={8} y={-24} fill={colors.tooltipText} fontSize={11}>
                  {pointLabel} {hoveredPoint.band}
                </text>
                <text x={8} y={-8} fill="#1565c0" fontSize={11}>
                  {resolvedValueLabel} = {hoveredPoint.energy.toFixed(valueDecimals)} {valueUnit}
                </text>
                {hoveredPoint.projectionWeight !== undefined &&
                  hoveredPoint.projectionWeightNormalized !== undefined && (
                    <text x={8} y={10} fill="#6d4c41" fontSize={11}>
                      Weight = {hoveredPoint.projectionWeight.toExponential(2)} (
                      {(hoveredPoint.projectionWeightNormalized * 100).toFixed(1)}%)
                    </text>
                  )}
              </g>
            )}
          </svg>
        </div>

        {/* Info panel */}
        <div className="band-plot-info">
          <span>
            {data.n_bands} {primaryCountLabel}
          </span>
          <span>
            {data.n_kpoints} {secondaryCountLabel}
          </span>
          {showFermiLevel && (
            <span>
              E_F ({activeFermiSourceLabel}) = {activeFermiDisplay}
            </span>
          )}
          {viewerType === "electronic" && (
            <span>
              Zero ({activeZeroReferenceLabel}) = {activeZeroReferenceDisplay}
            </span>
          )}
          {viewerType === "electronic" && (requestedZeroEnergyReferenceMode ?? zeroEnergyReferenceMode) === "vbm" && (
            <span>
              Valence band: {selectedValenceBandLabel}
            </span>
          )}
          {viewerType === "electronic" && (
            <span className={displayedBandGap ? "band-gap-info" : "metal-info"}>
              {displayedBandGap
                ? `${displayedBandGap.value.toFixed(3)} eV ${
                  displayedBandGap.is_direct ? "direct" : "indirect"
                } gap`
                : (requestedZeroEnergyReferenceMode ?? zeroEnergyReferenceMode) === "vbm" && selectedValenceBandIndex == null
                  ? "Select a valence band to compute the VBM-referenced gap"
                : `No gap detected for current ${activeFermiSourceLabel.toLowerCase()} Fermi reference`}
            </span>
          )}
          {showProjectionSummary && (
            <span className="band-plot-projection-pill">
              {projectionRenderMode === "color" ? "Color projection" : "Fat bands"}: {projectionDisplayLabel}
            </span>
          )}
          {selectedComparison && (
            <span className="band-plot-projection-pill">
              Overlay: {selectedComparison.label}
            </span>
          )}
          {overlayEditorEnabled && (
            <span className="band-plot-overlay-hint">{windowOverlayHint}</span>
          )}
        </div>
      </div>

      {showSidebar && (
        <div className="band-plot-sidebar">
        <div className="band-plot-controls">
          <button onClick={resetView} className="band-plot-reset">
            Reset View
          </button>
          {viewerType === "electronic" && (
            <button onClick={handleExportPlaceholder} className="band-plot-export">
              Export (stub)
            </button>
          )}
          <span className="band-plot-hint">{scrollHint}</span>
        </div>

        {exportNote && <div className="band-plot-export-note">{exportNote}</div>}

        <div className="band-plot-control-panel">
          {viewerType === "electronic" && hasComparisonControls && (
            <div className="band-control-section">
              <button
                type="button"
                className="band-control-header"
                onClick={() => setComparisonExpanded((prev) => !prev)}
              >
                <span className={`collapse-icon ${comparisonExpanded ? "expanded" : ""}`}>▶</span>
                {comparisonTitle}
              </button>
              {comparisonExpanded && (
                <div className="band-control-grid">
                  <div className="band-control-row">
                    <label>Overlay Bands</label>
                    <select
                      value={selectedComparisonId}
                      onChange={(event) => setSelectedComparisonId(event.target.value)}
                      disabled={availableComparisonOptions.length === 0}
                    >
                      <option value="">{availableComparisonOptions.length === 0 ? "No matching saved bands" : comparisonNoneLabel}</option>
                      {availableComparisonOptions.map((option) => (
                        <option key={option.id} value={option.id}>{option.label}</option>
                      ))}
                    </select>
                  </div>
                  <div className="band-control-note">
                    Comparison mode draws Wannier bands in blue and the selected saved band run in orange.
                  </div>
                </div>
              )}
            </div>
          )}

          <div className="band-control-section">
            <button
              type="button"
              className="band-control-header"
              onClick={() => setAppearanceExpanded((prev) => !prev)}
            >
              <span className={`collapse-icon ${appearanceExpanded ? "expanded" : ""}`}>▶</span>
              Appearance
            </button>
            {appearanceExpanded && (
              <div className="band-control-grid">
                {viewerType === "electronic" && (hasScfFermi || hasBandsFermi) && (
                  <div className="band-control-row">
                    <label>Fermi Reference</label>
                    <select
                      value={fermiReferenceMode}
                      onChange={(event) =>
                        setFermiReferenceMode(event.target.value as FermiReferenceMode)
                      }
                    >
                      {hasScfFermi && (
                        <option value="scf">
                          SCF ({(scfFermiEnergy as number).toFixed(3)} eV)
                        </option>
                      )}
                      {hasBandsFermi && (
                        <option value="bands">
                          Bands run ({data.fermi_energy.toFixed(3)} eV)
                        </option>
                      )}
                    </select>
                  </div>
                )}

                {viewerType === "electronic" && (
                  <>
                    <div className="band-control-row">
                      <label>Zero Energy</label>
                      <select
                        value={zeroEnergyReferenceMode}
                        onChange={(event) =>
                          setZeroEnergyReferenceMode(event.target.value as ZeroEnergyReferenceMode)
                        }
                      >
                        <option value="calculated-fermi">Calculated E_F</option>
                        <option value="vbm">VBM</option>
                      </select>
                    </div>
                    {(requestedZeroEnergyReferenceMode ?? zeroEnergyReferenceMode) === "vbm" && (
                      <div className="band-control-row">
                        <label>Valence Band</label>
                        <div className="band-control-stack">
                          <button
                            type="button"
                            className="band-control-apply"
                            onClick={startValenceBandSelection}
                          >
                            Change Valence Band
                          </button>
                          <div className="band-control-note">
                            Selected: {selectedValenceBandLabel}
                          </div>
                          {valenceBandSelectionStatus && (
                            <div className="band-control-vbm-status">
                              {valenceBandSelectionStatus}
                            </div>
                          )}
                        </div>
                      </div>
                    )}
                  </>
                )}

                <div className="band-control-row">
                  <label>Y Range ({valueUnit})</label>
                  <div className="band-control-range-inputs">
                    <input
                      type="number"
                      step="any"
                      value={manualYMinInput}
                      onChange={(event) => {
                        setManualYMinInput(event.target.value);
                        setManualRangeError(null);
                      }}
                      onKeyDown={handleManualRangeKeyDown}
                      aria-label="Y minimum"
                    />
                    <span className="band-control-range-separator">to</span>
                    <input
                      type="number"
                      step="any"
                      value={manualYMaxInput}
                      onChange={(event) => {
                        setManualYMaxInput(event.target.value);
                        setManualRangeError(null);
                      }}
                      onKeyDown={handleManualRangeKeyDown}
                      aria-label="Y maximum"
                    />
                  </div>
                  <button
                    type="button"
                    className="band-control-apply"
                    onClick={applyManualRange}
                  >
                    Apply
                  </button>
                  {manualRangeError && (
                    <span className="band-control-range-error">{manualRangeError}</span>
                  )}
                </div>

                <div className="band-control-row">
                  <label>Y Tick Step ({valueUnit})</label>
                  <input
                    type="number"
                    step="any"
                    min={0}
                    value={manualYTickStepInput}
                    onChange={(event) => setManualYTickStepInput(event.target.value)}
                    placeholder="Auto"
                    aria-label="Y tick step"
                  />
                  <div className="band-control-note">
                    Leave blank to use automatic tick spacing.
                  </div>
                </div>

                <div className="band-control-row">
                  <label>Line Thickness</label>
                  <input
                    type="range"
                    min={0.5}
                    max={5}
                    step={0.1}
                    value={lineWidth}
                    onChange={(event) => setLineWidth(Number.parseFloat(event.target.value))}
                  />
                  <span className="band-control-value">{lineWidth.toFixed(1)} px</span>
                </div>

                <div className="band-control-row">
                  <label>Line Opacity</label>
                  <input
                    type="range"
                    min={0.1}
                    max={1}
                    step={0.05}
                    value={lineOpacity}
                    onChange={(event) => setLineOpacity(Number.parseFloat(event.target.value))}
                  />
                  <span className="band-control-value">{lineOpacity.toFixed(2)}</span>
                </div>

                <div className="band-control-row">
                  <label>Plot Text Size</label>
                  <input
                    type="range"
                    min={0.7}
                    max={2}
                    step={0.05}
                    value={plotTextScale}
                    onChange={(event) => setPlotTextScale(Number.parseFloat(event.target.value))}
                  />
                  <span className="band-control-value">{plotTextScale.toFixed(2)}x</span>
                </div>

                <div className="band-control-row">
                  <label>Band Color Mode</label>
                  <select
                    value={colorMode}
                    onChange={(event) => setColorMode(event.target.value as ColorMode)}
                  >
                    <option value="single">Single color</option>
                    <option value="rainbow">Rainbow by band</option>
                  </select>
                </div>

                {colorMode === "single" && (
                  <div className="band-control-row">
                    <label>Band Color</label>
                    <input
                      type="color"
                      value={singleBandColor}
                      onChange={(event) => setSingleBandColor(event.target.value)}
                    />
                  </div>
                )}

                {colorMode === "rainbow" && (
                  <div className="band-control-row">
                    <label>Rainbow Palette</label>
                    <select
                      value={rainbowPalette}
                      onChange={(event) =>
                        setRainbowPalette(event.target.value as RainbowPalette)
                      }
                    >
                      <option value="jet">Jet-like</option>
                      <option value="sinebow">Sinebow</option>
                    </select>
                  </div>
                )}

                <div className="band-control-row">
                  <label>Plot Background</label>
                  <select
                    value={plotBgWhite ? "white" : "theme"}
                    onChange={(event) => setPlotBgWhite(event.target.value === "white")}
                  >
                    <option value="white">Always White</option>
                    <option value="theme">Match Theme</option>
                  </select>
                </div>
              </div>
            )}
          </div>

          {viewerType === "electronic" && (
            <div className="band-control-section">
              <button
                type="button"
                className="band-control-header"
                onClick={() => setBandGapExpanded((prev) => !prev)}
              >
                <span className={`collapse-icon ${bandGapExpanded ? "expanded" : ""}`}>▶</span>
                Band Gap
              </button>
              {bandGapExpanded && (
                <div className="band-control-grid">
                  <div className="band-control-row">
                    <label>Show Gap Overlay</label>
                    <input
                      type="checkbox"
                      checked={bandGapOverlayVisible}
                      disabled={!displayedBandGap}
                      onChange={(event) => setShowBandGapOverlay(event.target.checked)}
                    />
                  </div>

                  {displayedBandGap ? (
                    <>
                      <div className="band-control-readout">
                        <span>Gap Value</span>
                        <strong>{displayedBandGap.value.toFixed(3)} eV</strong>
                      </div>
                      <div className="band-control-readout">
                        <span>Gap Type</span>
                        <strong>{displayedBandGap.is_direct ? "Direct" : "Indirect"}</strong>
                      </div>
                      <div className="band-control-note">
                        Referenced to the current {activeZeroReferenceLabel} zero.
                      </div>
                    </>
                  ) : (
                    <div className="band-control-warning">
                      No band gap was detected for the current {activeFermiSourceLabel.toLowerCase()} Fermi reference.
                    </div>
                  )}
                </div>
              )}
            </div>
          )}

          {viewerType === "electronic" && (
            <div className="band-control-section">
              <button
                type="button"
                className="band-control-header"
                onClick={() => setProjectionExpanded((prev) => !prev)}
              >
                <span className={`collapse-icon ${projectionExpanded ? "expanded" : ""}`}>▶</span>
                Fat Bands & Projections
              </button>
              {projectionExpanded && (
                <div className="band-control-grid">
                  <div className="band-control-row">
                    <label>
                      Enable Projection
                      <span className="band-control-tech-name">projwfc.x</span>
                    </label>
                    <input
                      type="checkbox"
                      checked={fatBandsEnabled}
                      disabled={!hasProjectionData}
                      onChange={(event) => setFatBandsEnabled(event.target.checked)}
                    />
                  </div>

                  {!hasProjectionData && (
                    <div className="band-control-warning">
                      Projection data not found. Re-run with orbital projections enabled.
                    </div>
                  )}

                  {hasProjectionData && (
                    <>
                      <div className="band-control-row">
                        <label>Projection Style</label>
                        <select
                          value={projectionRenderMode}
                          onChange={(event) =>
                            setProjectionRenderMode(event.target.value as ProjectionRenderMode)
                          }
                        >
                          <option value="fat">Fat bands</option>
                          <option value="color">Red-blue color bands</option>
                        </select>
                      </div>

                      <div className="band-control-row">
                        <label>Projection</label>
                        <select
                          value={selectedProjectionEntryValue}
                          onChange={(event) => applyProjectionSelectionEntry(event.target.value)}
                          disabled={projectionSelectionEntries.length === 0}
                        >
                          {projectionSelectionEntries.map((entry) => (
                            <option key={entry.value} value={entry.value}>
                              {entry.label}
                            </option>
                          ))}
                        </select>
                      </div>

                      {projectionRenderMode === "fat" && (
                        <>
                          <div className="band-control-row">
                            <label>Weight Normalization</label>
                            <select
                              value={projectionNormalizeMode}
                              onChange={(event) =>
                                setProjectionNormalizeMode(
                                  event.target.value as ProjectionNormalizeMode,
                                )
                              }
                            >
                              <option value="global">
                                {projectionMode === "atom"
                                  ? "Global (all elements)"
                                  : hasElementResolvedOrbitals
                                    ? selectedProjectionId.startsWith("element-orbital-total-")
                                      ? "Global (all element totals)"
                                      : "Global (all element orbitals)"
                                    : "Global (all orbitals)"}
                              </option>
                              <option value="band">
                                {projectionMode === "atom"
                                  ? "Per band (all elements)"
                                  : hasElementResolvedOrbitals
                                    ? selectedProjectionId.startsWith("element-orbital-total-")
                                      ? "Per band (all element totals)"
                                      : "Per band (all element orbitals)"
                                    : "Per band (all orbitals)"}
                              </option>
                            </select>
                          </div>

                          <div className="band-control-row">
                            <label>Fat-Band Scale</label>
                            <input
                              type="range"
                              min={2}
                              max={20}
                              step={0.5}
                              value={fatScale}
                              onChange={(event) => setFatScale(Number.parseFloat(event.target.value))}
                            />
                            <span className="band-control-value">{fatScale.toFixed(1)}</span>
                          </div>

                          <div className="band-control-row">
                            <label>Fat-Band Opacity</label>
                            <input
                              type="range"
                              min={0.05}
                              max={1}
                              step={0.05}
                              value={fatOpacity}
                              onChange={(event) =>
                                setFatOpacity(Number.parseFloat(event.target.value))
                              }
                            />
                            <span className="band-control-value">{fatOpacity.toFixed(2)}</span>
                          </div>

                          <div className="band-control-row">
                            <label>Fat-Band Color Mode</label>
                            <select
                              value={fatColorMode}
                              onChange={(event) => setFatColorMode(event.target.value as FatColorMode)}
                            >
                              <option value="band">Match band color</option>
                              <option value="accent">Single accent color</option>
                            </select>
                          </div>
                        </>
                      )}

                      {projectionRenderMode === "fat" && fatColorMode === "accent" && (
                        <div className="band-control-row">
                          <label>Accent Color</label>
                          <input
                            type="color"
                            value={fatAccentColor}
                            onChange={(event) => setFatAccentColor(event.target.value)}
                          />
                        </div>
                      )}

                      {projectionRenderMode === "fat" && (
                        <div className="band-control-row">
                          <label>Keep Band Lines Visible</label>
                          <input
                            type="checkbox"
                            checked={showLinesWithFat}
                            onChange={(event) => setShowLinesWithFat(event.target.checked)}
                          />
                        </div>
                      )}

                      {projectionRenderMode === "color" && (
                        <>
                          <div className="band-control-row">
                            <label>Color Range</label>
                            <div className="band-control-range-inputs">
                              <input
                                type="number"
                                step="any"
                                value={projectionColorMin}
                                onChange={(event) => {
                                  const parsed = Number.parseFloat(event.target.value);
                                  if (Number.isFinite(parsed)) {
                                    setProjectionColorMin(parsed);
                                  }
                                }}
                                aria-label="Projection blue minimum"
                              />
                              <span className="band-control-range-separator">to</span>
                              <input
                                type="number"
                                step="any"
                                value={projectionColorMax}
                                onChange={(event) => {
                                  const parsed = Number.parseFloat(event.target.value);
                                  if (Number.isFinite(parsed)) {
                                    setProjectionColorMax(parsed);
                                  }
                                }}
                                aria-label="Projection red maximum"
                              />
                            </div>
                            <div className="band-control-note">
                              Raw projection weights: blue at the lower value, red at the upper value.
                            </div>
                            {!projectionColorRangeValid && (
                              <span className="band-control-range-error">
                                Upper projection value must be greater than lower value.
                              </span>
                            )}
                          </div>

                          <div className="band-control-row">
                            <label>Key Tick Count</label>
                            <input
                              type="number"
                              min={2}
                              max={12}
                              step={1}
                              value={projectionLegendTickCount}
                              onChange={(event) => {
                                const parsed = Number.parseInt(event.target.value, 10);
                                if (Number.isFinite(parsed)) {
                                  setProjectionLegendTickCount(
                                    Math.max(2, Math.min(12, parsed)),
                                  );
                                }
                              }}
                              aria-label="Projection key tick count"
                            />
                            <div className="band-control-note">
                              Includes the minimum and maximum labels.
                            </div>
                          </div>
                        </>
                      )}
                    </>
                  )}
                </div>
              )}
            </div>
          )}
        </div>
        </div>
      )}
    </div>
  );
}
