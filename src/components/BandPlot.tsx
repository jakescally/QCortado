import { useState, useRef, useMemo, useCallback, useEffect, useId, useLayoutEffect } from "react";
import { useTheme } from "../lib/ThemeContext";

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

interface BandPlotProps {
  data: BandData;
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

export type ColorMode = "single" | "rainbow";
export type RainbowPalette = "jet" | "sinebow";
type ProjectionMode = "atom" | "orbital";
type ProjectionNormalizeMode = "global" | "band";
type FatColorMode = "accent" | "band";
export type FermiReferenceMode = "scf" | "bands";

export interface BandPlotSharedSettings {
  fermiReferenceMode: FermiReferenceMode;
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

function getYAxisTicks(eMin: number, eMax: number): number[] {
  if (!Number.isFinite(eMin) || !Number.isFinite(eMax) || eMax < eMin) {
    return [0];
  }
  const step = getYAxisTickStep(eMax - eMin);
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

function isElectronicEFLabel(label: string): boolean {
  const normalized = label.replace(/\s+/g, "").replace(/−/g, "-").toLowerCase();
  return normalized === "e-e_f(ev)";
}

export function resolveBandPlotFermiContext(
  data: BandData,
  scfFermiEnergy?: number | null,
  requestedMode: FermiReferenceMode = "bands",
): {
  fermiEnergy: number;
  mode: FermiReferenceMode;
  hasScfFermi: boolean;
  hasBandsFermi: boolean;
} {
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
  data: BandData,
  scfFermiEnergy?: number | null,
  requestedMode: FermiReferenceMode = "bands",
): [number, number] {
  const { fermiEnergy } = resolveBandPlotFermiContext(data, scfFermiEnergy, requestedMode);
  let eMin = data.energy_range[0] - fermiEnergy;
  let eMax = data.energy_range[1] - fermiEnergy;
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

function calculateDisplayedBandGap(energies: number[][], kPoints: number[]): BandGap | null {
  if (energies.length === 0 || kPoints.length === 0) {
    return null;
  }

  let vbmEnergy = Number.NEGATIVE_INFINITY;
  let vbmK = 0;
  let cbmEnergy = Number.POSITIVE_INFINITY;
  let cbmK = 0;

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

      if (energy <= BAND_GAP_TOLERANCE_EV && energy > vbmEnergy) {
        vbmEnergy = energy;
        vbmK = kPoints[index];
      }

      if (energy >= -BAND_GAP_TOLERANCE_EV && energy < cbmEnergy) {
        cbmEnergy = energy;
        cbmK = kPoints[index];
      }
    }

    if (
      Number.isFinite(bandMin) &&
      Number.isFinite(bandMax) &&
      bandMin < -BAND_GAP_TOLERANCE_EV &&
      bandMax > BAND_GAP_TOLERANCE_EV
    ) {
      return null;
    }
  }

  if (!Number.isFinite(vbmEnergy) || !Number.isFinite(cbmEnergy)) {
    return null;
  }

  const gapValue = cbmEnergy - vbmEnergy;
  if (!(gapValue > BAND_GAP_TOLERANCE_EV)) {
    return null;
  }

  return {
    value: gapValue,
    is_direct: Math.abs(vbmK - cbmK) < DIRECT_GAP_K_TOLERANCE,
    vbm_k: vbmK,
    cbm_k: cbmK,
    vbm_energy: vbmEnergy,
    cbm_energy: cbmEnergy,
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

export function buildBandPlotProjectionOptions(data: BandData): BandPlotProjectionOption[] {
  return [
    { value: "none", label: "none" },
    ...buildProjectionSelectionEntries(data).map(({ value, label }) => ({ value, label })),
  ];
}

export function BandPlot({
  data,
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
}: BandPlotProps) {
  const { isDark } = useTheme();
  const colors = useMemo(() => isDark
    ? { bg: "#1e1e2e", axis: "#718096", grid: "#4a5568", text: "#e2e8f0", tooltip: "#2d3748", tooltipBorder: "#4a5568", tooltipText: "#e2e8f0" }
    : { bg: "#ffffff", axis: "#333", grid: "#999", text: "#000", tooltip: "#fff", tooltipBorder: "#ccc", tooltipText: "#333" },
  [isDark]);

  const svgRef = useRef<SVGSVGElement>(null);
  const plotCanvasRef = useRef<HTMLDivElement>(null);
  const clipPathId = useId();
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

  // Plot background toggle
  const [plotBgWhite, setPlotBgWhite] = useState(true);
  const [showBandGapOverlay, setShowBandGapOverlay] = useState(true);

  // UI section toggles
  const [appearanceExpanded, setAppearanceExpanded] = useState(true);
  const [bandGapExpanded, setBandGapExpanded] = useState(true);
  const [projectionExpanded, setProjectionExpanded] = useState(false);
  const [exportNote, setExportNote] = useState("");

  const requestedFermiReferenceMode = sharedSettings?.fermiReferenceMode ?? null;
  const fallbackFermiReferenceMode = (
    scfFermiEnergy != null && Number.isFinite(scfFermiEnergy)
  ) ? "scf" : "bands";
  const [fermiReferenceMode, setFermiReferenceMode] = useState<FermiReferenceMode>(
    requestedFermiReferenceMode ?? fallbackFermiReferenceMode,
  );

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

  const {
    fermiEnergy,
    mode: resolvedFermiReferenceMode,
    hasScfFermi,
    hasBandsFermi,
  } = useMemo(
    () =>
      resolveBandPlotFermiContext(
        data,
        scfFermiEnergy,
        requestedFermiReferenceMode ?? fermiReferenceMode,
      ),
    [data, fermiReferenceMode, requestedFermiReferenceMode, scfFermiEnergy],
  );

  const effectiveLineWidth = sharedSettings?.lineWidth ?? lineWidth;
  const effectiveLineOpacity = sharedSettings?.lineOpacity ?? lineOpacity;
  const effectivePlotTextScale = sharedSettings?.plotTextScale ?? plotTextScale;
  const effectiveColorMode = sharedSettings?.colorMode ?? colorMode;
  const effectiveSingleBandColor = sharedSettings?.singleBandColor ?? singleBandColor;
  const effectiveRainbowPalette = sharedSettings?.rainbowPalette ?? rainbowPalette;
  const effectivePlotBgWhite = sharedSettings?.plotBgWhite ?? plotBgWhite;
  const effectiveShowBandGapOverlay = sharedSettings?.showBandGapOverlay ?? showBandGapOverlay;

  const activeFermiSourceLabel = resolvedFermiReferenceMode === "scf" ? "SCF" : "Bands run";
  const activeFermiDisplay = Number.isFinite(fermiEnergy) ? `${fermiEnergy.toFixed(3)} eV` : "N/A";

  // Shift all energies relative to Fermi level (E - E_F)
  const shiftedEnergies = useMemo(() => {
    return data.energies.map((band) => band.map((e) => e - fermiEnergy));
  }, [data.energies, fermiEnergy]);

  const displayedBandGap = useMemo(() => {
    if (viewerType !== "electronic") return null;
    return calculateDisplayedBandGap(shiftedEnergies, data.k_points);
  }, [data.k_points, shiftedEnergies, viewerType]);

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
    );
  }, [
    data,
    energyRange,
    fermiReferenceMode,
    requestedFermiReferenceMode,
    scfFermiEnergy,
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

  const axisTickFontSize = Math.max(8, 11 * effectivePlotTextScale);
  const axisLabelFontSize = Math.max(10, 14 * effectivePlotTextScale);
  const symmetryLabelFontSize = axisLabelFontSize;
  const symmetryLabelYOffset = Math.max(20, symmetryLabelFontSize * 1.35);
  const yTickLabelYOffset = axisTickFontSize * 0.35;
  const yTicks = useMemo(() => getYAxisTicks(yDomain[0], yDomain[1]), [yDomain]);
  const maxYTickLabelChars = useMemo(() => {
    if (yTicks.length === 0) return 3;
    let maxChars = 0;
    for (const tick of yTicks) {
      maxChars = Math.max(maxChars, tick.toFixed(1).length);
    }
    return maxChars;
  }, [yTicks]);
  const estimatedYTickLabelWidth = Math.max(
    axisTickFontSize * 3,
    maxYTickLabelChars * axisTickFontSize * 0.62,
  );
  const fermiLabelFontSize = axisTickFontSize;
  const fermiLabelYOffset = fermiLabelFontSize * 0.35;

  // Margins & dimensions — scale from actual axis text footprint.
  const margin = useMemo(() => ({
    top: 30,
    right: Math.round(Math.max(30, 16 + fermiLabelFontSize * 2.2)),
    bottom: Math.round(Math.max(50, 18 + symmetryLabelFontSize * 2.2)),
    left: Math.round(
      Math.max(70, 20 + estimatedYTickLabelWidth + axisLabelFontSize * 1.6),
    ),
  }), [axisLabelFontSize, estimatedYTickLabelWidth, fermiLabelFontSize, symmetryLabelFontSize]);

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

  const fatBandsActive =
    viewerType === "electronic" &&
    fatBandsEnabled &&
    normalizedProjectionWeights !== null &&
    selectedProjectionGroup !== null;

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
          fatColorMode === "accent" ? fatAccentColor : bandColors[bandIdx] || fatAccentColor;
        const opacity = Math.max(0.06, Math.min(1, fatOpacity * (0.3 + normalizedWeight * 0.8)));

        points.push({
          key: `fat-${bandIdx}-${kIdx}`,
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
    bandColors,
    data.k_points,
    normalizedProjectionWeights,
    scales,
    shiftedEnergies,
  ]);

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

  const drawBandLines = !fatBandsActive || showLinesWithFat;
  const projectionLabel = selectedProjectionGroup?.label || "None";
  const showProjectionSummary = fatBandsActive && selectedProjectionGroup !== null;

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

              {/* Fermi level at E - E_F = 0 */}
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
                    E_F
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
                    <path
                      key={bandIdx}
                      d={bandToPath(band, data.k_points)}
                      fill="none"
                      stroke={bandColors[bandIdx]}
                      strokeWidth={effectiveLineWidth}
                      opacity={effectiveLineOpacity}
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
                      {tick.toFixed(1)}
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
                  {isElectronicEFLabel(yAxisLabel) ? (
                    <>
                      <tspan>E − E</tspan>
                      <tspan
                        baselineShift="sub"
                        fontSize={Math.max(8, axisLabelFontSize * 0.72)}
                      >
                        F
                      </tspan>
                      <tspan baselineShift="baseline"> (eV)</tspan>
                    </>
                  ) : (
                    yAxisLabel
                  )}
                </text>
              </g>

              {/* X-axis line */}
              <line x1={0} y1={plotHeight} x2={plotWidth} y2={plotHeight} stroke={colors.axis} />
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
                  {valueLabel} = {hoveredPoint.energy.toFixed(valueDecimals)} {valueUnit}
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
            <span className={displayedBandGap ? "band-gap-info" : "metal-info"}>
              {displayedBandGap
                ? `${displayedBandGap.value.toFixed(3)} eV ${
                  displayedBandGap.is_direct ? "direct" : "indirect"
                } gap`
                : "No gap detected at current E_F"}
            </span>
          )}
          {showProjectionSummary && (
            <span className="band-plot-projection-pill">
              Fat bands: {projectionLabel}
            </span>
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
                        Referenced to the current {activeFermiSourceLabel.toLowerCase()} Fermi
                        level.
                      </div>
                    </>
                  ) : (
                    <div className="band-control-warning">
                      No band gap was detected for the current Fermi reference.
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
                      Enable Fat Bands
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
                        <label>Projection Group Type</label>
                        <select
                          value={projectionMode}
                          onChange={(event) =>
                            setProjectionMode(event.target.value as ProjectionMode)
                          }
                        >
                          <option value="atom">By element</option>
                          <option value="orbital">By orbital</option>
                        </select>
                      </div>

                      {projectionMode === "orbital" && hasElementResolvedOrbitals && (
                        <div className="band-control-row">
                          <label>Element for Orbital</label>
                          <select
                            value={selectedOrbitalElementKey}
                            onChange={(event) => setSelectedOrbitalElementKey(event.target.value)}
                          >
                            {orbitalElementOptions.map((option) => (
                              <option key={option.key} value={option.key}>
                                {option.label}
                              </option>
                            ))}
                          </select>
                        </div>
                      )}

                      {projectionMode === "orbital" && !hasElementResolvedOrbitals && (
                        <div className="band-control-warning">
                          Element-resolved orbital channels are unavailable for this saved result.
                          Showing global orbital totals.
                        </div>
                      )}

                      <div className="band-control-row">
                        <label>
                          {projectionMode === "orbital"
                            ? "Orbital Channel"
                            : "Contribution"}
                        </label>
                        <select
                          value={selectedProjectionId}
                          onChange={(event) => setSelectedProjectionId(event.target.value)}
                        >
                          {projectionGroups.map((group) => (
                            <option key={group.id} value={group.id}>
                              {group.label}
                            </option>
                          ))}
                        </select>
                      </div>

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

                      {fatColorMode === "accent" && (
                        <div className="band-control-row">
                          <label>Accent Color</label>
                          <input
                            type="color"
                            value={fatAccentColor}
                            onChange={(event) => setFatAccentColor(event.target.value)}
                          />
                        </div>
                      )}

                      <div className="band-control-row">
                        <label>Keep Band Lines Visible</label>
                        <input
                          type="checkbox"
                          checked={showLinesWithFat}
                          onChange={(event) => setShowLinesWithFat(event.target.checked)}
                        />
                      </div>
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
