import { useState, useRef, useMemo, useCallback, useEffect, useId } from "react";

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

type ColorMode = "single" | "rainbow";
type RainbowPalette = "jet" | "sinebow";
type ProjectionMode = "atom" | "orbital";
type ProjectionNormalizeMode = "global" | "band";
type FatColorMode = "accent" | "band";

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

function Tooltip({ text }: { text: string }) {
  return (
    <span className="tooltip-container">
      <span className="tooltip-icon">?</span>
      <span className="tooltip-text">{text}</span>
    </span>
  );
}

export function BandPlot({
  data,
  width = 700,
  height = 500,
  energyRange,
  showFermiLevel = true,
  scfFermiEnergy,
  yAxisLabel = "E − E_F (eV)",
  pointLabel = "Band",
  valueLabel = "E − E_F",
  valueUnit = "eV",
  valueDecimals = 3,
  primaryCountLabel = "bands",
  secondaryCountLabel = "k-points",
  scrollHint = "Scroll: zoom Y | Shift+Scroll: pan energy",
  yClampRange = [-25, 25],
  viewerType = "electronic",
}: BandPlotProps) {
  const svgRef = useRef<SVGSVGElement>(null);
  const clipPathId = useId();
  const [hoveredPoint, setHoveredPoint] = useState<HoveredPoint | null>(null);
  const [isHoveringPlot, setIsHoveringPlot] = useState(false);

  // Y-axis energy window (adjustable via scroll)
  const [yMin, setYMin] = useState<number | null>(null);
  const [yMax, setYMax] = useState<number | null>(null);

  // Appearance controls
  const [lineWidth, setLineWidth] = useState(1.5);
  const [lineOpacity, setLineOpacity] = useState(0.85);
  const [colorMode, setColorMode] = useState<ColorMode>("single");
  const [singleBandColor, setSingleBandColor] = useState("#1565c0");
  const [rainbowPalette, setRainbowPalette] = useState<RainbowPalette>("jet");

  // Fat-band controls
  const [fatBandsEnabled, setFatBandsEnabled] = useState(false);
  const [projectionMode, setProjectionMode] = useState<ProjectionMode>("atom");
  const [selectedProjectionId, setSelectedProjectionId] = useState("");
  const [projectionNormalizeMode, setProjectionNormalizeMode] =
    useState<ProjectionNormalizeMode>("global");
  const [fatScale, setFatScale] = useState(8);
  const [fatOpacity, setFatOpacity] = useState(0.45);
  const [fatColorMode, setFatColorMode] = useState<FatColorMode>("band");
  const [fatAccentColor, setFatAccentColor] = useState("#f57c00");
  const [showLinesWithFat, setShowLinesWithFat] = useState(true);

  // UI section toggles
  const [appearanceExpanded, setAppearanceExpanded] = useState(true);
  const [projectionExpanded, setProjectionExpanded] = useState(false);
  const [exportNote, setExportNote] = useState("");

  // Margins
  const margin = { top: 30, right: 30, bottom: 50, left: 70 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;

  // Use SCF Fermi energy if available, otherwise fall back to data.fermi_energy
  const fermiEnergy = scfFermiEnergy ?? data.fermi_energy;

  // Shift all energies relative to Fermi level (E - E_F)
  const shiftedEnergies = useMemo(() => {
    return data.energies.map((band) => band.map((e) => e - fermiEnergy));
  }, [data.energies, fermiEnergy]);

  // Calculate shifted energy range
  const shiftedEnergyRange: [number, number] = useMemo(() => {
    return [
      data.energy_range[0] - fermiEnergy,
      data.energy_range[1] - fermiEnergy,
    ];
  }, [data.energy_range, fermiEnergy]);

  // Calculate scales - X is always fixed, Y is adjustable
  const scales = useMemo(() => {
    const kMin = data.k_points.length > 0 ? data.k_points[0] : 0;
    const kMax = data.k_points[data.k_points.length - 1] ?? 1;
    const kSpan = Math.abs(kMax - kMin) > 1e-12 ? kMax - kMin : 1;

    // Use provided energy range, custom Y range, or calculate from shifted data with padding
    let eMin: number;
    let eMax: number;

    if (yMin !== null && yMax !== null) {
      eMin = yMin;
      eMax = yMax;
    } else if (energyRange) {
      [eMin, eMax] = energyRange;
    } else {
      [eMin, eMax] = shiftedEnergyRange;
      const span = Math.max(eMax - eMin, 1e-6);
      const padding = span * 0.1;
      eMin -= padding;
      eMax += padding;

      // Clamp to reasonable range around Fermi level (now at 0) if too wide
      const maxRange = 20; // eV
      if (eMax - eMin > maxRange * 2) {
        eMin = -maxRange;
        eMax = maxRange;
      }
    }

    const eSpan = Math.abs(eMax - eMin) > 1e-12 ? eMax - eMin : 1;

    return {
      kMin,
      kMax,
      eMin,
      eMax,
      xScale: (k: number) => ((k - kMin) / kSpan) * plotWidth,
      yScale: (e: number) => plotHeight - ((e - eMin) / eSpan) * plotHeight,
    };
  }, [data.k_points, energyRange, plotWidth, plotHeight, yMin, yMax, shiftedEnergyRange]);

  const bandColors = useMemo(
    () =>
      shiftedEnergies.map((_, idx) =>
        bandColorForIndex(
          idx,
          shiftedEnergies.length,
          colorMode,
          singleBandColor,
          rainbowPalette,
        ),
      ),
    [shiftedEnergies, colorMode, singleBandColor, rainbowPalette],
  );

  const hasProjectionData = useMemo(() => {
    const atomCount = data.projections?.atom_groups?.length ?? 0;
    const orbitalCount = data.projections?.orbital_groups?.length ?? 0;
    return atomCount > 0 || orbitalCount > 0;
  }, [data.projections]);

  const projectionGroups = useMemo(() => {
    if (!hasProjectionData) return [];
    if (projectionMode === "atom") {
      const atomGroups = data.projections?.atom_groups ?? [];
      return aggregateElementProjectionGroups(atomGroups, data.n_bands, data.n_kpoints);
    }
    return data.projections?.orbital_groups ?? [];
  }, [data.projections, data.n_bands, data.n_kpoints, hasProjectionData, projectionMode]);

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

  const selectedProjectionGroup = useMemo(() => {
    if (!selectedProjectionId) return null;
    return projectionGroups.find((group) => group.id === selectedProjectionId) || null;
  }, [projectionGroups, selectedProjectionId]);

  const normalizedProjectionWeights = useMemo(() => {
    const weights = selectedProjectionGroup?.weights;
    if (!weights || weights.length === 0) return null;

    if (projectionNormalizeMode === "band") {
      return weights.map((bandWeights) => {
        const maxBandWeight = bandWeights.reduce((max, value) => {
          if (!Number.isFinite(value)) return max;
          return Math.max(max, value);
        }, 0);
        const denom = maxBandWeight > 0 ? maxBandWeight : 1;
        return bandWeights.map((value) => clamp01((Number.isFinite(value) ? value : 0) / denom));
      });
    }

    let globalMax = 0;
    for (const bandWeights of weights) {
      for (const value of bandWeights) {
        if (Number.isFinite(value) && value > globalMax) {
          globalMax = value;
        }
      }
    }
    const denom = globalMax > 0 ? globalMax : 1;
    return weights.map((bandWeights) =>
      bandWeights.map((value) => clamp01((Number.isFinite(value) ? value : 0) / denom)),
    );
  }, [selectedProjectionGroup, projectionNormalizeMode]);

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
        path += ` L ${scales.xScale(kPoints[i])},${scales.yScale(band[i])}`;
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

      // Clamp to limits if configured
      if (yClampRange) {
        const [minLimit, maxLimit] = yClampRange;
        if (newMin < minLimit) {
          newMin = minLimit;
        }
        if (newMax > maxLimit) {
          newMax = maxLimit;
        }
      }

      setYMin(newMin);
      setYMax(newMax);
    },
    [yMin, yMax, scales, yClampRange],
  );

  // Reset view
  const resetView = useCallback(() => {
    setYMin(null);
    setYMax(null);
  }, []);

  const handleExportPlaceholder = useCallback(() => {
    setExportNote("Export UI stub added. Wire export logic when ready.");
  }, []);

  // Y-axis ticks
  const yTicks = useMemo(() => {
    const range = scales.eMax - scales.eMin;
    const step =
      range > 500
        ? 100
        : range > 200
          ? 50
          : range > 100
            ? 25
            : range > 20
              ? 5
              : range > 10
                ? 2
                : range > 5
                  ? 1
                  : 0.5;
    const ticks: number[] = [];
    let tick = Math.ceil(scales.eMin / step) * step;
    while (tick <= scales.eMax) {
      ticks.push(tick);
      tick += step;
    }
    return ticks;
  }, [scales]);

  // Prevent page scroll while interacting with the plot area.
  useEffect(() => {
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
  }, [isHoveringPlot]);

  const drawBandLines = !fatBandsActive || showLinesWithFat;
  const projectionLabel = selectedProjectionGroup?.label || "None";
  const showProjectionSummary = fatBandsActive && selectedProjectionGroup !== null;

  return (
    <div className="band-plot-container">
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
              <div className="band-control-row">
                <label>
                  Line Thickness
                  <span className="band-control-tech-name">strokeWidth</span>
                  <Tooltip text="Controls the line width for each band trajectory. Thicker lines improve readability for dense plots but can obscure close crossings." />
                </label>
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
                <label>
                  Line Opacity
                  <span className="band-control-tech-name">opacity</span>
                  <Tooltip text="Lower opacity helps visualize many bands simultaneously and keeps fat-band overlays visible." />
                </label>
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
                <label>
                  Band Color Mode
                  <span className="band-control-tech-name">plot style</span>
                  <Tooltip text="Single color gives the traditional publication look. Rainbow assigns each band a distinct color, similar to your group notebook style." />
                </label>
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
                  <label>
                    Band Color
                    <Tooltip text="Primary color for all bands when single-color mode is active." />
                  </label>
                  <input
                    type="color"
                    value={singleBandColor}
                    onChange={(event) => setSingleBandColor(event.target.value)}
                  />
                </div>
              )}

              {colorMode === "rainbow" && (
                <div className="band-control-row">
                  <label>
                    Rainbow Palette
                    <Tooltip text="`jet` mirrors classic scientific plotting palettes; `sinebow` gives smoother hue transitions." />
                  </label>
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
            </div>
          )}
        </div>

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
                    <Tooltip text="Uses QE orbital projection weights to draw point thickness/intensity. Requires projections to be computed during the band run." />
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
                      <label>
                        Projection Group Type
                        <Tooltip text="Choose whether projection weights are grouped by chemical element totals (all sites summed) or by orbital angular momentum (s/p/d/f)." />
                      </label>
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

                    <div className="band-control-row">
                      <label>
                        Selected Contribution
                        <Tooltip text="The selected element-total or orbital channel whose projection weights are used for fat-band sizing." />
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
                      <label>
                        Weight Normalization
                        <span className="band-control-tech-name">weight scaling</span>
                        <Tooltip text="Global uses one maximum across all bands; per-band normalization boosts weak bands and emphasizes shape over absolute intensity." />
                      </label>
                      <select
                        value={projectionNormalizeMode}
                        onChange={(event) =>
                          setProjectionNormalizeMode(
                            event.target.value as ProjectionNormalizeMode,
                          )
                        }
                      >
                        <option value="global">Global</option>
                        <option value="band">Per band</option>
                      </select>
                    </div>

                    <div className="band-control-row">
                      <label>
                        Fat-Band Scale
                        <Tooltip text="Multiplies marker radius from projection weight. Larger values increase visual contrast but can overfill dense regions." />
                      </label>
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
                      <label>
                        Fat-Band Opacity
                        <Tooltip text="Opacity for projection markers. Lower opacity helps retain underlying band lines in crowded regions." />
                      </label>
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
                      <label>
                        Fat-Band Color Mode
                        <Tooltip text="Match-band uses each band's line color (including rainbow). Accent uses one projection color for all markers." />
                      </label>
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
                        <label>
                          Fat-Band Accent Color
                          <Tooltip text="Used when color mode is set to single accent." />
                        </label>
                        <input
                          type="color"
                          value={fatAccentColor}
                          onChange={(event) => setFatAccentColor(event.target.value)}
                        />
                      </div>
                    )}

                    <div className="band-control-row">
                      <label>
                        Keep Band Lines Visible
                        <Tooltip text="When enabled, the regular band lines stay visible underneath fat-band markers." />
                      </label>
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

      <svg
        ref={svgRef}
        width={width}
        height={height}
        onMouseEnter={() => setIsHoveringPlot(true)}
        onMouseMove={handleMouseMove}
        onMouseLeave={() => {
          setIsHoveringPlot(false);
          setHoveredPoint(null);
        }}
        onWheel={handleWheel}
        style={{ cursor: "crosshair" }}
      >
        <defs>
          <clipPath id={clipPathId}>
            <rect x={0} y={0} width={plotWidth} height={plotHeight} />
          </clipPath>
        </defs>

        {/* White background */}
        <rect width={width} height={height} fill="#ffffff" />

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
                stroke="#999"
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
                stroke="#333"
                strokeWidth={0.5}
              />
              <text
                x={scales.xScale(point.k_distance)}
                y={plotHeight + 20}
                textAnchor="middle"
                fill="#000"
                fontSize={14}
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
                x={plotWidth + 5}
                y={scales.yScale(0) + 4}
                fill="#d32f2f"
                fontSize={11}
              >
                E_F
              </text>
            </g>
          )}

          {/* Band lines + fat points */}
          <g clipPath={`url(#${clipPathId})`}>
            {drawBandLines &&
              shiftedEnergies.map((band, bandIdx) => (
                <path
                  key={bandIdx}
                  d={bandToPath(band, data.k_points)}
                  fill="none"
                  stroke={bandColors[bandIdx]}
                  strokeWidth={lineWidth}
                  opacity={lineOpacity}
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
            <line x1={0} y1={0} x2={0} y2={plotHeight} stroke="#333" />
            {yTicks.map((tick) => (
              <g key={tick}>
                <line
                  x1={-5}
                  x2={0}
                  y1={scales.yScale(tick)}
                  y2={scales.yScale(tick)}
                  stroke="#333"
                />
                <text
                  x={-10}
                  y={scales.yScale(tick) + 4}
                  textAnchor="end"
                  fill="#333"
                  fontSize={11}
                >
                  {tick.toFixed(1)}
                </text>
              </g>
            ))}
            <text
              transform={`translate(-50, ${plotHeight / 2}) rotate(-90)`}
              textAnchor="middle"
              fill="#333"
              fontSize={14}
            >
              {yAxisLabel}
            </text>
          </g>

          {/* X-axis line */}
          <line x1={0} y1={plotHeight} x2={plotWidth} y2={plotHeight} stroke="#333" />
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
              fill="#fff"
              stroke="#ccc"
              rx={4}
              filter="drop-shadow(0 2px 4px rgba(0,0,0,0.1))"
            />
            <text x={8} y={-24} fill="#333" fontSize={11}>
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
            E_F = {scfFermiEnergy != null ? `${scfFermiEnergy.toFixed(3)} eV` : "N/A"}
          </span>
        )}
        {showProjectionSummary && (
          <span className="band-plot-projection-pill">
            Fat bands: {projectionLabel}
          </span>
        )}
      </div>
    </div>
  );
}
