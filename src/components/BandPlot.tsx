import { useState, useRef, useMemo, useCallback } from "react";

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

export interface BandData {
  k_points: number[];
  energies: number[][];
  fermi_energy: number;
  high_symmetry_points: HighSymmetryMarker[];
  n_bands: number;
  n_kpoints: number;
  band_gap: BandGap | null;
  energy_range: [number, number];
}

interface BandPlotProps {
  data: BandData;
  width?: number;
  height?: number;
  energyRange?: [number, number];
  showFermiLevel?: boolean;
  /** Actual Fermi energy from SCF calculation (bands always reports 0) */
  scfFermiEnergy?: number;
}

interface HoveredPoint {
  band: number;
  k: number;
  energy: number;
  x: number;
  y: number;
}

// Format high-symmetry point labels (handle Greek letters)
function formatLabel(label: string): string {
  const greekMap: Record<string, string> = {
    "G": "Γ",
    "Gamma": "Γ",
    "GAMMA": "Γ",
    "Σ": "Σ",
    "Sigma": "Σ",
    "Delta": "Δ",
    "Lambda": "Λ",
  };
  return greekMap[label] || label;
}

export function BandPlot({
  data,
  width = 700,
  height = 500,
  energyRange,
  showFermiLevel = true,
  scfFermiEnergy,
}: BandPlotProps) {
  const svgRef = useRef<SVGSVGElement>(null);
  const [hoveredPoint, setHoveredPoint] = useState<HoveredPoint | null>(null);

  // Y-axis energy window (adjustable via scroll)
  const [yMin, setYMin] = useState<number | null>(null);
  const [yMax, setYMax] = useState<number | null>(null);

  // Margins
  const margin = { top: 30, right: 30, bottom: 50, left: 70 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;

  // Calculate scales - X is always fixed, Y is adjustable
  const scales = useMemo(() => {
    const kMin = 0;
    const kMax = data.k_points[data.k_points.length - 1] || 1;

    // Use provided energy range, custom Y range, or calculate from data with padding
    let eMin: number, eMax: number;

    if (yMin !== null && yMax !== null) {
      eMin = yMin;
      eMax = yMax;
    } else if (energyRange) {
      [eMin, eMax] = energyRange;
    } else {
      [eMin, eMax] = data.energy_range;
      const padding = (eMax - eMin) * 0.1;
      eMin -= padding;
      eMax += padding;

      // Clamp to reasonable range around Fermi level if too wide
      const maxRange = 20; // eV
      if (eMax - eMin > maxRange * 2) {
        const center = data.fermi_energy;
        eMin = center - maxRange;
        eMax = center + maxRange;
      }
    }

    return {
      kMin,
      kMax,
      eMin,
      eMax,
      xScale: (k: number) => ((k - kMin) / (kMax - kMin)) * plotWidth,
      yScale: (e: number) => plotHeight - ((e - eMin) / (eMax - eMin)) * plotHeight,
    };
  }, [data, energyRange, plotWidth, plotHeight, yMin, yMax]);

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
    [scales]
  );

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
      const k = (x / plotWidth) * scales.kMax;
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

      // Find closest band at that k-point
      let closestBandIdx = 0;
      minDist = Infinity;
      for (let b = 0; b < data.energies.length; b++) {
        const dist = Math.abs(data.energies[b][closestKIdx] - energy);
        if (dist < minDist) {
          minDist = dist;
          closestBandIdx = b;
        }
      }

      // Only show tooltip if close enough
      const bandEnergy = data.energies[closestBandIdx][closestKIdx];
      const pixelDist = Math.abs(scales.yScale(bandEnergy) - y);
      if (pixelDist < 15) {
        setHoveredPoint({
          band: closestBandIdx + 1,
          k: data.k_points[closestKIdx],
          energy: bandEnergy,
          x: scales.xScale(data.k_points[closestKIdx]),
          y: scales.yScale(bandEnergy),
        });
      } else {
        setHoveredPoint(null);
      }
    },
    [data, scales, plotWidth, plotHeight, margin]
  );

  // Handle scroll to adjust Y-axis range
  // Limit the vertical range to -25 to 25 eV
  const Y_LIMIT_MIN = -25;
  const Y_LIMIT_MAX = 25;

  const handleWheel = useCallback((e: React.WheelEvent) => {
    e.preventDefault();

    const currentMin = yMin ?? scales.eMin;
    const currentMax = yMax ?? scales.eMax;
    const range = currentMax - currentMin;

    let newMin: number, newMax: number;

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

    // Clamp to limits
    if (newMin < Y_LIMIT_MIN) {
      newMin = Y_LIMIT_MIN;
    }
    if (newMax > Y_LIMIT_MAX) {
      newMax = Y_LIMIT_MAX;
    }

    setYMin(newMin);
    setYMax(newMax);
  }, [yMin, yMax, scales]);

  // Reset view
  const resetView = useCallback(() => {
    setYMin(null);
    setYMax(null);
  }, []);

  // Y-axis ticks
  const yTicks = useMemo(() => {
    const range = scales.eMax - scales.eMin;
    const step = range > 10 ? 2 : range > 5 ? 1 : 0.5;
    const ticks: number[] = [];
    let tick = Math.ceil(scales.eMin / step) * step;
    while (tick <= scales.eMax) {
      ticks.push(tick);
      tick += step;
    }
    return ticks;
  }, [scales]);

  return (
    <div className="band-plot-container">
      <div className="band-plot-controls">
        <button onClick={resetView} className="band-plot-reset">
          Reset View
        </button>
        <span className="band-plot-hint">Scroll: zoom Y | Shift+Scroll: pan energy</span>
      </div>

      <svg
        ref={svgRef}
        width={width}
        height={height}
        onMouseMove={handleMouseMove}
        onMouseLeave={() => setHoveredPoint(null)}
        onWheel={handleWheel}
        style={{ cursor: "crosshair" }}
      >
        <defs>
          <clipPath id="plot-area">
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

          {/* Fermi level at E=0 (band calculations reference to Fermi) */}
          {showFermiLevel && (
            <g>
              <line
                x1={0}
                x2={plotWidth}
                y1={scales.yScale(data.fermi_energy)}
                y2={scales.yScale(data.fermi_energy)}
                stroke="#d32f2f"
                strokeWidth={1}
                strokeDasharray="4,4"
              />
              <text
                x={plotWidth + 5}
                y={scales.yScale(data.fermi_energy) + 4}
                fill="#d32f2f"
                fontSize={11}
              >
                E_F
              </text>
            </g>
          )}

          {/* Band lines */}
          <g clipPath="url(#plot-area)">
            {data.energies.map((band, bandIdx) => (
              <path
                key={bandIdx}
                d={bandToPath(band, data.k_points)}
                fill="none"
                stroke="#1565c0"
                strokeWidth={1.5}
                opacity={0.85}
              />
            ))}
          </g>

          {/* Band gap annotation - disabled for now
          {data.band_gap && (
            <g>
              <line
                x1={scales.xScale(data.band_gap.vbm_k)}
                x2={scales.xScale(data.band_gap.cbm_k)}
                y1={scales.yScale(data.band_gap.vbm_energy)}
                y2={scales.yScale(data.band_gap.cbm_energy)}
                stroke="#2e7d32"
                strokeWidth={2}
                strokeDasharray="6,3"
              />
              <circle
                cx={scales.xScale(data.band_gap.vbm_k)}
                cy={scales.yScale(data.band_gap.vbm_energy)}
                r={4}
                fill="#2e7d32"
              />
              <circle
                cx={scales.xScale(data.band_gap.cbm_k)}
                cy={scales.yScale(data.band_gap.cbm_energy)}
                r={4}
                fill="#2e7d32"
              />
              <text
                x={(scales.xScale(data.band_gap.vbm_k) + scales.xScale(data.band_gap.cbm_k)) / 2}
                y={(scales.yScale(data.band_gap.vbm_energy) + scales.yScale(data.band_gap.cbm_energy)) / 2 - 10}
                textAnchor="middle"
                fill="#2e7d32"
                fontSize={12}
                fontWeight="bold"
              >
                {data.band_gap.value.toFixed(2)} eV
              </text>
            </g>
          )}
          */}

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
              E − E_F (eV)
            </text>
          </g>

          {/* X-axis line */}
          <line x1={0} y1={plotHeight} x2={plotWidth} y2={plotHeight} stroke="#333" />
        </g>

        {/* Tooltip */}
        {hoveredPoint && (
          <g transform={`translate(${margin.left + hoveredPoint.x + 10}, ${margin.top + hoveredPoint.y - 10})`}>
            <rect
              x={0}
              y={-30}
              width={130}
              height={40}
              fill="#fff"
              stroke="#ccc"
              rx={4}
              filter="drop-shadow(0 2px 4px rgba(0,0,0,0.1))"
            />
            <text x={8} y={-12} fill="#333" fontSize={11}>
              Band {hoveredPoint.band}
            </text>
            <text x={8} y={4} fill="#1565c0" fontSize={11}>
              E − E_F = {hoveredPoint.energy.toFixed(3)} eV
            </text>
          </g>
        )}
      </svg>

      {/* Info panel */}
      <div className="band-plot-info">
        <span>{data.n_bands} bands</span>
        <span>{data.n_kpoints} k-points</span>
        <span>E_F = {scfFermiEnergy != null ? `${scfFermiEnergy.toFixed(3)} eV` : "N/A"}</span>
        {/* Band gap info - disabled for now
        {data.band_gap ? (
          <span className="band-gap-info">
            Gap: {data.band_gap.value.toFixed(3)} eV ({data.band_gap.is_direct ? "direct" : "indirect"})
          </span>
        ) : (
          <span className="metal-info">Metallic</span>
        )}
        */}
      </div>
    </div>
  );
}
