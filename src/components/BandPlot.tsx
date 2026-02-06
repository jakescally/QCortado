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
}

interface HoveredPoint {
  band: number;
  k: number;
  energy: number;
  x: number;
  y: number;
}

export function BandPlot({
  data,
  width = 700,
  height = 500,
  energyRange,
  showFermiLevel = true,
}: BandPlotProps) {
  const svgRef = useRef<SVGSVGElement>(null);
  const [hoveredPoint, setHoveredPoint] = useState<HoveredPoint | null>(null);
  const [zoom, setZoom] = useState(1);
  const [pan, setPan] = useState({ x: 0, y: 0 });
  const [isDragging, setIsDragging] = useState(false);
  const [dragStart, setDragStart] = useState({ x: 0, y: 0 });

  // Margins
  const margin = { top: 30, right: 30, bottom: 50, left: 70 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;

  // Calculate scales
  const scales = useMemo(() => {
    const kMin = 0;
    const kMax = data.k_points[data.k_points.length - 1] || 1;

    // Use provided energy range or calculate from data with padding
    let [eMin, eMax] = energyRange || data.energy_range;
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

    return {
      kMin,
      kMax,
      eMin,
      eMax,
      xScale: (k: number) => ((k - kMin) / (kMax - kMin)) * plotWidth,
      yScale: (e: number) => plotHeight - ((e - eMin) / (eMax - eMin)) * plotHeight,
    };
  }, [data, energyRange, plotWidth, plotHeight]);

  // Generate SVG path for a band
  const bandToPath = useCallback(
    (band: number[], kPoints: number[]) => {
      if (band.length === 0) return "";

      const points = band.map((e, i) => {
        const x = scales.xScale(kPoints[i]);
        const y = scales.yScale(e);
        return `${x},${y}`;
      });

      return `M ${points.join(" L ")}`;
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

  // Handle zoom
  const handleWheel = useCallback((e: React.WheelEvent) => {
    e.preventDefault();
    const delta = e.deltaY > 0 ? 0.9 : 1.1;
    setZoom((z) => Math.max(0.5, Math.min(5, z * delta)));
  }, []);

  // Handle pan
  const handleMouseDown = useCallback((e: React.MouseEvent) => {
    if (e.button === 0) {
      setIsDragging(true);
      setDragStart({ x: e.clientX - pan.x, y: e.clientY - pan.y });
    }
  }, [pan]);

  const handleMouseUp = useCallback(() => {
    setIsDragging(false);
  }, []);

  const handleDrag = useCallback(
    (e: React.MouseEvent) => {
      if (isDragging) {
        setPan({
          x: e.clientX - dragStart.x,
          y: e.clientY - dragStart.y,
        });
      }
    },
    [isDragging, dragStart]
  );

  // Reset view
  const resetView = useCallback(() => {
    setZoom(1);
    setPan({ x: 0, y: 0 });
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
        <span className="band-plot-zoom-label">Zoom: {zoom.toFixed(1)}x</span>
      </div>

      <svg
        ref={svgRef}
        width={width}
        height={height}
        onMouseMove={(e) => {
          handleMouseMove(e);
          handleDrag(e);
        }}
        onMouseDown={handleMouseDown}
        onMouseUp={handleMouseUp}
        onMouseLeave={() => {
          setHoveredPoint(null);
          setIsDragging(false);
        }}
        onWheel={handleWheel}
        style={{ cursor: isDragging ? "grabbing" : "crosshair" }}
      >
        <defs>
          <clipPath id="plot-area">
            <rect x={0} y={0} width={plotWidth} height={plotHeight} />
          </clipPath>
        </defs>

        {/* Background */}
        <rect width={width} height={height} fill="#1e1e2e" />

        {/* Plot area */}
        <g transform={`translate(${margin.left}, ${margin.top})`}>
          {/* Grid lines */}
          <g className="grid-lines" opacity={0.2}>
            {yTicks.map((tick) => (
              <line
                key={tick}
                x1={0}
                x2={plotWidth}
                y1={scales.yScale(tick)}
                y2={scales.yScale(tick)}
                stroke="#888"
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
                stroke="#666"
                strokeWidth={1}
              />
              <text
                x={scales.xScale(point.k_distance)}
                y={plotHeight + 25}
                textAnchor="middle"
                fill="#cdd6f4"
                fontSize={14}
                fontFamily="serif"
              >
                {point.label}
              </text>
            </g>
          ))}

          {/* Fermi level */}
          {showFermiLevel && (
            <g>
              <line
                x1={0}
                x2={plotWidth}
                y1={scales.yScale(data.fermi_energy)}
                y2={scales.yScale(data.fermi_energy)}
                stroke="#f38ba8"
                strokeWidth={1.5}
                strokeDasharray="4,4"
              />
              <text
                x={plotWidth + 5}
                y={scales.yScale(data.fermi_energy) + 4}
                fill="#f38ba8"
                fontSize={11}
              >
                E_F
              </text>
            </g>
          )}

          {/* Band lines */}
          <g clipPath="url(#plot-area)">
            <g transform={`translate(${pan.x}, ${pan.y}) scale(${zoom})`}>
              {data.energies.map((band, bandIdx) => (
                <path
                  key={bandIdx}
                  d={bandToPath(band, data.k_points)}
                  fill="none"
                  stroke="#89b4fa"
                  strokeWidth={1.5 / zoom}
                  opacity={0.9}
                />
              ))}
            </g>
          </g>

          {/* Band gap annotation */}
          {data.band_gap && (
            <g>
              <line
                x1={scales.xScale(data.band_gap.vbm_k)}
                x2={scales.xScale(data.band_gap.cbm_k)}
                y1={scales.yScale(data.band_gap.vbm_energy)}
                y2={scales.yScale(data.band_gap.cbm_energy)}
                stroke="#a6e3a1"
                strokeWidth={2}
                strokeDasharray="6,3"
              />
              <circle
                cx={scales.xScale(data.band_gap.vbm_k)}
                cy={scales.yScale(data.band_gap.vbm_energy)}
                r={4}
                fill="#a6e3a1"
              />
              <circle
                cx={scales.xScale(data.band_gap.cbm_k)}
                cy={scales.yScale(data.band_gap.cbm_energy)}
                r={4}
                fill="#a6e3a1"
              />
              <text
                x={(scales.xScale(data.band_gap.vbm_k) + scales.xScale(data.band_gap.cbm_k)) / 2}
                y={(scales.yScale(data.band_gap.vbm_energy) + scales.yScale(data.band_gap.cbm_energy)) / 2 - 10}
                textAnchor="middle"
                fill="#a6e3a1"
                fontSize={12}
                fontWeight="bold"
              >
                {data.band_gap.value.toFixed(2)} eV
              </text>
            </g>
          )}

          {/* Hover point */}
          {hoveredPoint && (
            <g>
              <circle
                cx={hoveredPoint.x}
                cy={hoveredPoint.y}
                r={5}
                fill="#fab387"
                stroke="#1e1e2e"
                strokeWidth={2}
              />
            </g>
          )}

          {/* Y-axis */}
          <g>
            <line x1={0} y1={0} x2={0} y2={plotHeight} stroke="#888" />
            {yTicks.map((tick) => (
              <g key={tick}>
                <line
                  x1={-5}
                  x2={0}
                  y1={scales.yScale(tick)}
                  y2={scales.yScale(tick)}
                  stroke="#888"
                />
                <text
                  x={-10}
                  y={scales.yScale(tick) + 4}
                  textAnchor="end"
                  fill="#cdd6f4"
                  fontSize={11}
                >
                  {tick.toFixed(1)}
                </text>
              </g>
            ))}
            <text
              transform={`translate(-50, ${plotHeight / 2}) rotate(-90)`}
              textAnchor="middle"
              fill="#cdd6f4"
              fontSize={14}
            >
              Energy (eV)
            </text>
          </g>

          {/* X-axis */}
          <line x1={0} y1={plotHeight} x2={plotWidth} y2={plotHeight} stroke="#888" />
        </g>

        {/* Tooltip */}
        {hoveredPoint && (
          <g transform={`translate(${margin.left + hoveredPoint.x + 10}, ${margin.top + hoveredPoint.y - 10})`}>
            <rect
              x={0}
              y={-30}
              width={120}
              height={40}
              fill="#313244"
              stroke="#45475a"
              rx={4}
            />
            <text x={8} y={-12} fill="#cdd6f4" fontSize={11}>
              Band {hoveredPoint.band}
            </text>
            <text x={8} y={4} fill="#89b4fa" fontSize={11}>
              E = {hoveredPoint.energy.toFixed(3)} eV
            </text>
          </g>
        )}
      </svg>

      {/* Info panel */}
      <div className="band-plot-info">
        <span>{data.n_bands} bands</span>
        <span>{data.n_kpoints} k-points</span>
        {data.band_gap ? (
          <span className="band-gap-info">
            Gap: {data.band_gap.value.toFixed(3)} eV ({data.band_gap.is_direct ? "direct" : "indirect"})
          </span>
        ) : (
          <span className="metal-info">Metallic</span>
        )}
      </div>
    </div>
  );
}
