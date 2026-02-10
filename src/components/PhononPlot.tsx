import { useMemo, useCallback, useRef, useState } from "react";

// Phonon data interfaces
interface PhononDOS {
  frequencies: number[];
  dos: number[];
  omega_max: number;
  omega_min: number;
}

interface PhononHighSymmetryMarker {
  q_distance: number;
  label: string;
}

interface PhononDispersion {
  q_points: number[];
  frequencies: number[][];
  high_symmetry_points: PhononHighSymmetryMarker[];
  n_modes: number;
  n_qpoints: number;
  frequency_range: [number, number];
}

// Format high-symmetry point labels (handle Greek letters)
function formatLabel(label: string): string {
  const greekMap: Record<string, string> = {
    G: "Γ",
    Gamma: "Γ",
    GAMMA: "Γ",
    Sigma: "Σ",
    Delta: "Δ",
    Lambda: "Λ",
  };
  return greekMap[label] || label;
}

// ============================================================================
// Phonon DOS Plot (vertical orientation: frequency on Y, DOS on X)
// ============================================================================

interface PhononDOSPlotProps {
  data: PhononDOS;
  width?: number;
  height?: number;
  sharedYScale?: (freq: number) => number;
  freqMin?: number;
  freqMax?: number;
}

export function PhononDOSPlot({
  data,
  width = 200,
  height = 500,
  sharedYScale,
  freqMin,
  freqMax,
}: PhononDOSPlotProps) {
  const margin = { top: 30, right: 20, bottom: 50, left: 50 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;

  // Calculate scales
  const scales = useMemo(() => {
    const dosMax = Math.max(...data.dos) * 1.1;

    // Use provided frequency range or calculate from data
    const fMin = freqMin !== undefined ? freqMin : data.omega_min - (data.omega_max - data.omega_min) * 0.05;
    const fMax = freqMax !== undefined ? freqMax : data.omega_max + (data.omega_max - data.omega_min) * 0.05;

    return {
      xScale: (d: number) => (d / dosMax) * plotWidth,
      yScale: sharedYScale || ((f: number) => plotHeight - ((f - fMin) / (fMax - fMin)) * plotHeight),
      fMin,
      fMax,
      dosMax,
    };
  }, [data, plotWidth, plotHeight, sharedYScale, freqMin, freqMax]);

  // Generate path for DOS curve
  const dosPath = useMemo(() => {
    if (data.frequencies.length === 0) return "";

    let path = `M 0,${scales.yScale(data.frequencies[0])}`;
    for (let i = 0; i < data.frequencies.length; i++) {
      const x = scales.xScale(data.dos[i]);
      const y = scales.yScale(data.frequencies[i]);
      path += ` L ${x},${y}`;
    }
    // Close the path for filling
    path += ` L 0,${scales.yScale(data.frequencies[data.frequencies.length - 1])}`;
    path += " Z";
    return path;
  }, [data, scales]);

  // Y-axis ticks
  const yTicks = useMemo(() => {
    const range = scales.fMax - scales.fMin;
    const step = range > 500 ? 100 : range > 200 ? 50 : range > 100 ? 25 : 10;
    const ticks: number[] = [];
    let tick = Math.ceil(scales.fMin / step) * step;
    while (tick <= scales.fMax) {
      ticks.push(tick);
      tick += step;
    }
    return ticks;
  }, [scales]);

  return (
    <svg width={width} height={height}>
      <rect width={width} height={height} fill="#ffffff" />

      <g transform={`translate(${margin.left}, ${margin.top})`}>
        {/* Grid lines */}
        <g opacity={0.3}>
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

        {/* Zero frequency line */}
        {scales.fMin <= 0 && scales.fMax >= 0 && (
          <line
            x1={0}
            x2={plotWidth}
            y1={scales.yScale(0)}
            y2={scales.yScale(0)}
            stroke="#d32f2f"
            strokeWidth={1}
            strokeDasharray="4,4"
          />
        )}

        {/* DOS curve */}
        <path d={dosPath} fill="#1976d2" fillOpacity={0.3} stroke="#1976d2" strokeWidth={1.5} />

        {/* Y-axis */}
        <g>
          <line x1={0} y1={0} x2={0} y2={plotHeight} stroke="#333" />
          {!sharedYScale &&
            yTicks.map((tick) => (
              <g key={tick}>
                <line x1={-5} x2={0} y1={scales.yScale(tick)} y2={scales.yScale(tick)} stroke="#333" />
                <text x={-10} y={scales.yScale(tick) + 4} textAnchor="end" fill="#333" fontSize={11}>
                  {tick}
                </text>
              </g>
            ))}
          {!sharedYScale && (
            <text
              transform={`translate(-40, ${plotHeight / 2}) rotate(-90)`}
              textAnchor="middle"
              fill="#333"
              fontSize={12}
            >
              Frequency (cm⁻¹)
            </text>
          )}
        </g>

        {/* X-axis */}
        <g>
          <line x1={0} y1={plotHeight} x2={plotWidth} y2={plotHeight} stroke="#333" />
          <text x={plotWidth / 2} y={plotHeight + 35} textAnchor="middle" fill="#333" fontSize={12}>
            DOS
          </text>
        </g>
      </g>
    </svg>
  );
}

// ============================================================================
// Phonon Dispersion Plot (same layout as BandPlot)
// ============================================================================

interface PhononDispersionPlotProps {
  data: PhononDispersion;
  width?: number;
  height?: number;
  sharedYScale?: (freq: number) => number;
  freqMin?: number;
  freqMax?: number;
}

export function PhononDispersionPlot({
  data,
  width = 700,
  height = 500,
  sharedYScale,
  freqMin,
  freqMax,
}: PhononDispersionPlotProps) {
  const svgRef = useRef<SVGSVGElement>(null);
  const [hoveredPoint, setHoveredPoint] = useState<{
    mode: number;
    q: number;
    freq: number;
    x: number;
    y: number;
  } | null>(null);

  // Y-axis range (adjustable via scroll)
  const [yMin, setYMin] = useState<number | null>(null);
  const [yMax, setYMax] = useState<number | null>(null);

  const margin = { top: 30, right: 30, bottom: 50, left: 70 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;

  // Calculate scales
  const scales = useMemo(() => {
    const qMin = 0;
    const qMax = data.q_points[data.q_points.length - 1] || 1;

    // Use provided frequency range, custom Y range, or calculate from data
    let fMin: number, fMax: number;

    if (yMin !== null && yMax !== null) {
      fMin = yMin;
      fMax = yMax;
    } else if (freqMin !== undefined && freqMax !== undefined) {
      fMin = freqMin;
      fMax = freqMax;
    } else {
      const rawMax = Number(data.frequency_range?.[1]);
      const upperBase = Number.isFinite(rawMax) ? Math.max(rawMax, 1) : 1;
      const padding = Math.max(5, upperBase * 0.08);
      fMin = 0;
      fMax = upperBase + padding;
    }

    return {
      qMin,
      qMax,
      fMin,
      fMax,
      xScale: (q: number) => ((q - qMin) / (qMax - qMin)) * plotWidth,
      yScale:
        sharedYScale ||
        ((f: number) => plotHeight - ((f - fMin) / (fMax - fMin)) * plotHeight),
    };
  }, [data, plotWidth, plotHeight, yMin, yMax, sharedYScale, freqMin, freqMax]);

  // Generate path for a mode
  const modeToPath = useCallback(
    (mode: number[], qPoints: number[]) => {
      if (mode.length === 0 || qPoints.length === 0) return "";

      let path = `M ${scales.xScale(qPoints[0])},${scales.yScale(mode[0])}`;
      for (let i = 1; i < mode.length && i < qPoints.length; i++) {
        path += ` L ${scales.xScale(qPoints[i])},${scales.yScale(mode[i])}`;
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

      if (x < 0 || x > plotWidth || y < 0 || y > plotHeight) {
        setHoveredPoint(null);
        return;
      }

      // Convert to data coordinates
      const q = (x / plotWidth) * scales.qMax;
      const freq = scales.fMin + (1 - y / plotHeight) * (scales.fMax - scales.fMin);

      // Find closest q-point index
      let closestQIdx = 0;
      let minDist = Infinity;
      for (let i = 0; i < data.q_points.length; i++) {
        const dist = Math.abs(data.q_points[i] - q);
        if (dist < minDist) {
          minDist = dist;
          closestQIdx = i;
        }
      }

      // Find closest mode at that q-point
      let closestModeIdx = 0;
      minDist = Infinity;
      for (let m = 0; m < data.frequencies.length; m++) {
        const dist = Math.abs(data.frequencies[m][closestQIdx] - freq);
        if (dist < minDist) {
          minDist = dist;
          closestModeIdx = m;
        }
      }

      // Only show tooltip if close enough
      const modeFreq = data.frequencies[closestModeIdx][closestQIdx];
      const pixelDist = Math.abs(scales.yScale(modeFreq) - y);
      if (pixelDist < 15) {
        setHoveredPoint({
          mode: closestModeIdx + 1,
          q: data.q_points[closestQIdx],
          freq: modeFreq,
          x: scales.xScale(data.q_points[closestQIdx]),
          y: scales.yScale(modeFreq),
        });
      } else {
        setHoveredPoint(null);
      }
    },
    [data, scales, plotWidth, plotHeight, margin]
  );

  // Handle wheel for zoom/pan
  const handleWheel = useCallback(
    (e: React.WheelEvent) => {
      e.preventDefault();

      const currentMin = yMin ?? scales.fMin;
      const currentMax = yMax ?? scales.fMax;
      const range = currentMax - currentMin;

      let newMin: number, newMax: number;

      if (e.shiftKey) {
        // Shift + scroll = pan
        const shift = e.deltaY > 0 ? range * 0.02 : -range * 0.02;
        newMin = currentMin + shift;
        newMax = currentMax + shift;
      } else {
        // Normal scroll = zoom
        const factor = e.deltaY > 0 ? 1.03 : 0.97;
        const center = (currentMax + currentMin) / 2;
        const newRange = range * factor;
        newMin = center - newRange / 2;
        newMax = center + newRange / 2;
      }

      setYMin(newMin);
      setYMax(newMax);
    },
    [yMin, yMax, scales]
  );

  // Reset view
  const resetView = useCallback(() => {
    setYMin(null);
    setYMax(null);
  }, []);

  // Y-axis ticks
  const yTicks = useMemo(() => {
    const range = scales.fMax - scales.fMin;
    const step = range > 500 ? 100 : range > 200 ? 50 : range > 100 ? 25 : 10;
    const ticks: number[] = [];
    let tick = Math.ceil(scales.fMin / step) * step;
    while (tick <= scales.fMax) {
      ticks.push(tick);
      tick += step;
    }
    return ticks;
  }, [scales]);

  return (
    <div className="phonon-dispersion-container">
      {!sharedYScale && (
        <div className="phonon-plot-controls">
          <button onClick={resetView} className="phonon-plot-reset">
            Reset View
          </button>
          <span className="phonon-plot-hint">Scroll: zoom | Shift+Scroll: pan</span>
        </div>
      )}

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
          <clipPath id="phonon-plot-area">
            <rect x={0} y={0} width={plotWidth} height={plotHeight} />
          </clipPath>
        </defs>

        <rect width={width} height={height} fill="#ffffff" />

        <g transform={`translate(${margin.left}, ${margin.top})`}>
          {/* Grid lines */}
          <g opacity={0.3}>
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
                x1={scales.xScale(point.q_distance)}
                x2={scales.xScale(point.q_distance)}
                y1={0}
                y2={plotHeight}
                stroke="#333"
                strokeWidth={0.5}
              />
              <text
                x={scales.xScale(point.q_distance)}
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

          {/* Zero frequency line (imaginary modes indicator) */}
          {scales.fMin <= 0 && scales.fMax >= 0 && (
            <line
              x1={0}
              x2={plotWidth}
              y1={scales.yScale(0)}
              y2={scales.yScale(0)}
              stroke="#d32f2f"
              strokeWidth={1}
              strokeDasharray="4,4"
            />
          )}

          {/* Mode lines */}
          <g clipPath="url(#phonon-plot-area)">
            {data.frequencies.map((mode, modeIdx) => (
              <path
                key={modeIdx}
                d={modeToPath(mode, data.q_points)}
                fill="none"
                stroke="#1565c0"
                strokeWidth={1.5}
                opacity={0.85}
              />
            ))}
          </g>

          {/* Hover point */}
          {hoveredPoint && (
            <circle
              cx={hoveredPoint.x}
              cy={hoveredPoint.y}
              r={5}
              fill="#ff9800"
              stroke="#fff"
              strokeWidth={2}
            />
          )}

          {/* Y-axis */}
          <g>
            <line x1={0} y1={0} x2={0} y2={plotHeight} stroke="#333" />
            {yTicks.map((tick) => (
              <g key={tick}>
                <line x1={-5} x2={0} y1={scales.yScale(tick)} y2={scales.yScale(tick)} stroke="#333" />
                <text x={-10} y={scales.yScale(tick) + 4} textAnchor="end" fill="#333" fontSize={11}>
                  {tick}
                </text>
              </g>
            ))}
            <text
              transform={`translate(-50, ${plotHeight / 2}) rotate(-90)`}
              textAnchor="middle"
              fill="#333"
              fontSize={14}
            >
              Frequency (cm⁻¹)
            </text>
          </g>

          {/* X-axis line */}
          <line x1={0} y1={plotHeight} x2={plotWidth} y2={plotHeight} stroke="#333" />
        </g>

        {/* Tooltip */}
        {hoveredPoint && (
          <g
            transform={`translate(${margin.left + hoveredPoint.x + 10}, ${margin.top + hoveredPoint.y - 10})`}
          >
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
              Mode {hoveredPoint.mode}
            </text>
            <text x={8} y={4} fill="#1565c0" fontSize={11}>
              {hoveredPoint.freq.toFixed(1)} cm⁻¹
            </text>
          </g>
        )}
      </svg>

      {/* Info panel */}
      {!sharedYScale && (
        <div className="phonon-plot-info">
          <span>{data.n_modes} modes</span>
          <span>{data.n_qpoints} q-points</span>
          <span>
            {data.frequency_range[0].toFixed(0)} - {data.frequency_range[1].toFixed(0)} cm⁻¹
          </span>
          {data.frequency_range[0] < 0 && (
            <span className="instability-warning">Imaginary modes present</span>
          )}
        </div>
      )}
    </div>
  );
}

// ============================================================================
// Combined Phonon Plot (Dispersion 70% | DOS 30% with shared Y-axis)
// ============================================================================

interface PhononPlotProps {
  dos: PhononDOS;
  dispersion: PhononDispersion;
  width?: number;
  height?: number;
}

export function PhononPlot({ dos, dispersion, width = 900, height = 500 }: PhononPlotProps) {
  // Y-axis range (adjustable via scroll)
  const [yMin, setYMin] = useState<number | null>(null);
  const [yMax, setYMax] = useState<number | null>(null);

  // Calculate common frequency range
  const freqRange = useMemo(() => {
    const allMax = Math.max(dos.omega_max, dispersion.frequency_range[1]);
    const upperBase = Math.max(allMax, 1);
    const padding = Math.max(5, upperBase * 0.08);

    return {
      min: yMin ?? 0,
      max: yMax ?? upperBase + padding,
    };
  }, [dos, dispersion, yMin, yMax]);

  const margin = { top: 30, right: 20, bottom: 50, left: 70 };
  const dispersionWidth = Math.floor(width * 0.7);
  const dosWidth = width - dispersionWidth;
  const plotHeight = height - margin.top - margin.bottom;

  // Shared Y-scale for both plots
  const sharedYScale = useCallback(
    (freq: number) => {
      return (
        plotHeight -
        ((freq - freqRange.min) / (freqRange.max - freqRange.min)) * plotHeight
      );
    },
    [freqRange, plotHeight]
  );

  // Handle wheel for zoom/pan (applied to container)
  const handleWheel = useCallback(
    (e: React.WheelEvent) => {
      e.preventDefault();

      const currentMin = freqRange.min;
      const currentMax = freqRange.max;
      const range = currentMax - currentMin;

      let newMin: number, newMax: number;

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

      setYMin(newMin);
      setYMax(newMax);
    },
    [freqRange]
  );

  const resetView = useCallback(() => {
    setYMin(null);
    setYMax(null);
  }, []);

  return (
    <div className="phonon-combined-plot" onWheel={handleWheel}>
      <div className="phonon-plot-controls">
        <button onClick={resetView} className="phonon-plot-reset">
          Reset View
        </button>
        <span className="phonon-plot-hint">Scroll: zoom | Shift+Scroll: pan</span>
      </div>

      <div className="phonon-combined-container">
        {/* Dispersion plot (left, 70%) */}
        <PhononDispersionPlot
          data={dispersion}
          width={dispersionWidth}
          height={height}
          sharedYScale={sharedYScale}
          freqMin={freqRange.min}
          freqMax={freqRange.max}
        />

        {/* DOS plot (right, 30%) */}
        <PhononDOSPlot
          data={dos}
          width={dosWidth}
          height={height}
          sharedYScale={sharedYScale}
          freqMin={freqRange.min}
          freqMax={freqRange.max}
        />
      </div>

      {/* Info panel */}
      <div className="phonon-plot-info">
        <span>{dispersion.n_modes} modes</span>
        <span>{dispersion.n_qpoints} q-points</span>
        <span>
          {dispersion.frequency_range[0].toFixed(0)} - {dispersion.frequency_range[1].toFixed(0)}{" "}
          cm⁻¹
        </span>
        {dispersion.frequency_range[0] < 0 && (
          <span className="instability-warning">Imaginary modes present</span>
        )}
      </div>
    </div>
  );
}
