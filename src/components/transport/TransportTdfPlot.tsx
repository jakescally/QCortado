import { useCallback, useLayoutEffect, useMemo, useRef, useState } from "react";
import { useTheme } from "../../lib/ThemeContext";
import {
  TransportAxisSettings,
  TransportMetricDefinition,
  TransportTdfData,
  formatTransportNumber,
  getTdfComponentValue,
} from "../../lib/transport";

interface TransportTdfPlotProps {
  data: TransportTdfData;
  component: string;
  axisSettings: TransportAxisSettings;
  metricDefinition: TransportMetricDefinition;
}

interface HoveredPoint {
  energy: number;
  value: number;
  x: number;
  y: number;
}

function clamp(value: number, min: number, max: number): number {
  return Math.min(Math.max(value, min), max);
}

export function TransportTdfPlot({
  data,
  component,
  axisSettings,
  metricDefinition,
}: TransportTdfPlotProps) {
  const { isDark } = useTheme();
  const canvasRef = useRef<HTMLDivElement>(null);
  const svgRef = useRef<SVGSVGElement>(null);
  const [containerSize, setContainerSize] = useState<{ width: number; height: number } | null>(null);
  const [hoveredPoint, setHoveredPoint] = useState<HoveredPoint | null>(null);

  useLayoutEffect(() => {
    const node = canvasRef.current;
    if (!node) {
      return;
    }
    const observer = new ResizeObserver((entries) => {
      const entry = entries[0];
      if (!entry) {
        return;
      }
      const width = Math.max(260, Math.floor(entry.contentRect.width));
      const height = Math.max(240, Math.floor(entry.contentRect.height));
      setContainerSize((previous) => (
        previous && previous.width === width && previous.height === height
          ? previous
          : { width, height }
      ));
    });
    observer.observe(node);
    return () => observer.disconnect();
  }, []);

  const colors = useMemo(() => (
    isDark
      ? {
          bg: "#1e1e2e",
          axis: "#718096",
          grid: "#4a5568",
          text: "#e2e8f0",
          tooltip: "#2d3748",
          tooltipBorder: "#4a5568",
          tooltipText: "#e2e8f0",
        }
      : {
          bg: "#ffffff",
          axis: "#333333",
          grid: "#999999",
          text: "#000000",
          tooltip: "#ffffff",
          tooltipBorder: "#cccccc",
          tooltipText: "#333333",
        }
  ), [isDark]);

  const points = useMemo(() => data.energy_values_ev.map((energy, index) => ({
    energy,
    value: getTdfComponentValue(data, component, index),
  })).filter((point): point is { energy: number; value: number } => point.value != null), [component, data]);

  const width = containerSize?.width ?? 320;
  const height = containerSize?.height ?? 280;
  const margin = { top: 26, right: 22, bottom: 58, left: 74 };
  const plotWidth = Math.max(1, width - margin.left - margin.right);
  const plotHeight = Math.max(1, height - margin.top - margin.bottom);
  const xValues = points.map((point) => point.energy);
  const yValues = points.map((point) => point.value);
  const xMin = xValues.length > 0 ? Math.min(...xValues) : -1;
  const xMax = xValues.length > 0 ? Math.max(...xValues) : 1;
  const yMin = yValues.length > 0 ? Math.min(...yValues) : 0;
  const yMax = yValues.length > 0 ? Math.max(...yValues) : 1;
  const xSpan = Math.abs(xMax - xMin) < 1e-12 ? 1 : xMax - xMin;
  const ySpanRaw = Math.abs(yMax - yMin) < 1e-12 ? Math.max(Math.abs(yMax), 1) : yMax - yMin;
  const yPad = ySpanRaw * 0.12;
  const yDomainMin = yMin - yPad;
  const yDomainMax = yMax + yPad;
  const ySpan = Math.abs(yDomainMax - yDomainMin) < 1e-12 ? 1 : yDomainMax - yDomainMin;

  const toX = useCallback(
    (value: number) => ((value - xMin) / xSpan) * plotWidth,
    [plotWidth, xMin, xSpan],
  );
  const toY = useCallback(
    (value: number) => plotHeight - ((value - yDomainMin) / ySpan) * plotHeight,
    [plotHeight, yDomainMin, ySpan],
  );

  const pathData = useMemo(() => points.map((point, index) => (
    `${index === 0 ? "M" : "L"} ${toX(point.energy).toFixed(2)} ${toY(point.value).toFixed(2)}`
  )).join(" "), [points, toX, toY]);

  const xTicks = useMemo(() => {
    const steps = 4;
    return Array.from({ length: steps + 1 }, (_, index) => xMin + (index / steps) * (xMax - xMin));
  }, [xMax, xMin]);
  const yTicks = useMemo(() => {
    const steps = 4;
    return Array.from({ length: steps + 1 }, (_, index) => yDomainMin + (index / steps) * (yDomainMax - yDomainMin));
  }, [yDomainMax, yDomainMin]);

  const handleMouseMove = useCallback((event: React.MouseEvent<SVGSVGElement>) => {
    if (!svgRef.current || points.length === 0) {
      return;
    }
    const rect = svgRef.current.getBoundingClientRect();
    const mouseX = event.clientX - rect.left - margin.left;
    const mouseY = event.clientY - rect.top - margin.top;
    let best: HoveredPoint | null = null;
    let bestDistance = Number.POSITIVE_INFINITY;
    points.forEach((point) => {
      const x = toX(point.energy);
      const y = toY(point.value);
      const distance = Math.hypot(x - mouseX, y - mouseY);
      if (distance < bestDistance) {
        bestDistance = distance;
        best = { energy: point.energy, value: point.value, x, y };
      }
    });
    setHoveredPoint(bestDistance <= 20 ? best : null);
  }, [margin.left, margin.top, points, toX, toY]);

  const tooltipX = hoveredPoint ? clamp(margin.left + hoveredPoint.x + 10, margin.left + 4, width - 180) : 0;
  const tooltipY = hoveredPoint ? clamp(margin.top + hoveredPoint.y - 68, margin.top + 2, height - 76) : 0;

  return (
    <div className="transport-tdf-plot-shell">
      <div className="transport-tdf-plot-canvas" ref={canvasRef}>
        <svg
          ref={svgRef}
          width={width}
          height={height}
          viewBox={`0 0 ${width} ${height}`}
          onMouseMove={handleMouseMove}
          onMouseLeave={() => setHoveredPoint(null)}
          style={{ display: "block" }}
        >
          <rect x={0} y={0} width={width} height={height} fill={colors.bg} />

          <g transform={`translate(${margin.left}, ${margin.top})`}>
            <g opacity={0.3}>
              {xTicks.map((tick) => (
                <line
                  key={`x-grid-${tick}`}
                  x1={toX(tick)}
                  x2={toX(tick)}
                  y1={0}
                  y2={plotHeight}
                  stroke={colors.grid}
                  strokeDasharray="2,2"
                />
              ))}
              {yTicks.map((tick) => (
                <line
                  key={`y-grid-${tick}`}
                  x1={0}
                  x2={plotWidth}
                  y1={toY(tick)}
                  y2={toY(tick)}
                  stroke={colors.grid}
                  strokeDasharray="2,2"
                />
              ))}
            </g>

            <line x1={0} x2={plotWidth} y1={plotHeight} y2={plotHeight} stroke={colors.axis} />
            <line x1={0} x2={0} y1={0} y2={plotHeight} stroke={colors.axis} />

            {pathData && (
              <path
                d={pathData}
                fill="none"
                stroke={metricDefinition.accentColor}
                strokeWidth={2.6}
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            )}

            {points.map((point) => (
              <circle
                key={`${point.energy}-${point.value}`}
                cx={toX(point.energy)}
                cy={toY(point.value)}
                r={3.4}
                fill={metricDefinition.accentColor}
              />
            ))}

            {xTicks.map((tick) => (
              <g key={`x-tick-${tick}`}>
                <line
                  x1={toX(tick)}
                  x2={toX(tick)}
                  y1={plotHeight}
                  y2={plotHeight + 6}
                  stroke={colors.axis}
                />
                <text
                  x={toX(tick)}
                  y={plotHeight + 22}
                  textAnchor="middle"
                  fill={colors.text}
                  fontSize={11}
                >
                  {formatTransportNumber(tick, axisSettings.precisionMode, axisSettings.precision)}
                </text>
              </g>
            ))}
            {yTicks.map((tick) => (
              <g key={`y-tick-${tick}`}>
                <line
                  x1={-6}
                  x2={0}
                  y1={toY(tick)}
                  y2={toY(tick)}
                  stroke={colors.axis}
                />
                <text
                  x={-10}
                  y={toY(tick) + 4}
                  textAnchor="end"
                  fill={colors.text}
                  fontSize={11}
                >
                  {formatTransportNumber(tick, axisSettings.precisionMode, axisSettings.precision)}
                </text>
              </g>
            ))}

            <text
              x={plotWidth / 2}
              y={plotHeight + 44}
              textAnchor="middle"
              fill={colors.text}
              fontSize={13}
              fontFamily="serif"
              fontStyle="italic"
            >
              Energy (eV)
            </text>
            <text
              x={-54}
              y={plotHeight / 2}
              textAnchor="middle"
              fill={colors.text}
              fontSize={13}
              fontFamily="serif"
              fontStyle="italic"
              transform={`rotate(-90 -54 ${plotHeight / 2})`}
            >
              {metricDefinition.symbol} ({metricDefinition.defaultUnit})
            </text>

            {hoveredPoint && (
              <circle
                cx={hoveredPoint.x}
                cy={hoveredPoint.y}
                r={5.2}
                fill="#ff9800"
                stroke={colors.bg}
                strokeWidth={2}
              />
            )}
          </g>

          {hoveredPoint && (
            <g transform={`translate(${tooltipX}, ${tooltipY})`}>
              <rect
                x={0}
                y={0}
                width={168}
                height={56}
                fill={colors.tooltip}
                stroke={colors.tooltipBorder}
                rx={4}
              />
              <text x={10} y={18} fill={colors.tooltipText} fontSize={11}>
                {metricDefinition.symbol} · {component}
              </text>
              <text x={10} y={34} fill={colors.tooltipText} fontSize={11}>
                E = {formatTransportNumber(
                  hoveredPoint.energy,
                  axisSettings.precisionMode,
                  axisSettings.precision,
                )} eV
              </text>
              <text x={10} y={50} fill={colors.tooltipText} fontSize={11}>
                value = {formatTransportNumber(
                  hoveredPoint.value,
                  axisSettings.precisionMode,
                  axisSettings.precision,
                )}
              </text>
            </g>
          )}
        </svg>
      </div>
    </div>
  );
}
