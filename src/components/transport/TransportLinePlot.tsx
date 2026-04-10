import { useCallback, useEffect, useId, useLayoutEffect, useMemo, useRef, useState } from "react";
import { useTheme } from "../../lib/ThemeContext";
import {
  TransportAxisMode,
  TransportAxisSettings,
  TransportMetricDefinition,
  TransportSeriesPoint,
  buildTransportTicks,
  formatTransportNumber,
} from "../../lib/transport";

interface TransportLinePlotProps {
  lockScopeKey: string;
  metricDefinition: TransportMetricDefinition;
  axisMode: TransportAxisMode;
  axisSettings: TransportAxisSettings;
  component: string;
  series: TransportSeriesPoint[];
  xAxisLabel: string;
  yAxisLabel: string;
  valueUnit: string;
  fixedContextLabel: string;
  onCustomYRangeChange: (range: [number, number] | null) => void;
}

interface HoveredPoint {
  point: TransportSeriesPoint;
  x: number;
  y: number;
}

function clamp(value: number, min: number, max: number): number {
  return Math.min(Math.max(value, min), max);
}

function getUniqueSeriesValues(values: number[]): number[] {
  const seen = new Set<string>();
  return values.filter((value) => {
    const key = value.toFixed(6);
    if (seen.has(key)) {
      return false;
    }
    seen.add(key);
    return true;
  });
}

export function TransportLinePlot({
  lockScopeKey,
  metricDefinition,
  axisMode,
  axisSettings,
  component,
  series,
  xAxisLabel,
  yAxisLabel,
  valueUnit,
  fixedContextLabel,
  onCustomYRangeChange,
}: TransportLinePlotProps) {
  const { isDark } = useTheme();
  const canvasRef = useRef<HTMLDivElement>(null);
  const svgRef = useRef<SVGSVGElement>(null);
  const clipPathId = useId();
  const [containerSize, setContainerSize] = useState<{ width: number; height: number } | null>(null);
  const [hoveredPoint, setHoveredPoint] = useState<HoveredPoint | null>(null);
  const [lockedAutoYRange, setLockedAutoYRange] = useState<[number, number] | null>(null);
  const previousLockScopeKeyRef = useRef<string | null>(null);

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
      const nextWidth = Math.max(320, Math.floor(entry.contentRect.width));
      const nextHeight = Math.max(320, Math.floor(entry.contentRect.height));
      setContainerSize((previous) => {
        if (previous && previous.width === nextWidth && previous.height === nextHeight) {
          return previous;
        }
        return { width: nextWidth, height: nextHeight };
      });
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

  const width = containerSize?.width ?? 780;
  const height = containerSize?.height ?? 500;
  const margin = { top: 28, right: 36, bottom: 68, left: 88 };
  const plotWidth = Math.max(1, width - margin.left - margin.right);
  const plotHeight = Math.max(1, height - margin.top - margin.bottom);

  const validSeries = useMemo(
    () => series.filter((point) => Number.isFinite(point.x) && point.y != null && Number.isFinite(point.y)),
    [series],
  );

  const xValues = validSeries.map((point) => point.x);
  const yValues = validSeries.map((point) => Number(point.y));
  const xMin = xValues.length > 0 ? Math.min(...xValues) : 0;
  const xMax = xValues.length > 0 ? Math.max(...xValues) : 1;
  const yMinAuto = yValues.length > 0 ? Math.min(...yValues) : 0;
  const yMaxAuto = yValues.length > 0 ? Math.max(...yValues) : 1;
  const ySpanAuto = Math.abs(yMaxAuto - yMinAuto) < 1e-12
    ? Math.max(Math.abs(yMaxAuto), 1)
    : yMaxAuto - yMinAuto;
  const autoPadding = ySpanAuto * 0.1;
  const autoYRange = useMemo<[number, number]>(() => ([
    yMinAuto - autoPadding,
    yMaxAuto + autoPadding,
  ]), [autoPadding, yMaxAuto, yMinAuto]);

  useEffect(() => {
    if (!axisSettings.keepYRangeAcrossSlices) {
      previousLockScopeKeyRef.current = lockScopeKey;
      setLockedAutoYRange((previous) => (previous == null ? previous : null));
      return;
    }

    const plotChanged = previousLockScopeKeyRef.current !== lockScopeKey;
    previousLockScopeKeyRef.current = lockScopeKey;
    setLockedAutoYRange((previous) => {
      if (!plotChanged && previous != null) {
        return previous;
      }
      if (
        previous
        && previous[0] === autoYRange[0]
        && previous[1] === autoYRange[1]
      ) {
        return previous;
      }
      return autoYRange;
    });
  }, [autoYRange, axisSettings.keepYRangeAcrossSlices, lockScopeKey]);

  const resolvedAutoYRange = axisSettings.keepYRangeAcrossSlices
    ? lockedAutoYRange ?? autoYRange
    : autoYRange;
  const yMin = axisSettings.customYRange?.[0] ?? resolvedAutoYRange[0];
  const yMax = axisSettings.customYRange?.[1] ?? resolvedAutoYRange[1];
  const xSpan = Math.abs(xMax - xMin) < 1e-12 ? 1 : xMax - xMin;
  const ySpan = Math.abs(yMax - yMin) < 1e-12 ? 1 : yMax - yMin;

  const toX = useCallback(
    (value: number) => ((value - xMin) / xSpan) * plotWidth,
    [plotWidth, xMin, xSpan],
  );
  const toY = useCallback(
    (value: number) => plotHeight - ((value - yMin) / ySpan) * plotHeight,
    [plotHeight, yMin, ySpan],
  );

  const plotPoints = useMemo(
    () => validSeries.map((point) => ({
      point,
      x: toX(point.x),
      y: toY(Number(point.y)),
    })),
    [toX, toY, validSeries],
  );

  const pathData = useMemo(() => {
    if (plotPoints.length === 0) {
      return "";
    }
    return plotPoints
      .map((point, index) => `${index === 0 ? "M" : "L"} ${point.x.toFixed(2)} ${point.y.toFixed(2)}`)
      .join(" ");
  }, [plotPoints]);

  const xTicks = useMemo(() => {
    const uniqueValues = getUniqueSeriesValues(xValues);
    if (uniqueValues.length > 1 && uniqueValues.length <= 8) {
      return uniqueValues;
    }
    return buildTransportTicks(xMin, xMax, axisSettings.tickDensity);
  }, [axisSettings.tickDensity, xMax, xMin, xValues]);

  const yTicks = useMemo(
    () => buildTransportTicks(yMin, yMax, axisSettings.tickDensity),
    [axisSettings.tickDensity, yMax, yMin],
  );

  const handleMouseMove = useCallback((event: React.MouseEvent<SVGSVGElement>) => {
    if (!svgRef.current || plotPoints.length === 0) {
      return;
    }
    const rect = svgRef.current.getBoundingClientRect();
    const mouseX = event.clientX - rect.left - margin.left;
    const mouseY = event.clientY - rect.top - margin.top;
    if (
      mouseX < 0
      || mouseX > plotWidth
      || mouseY < 0
      || mouseY > plotHeight
    ) {
      setHoveredPoint(null);
      return;
    }

    let closest: HoveredPoint | null = null;
    let bestDistance = Number.POSITIVE_INFINITY;
    plotPoints.forEach((point) => {
      const distance = Math.hypot(point.x - mouseX, point.y - mouseY);
      if (distance < bestDistance) {
        bestDistance = distance;
        closest = point;
      }
    });
    setHoveredPoint(bestDistance <= 24 ? closest : null);
  }, [margin.left, margin.top, plotHeight, plotPoints, plotWidth]);

  const handleMouseLeave = useCallback(() => {
    setHoveredPoint(null);
  }, []);

  const handleWheel = useCallback((event: React.WheelEvent<SVGSVGElement>) => {
    event.preventDefault();
    const rect = svgRef.current?.getBoundingClientRect();
    if (!rect) {
      return;
    }
    const cursorY = clamp(event.clientY - rect.top - margin.top, 0, plotHeight);
    const cursorRatio = 1 - cursorY / plotHeight;
    const cursorValue = yMin + cursorRatio * ySpan;

    if (event.shiftKey) {
      const shift = (event.deltaY / Math.max(plotHeight, 1)) * ySpan * 0.65;
      onCustomYRangeChange([yMin + shift, yMax + shift]);
      return;
    }

    const scale = event.deltaY < 0 ? 0.88 : 1.12;
    const nextMin = cursorValue - (cursorValue - yMin) * scale;
    const nextMax = cursorValue + (yMax - cursorValue) * scale;
    if (Math.abs(nextMax - nextMin) < 1e-12) {
      return;
    }
    onCustomYRangeChange([nextMin, nextMax]);
  }, [margin.top, onCustomYRangeChange, plotHeight, yMax, yMin, ySpan]);

  const tooltipX = hoveredPoint
    ? clamp(margin.left + hoveredPoint.x + 12, margin.left + 8, width - 220)
    : 0;
  const tooltipY = hoveredPoint
    ? clamp(margin.top + hoveredPoint.y - 84, margin.top + 8, height - 96)
    : 0;
  const sliceInsetWidth = Math.max(112, Math.min(210, fixedContextLabel.length * 7.1 + 18));
  const sliceInsetX = plotWidth - sliceInsetWidth - 10;
  const sliceInsetY = plotHeight - 34;

  return (
    <div className="transport-lineplot-shell">
      <div className="transport-lineplot-canvas" ref={canvasRef}>
        <svg
          ref={svgRef}
          width={width}
          height={height}
          viewBox={`0 0 ${width} ${height}`}
          onMouseMove={handleMouseMove}
          onMouseLeave={handleMouseLeave}
          onWheel={handleWheel}
          style={{ cursor: "crosshair", display: "block" }}
        >
          <defs>
            <clipPath id={clipPathId}>
              <rect x={0} y={0} width={plotWidth} height={plotHeight} />
            </clipPath>
          </defs>

          <rect width={width} height={height} fill={colors.bg} />

          <g transform={`translate(${margin.left}, ${margin.top})`}>
            {axisSettings.showGrid && (
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
            )}

            {axisSettings.showZeroLine && yMin <= 0 && yMax >= 0 && (
              <line
                x1={0}
                x2={plotWidth}
                y1={toY(0)}
                y2={toY(0)}
                stroke="#d32f2f"
                strokeDasharray="4,4"
              />
            )}

            <g clipPath={`url(#${clipPathId})`}>
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

              {plotPoints.map((point) => (
                <circle
                  key={`${point.point.index}-${point.point.x}`}
                  cx={point.x}
                  cy={point.y}
                  r={3.6}
                  fill={metricDefinition.accentColor}
                  opacity={0.9}
                />
              ))}

              {hoveredPoint && (
                <>
                  <line
                    x1={hoveredPoint.x}
                    x2={hoveredPoint.x}
                    y1={0}
                    y2={plotHeight}
                    stroke={colors.grid}
                    strokeDasharray="4,4"
                    opacity={0.55}
                  />
                  <line
                    x1={0}
                    x2={plotWidth}
                    y1={hoveredPoint.y}
                    y2={hoveredPoint.y}
                    stroke={colors.grid}
                    strokeDasharray="4,4"
                    opacity={0.55}
                  />
                </>
              )}
            </g>

            <line x1={0} y1={plotHeight} x2={plotWidth} y2={plotHeight} stroke={colors.axis} />
            <line x1={0} y1={0} x2={0} y2={plotHeight} stroke={colors.axis} />

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
                  y={plotHeight + 24}
                  textAnchor="middle"
                  fill={colors.text}
                  fontSize={12}
                >
                  {formatTransportNumber(
                    tick,
                    axisSettings.precisionMode,
                    axisSettings.precision,
                  )}
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
                  fontSize={12}
                >
                  {formatTransportNumber(
                    tick,
                    axisSettings.precisionMode,
                    axisSettings.precision,
                  )}
                </text>
              </g>
            ))}

            <text
              x={plotWidth / 2}
              y={plotHeight + 52}
              textAnchor="middle"
              fill={colors.text}
              fontSize={14}
              fontFamily="serif"
              fontStyle="italic"
            >
              {xAxisLabel}
            </text>
            <text
              x={-58}
              y={plotHeight / 2}
              transform={`rotate(-90 -58 ${plotHeight / 2})`}
              textAnchor="middle"
              dominantBaseline="middle"
              fill={colors.text}
              fontSize={14}
              fontFamily="serif"
              fontStyle="italic"
            >
              {yAxisLabel}
            </text>

            <g pointerEvents="none">
              <rect
                x={sliceInsetX}
                y={sliceInsetY}
                width={sliceInsetWidth}
                height={24}
                rx={6}
                fill={colors.tooltip}
                fillOpacity={0.92}
                stroke={colors.tooltipBorder}
              />
              <text
                x={sliceInsetX + sliceInsetWidth / 2}
                y={sliceInsetY + 16}
                textAnchor="middle"
                fill={colors.tooltipText}
                fontSize={11}
              >
                {fixedContextLabel}
              </text>
            </g>

            {hoveredPoint && (
              <circle
                cx={hoveredPoint.x}
                cy={hoveredPoint.y}
                r={5.5}
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
                width={208}
                height={80}
                fill={colors.tooltip}
                stroke={colors.tooltipBorder}
                rx={4}
              />
              <text x={10} y={18} fill={colors.tooltipText} fontSize={11}>
                {metricDefinition.symbol} · {component}
              </text>
              <text x={10} y={36} fill={colors.tooltipText} fontSize={11}>
                x = {formatTransportNumber(
                  hoveredPoint.point.x,
                  axisSettings.precisionMode,
                  axisSettings.precision,
                )} {axisMode === "temperature" ? "K" : "eV"}
              </text>
              <text x={10} y={52} fill={colors.tooltipText} fontSize={11}>
                y = {formatTransportNumber(
                  hoveredPoint.point.y,
                  axisSettings.precisionMode,
                  axisSettings.precision,
                )} {valueUnit}
              </text>
              <text x={10} y={68} fill={colors.tooltipText} fontSize={11}>
                {fixedContextLabel}
              </text>
            </g>
          )}
        </svg>
      </div>
    </div>
  );
}
