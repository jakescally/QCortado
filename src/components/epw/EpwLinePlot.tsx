import { useMemo } from "react";
import { EpwSeriesPoint, formatEpwNumber } from "../../lib/epw";

interface EpwLinePlotProps {
  title: string;
  xLabel: string;
  yLabel: string;
  series: EpwSeriesPoint[];
  color?: string;
}

interface EpwCompactHeatmapProps {
  title: string;
  columnLabels: string[];
  rowLabels: string[];
  values: Array<Array<number | null>>;
  formatter?: (value: number) => string;
}

function finiteValues(series: EpwSeriesPoint[], key: "x" | "y"): number[] {
  return series
    .map((point) => point[key])
    .filter((value): value is number => typeof value === "number" && Number.isFinite(value));
}

function scale(value: number, min: number, max: number, size: number): number {
  const span = Math.abs(max - min) < 1e-12 ? 1 : max - min;
  return ((value - min) / span) * size;
}

export function EpwLinePlot({
  title,
  xLabel,
  yLabel,
  series,
  color = "#2563eb",
}: EpwLinePlotProps) {
  const validSeries = useMemo(
    () => series.filter((point) => Number.isFinite(point.x) && point.y != null && Number.isFinite(point.y)),
    [series],
  );

  const width = 720;
  const height = 260;
  const margin = { top: 28, right: 24, bottom: 48, left: 72 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  const xValues = finiteValues(validSeries, "x");
  const yValues = finiteValues(validSeries, "y");
  const xMin = xValues.length ? Math.min(...xValues) : 0;
  const xMax = xValues.length ? Math.max(...xValues) : 1;
  const yMinRaw = yValues.length ? Math.min(...yValues) : 0;
  const yMaxRaw = yValues.length ? Math.max(...yValues) : 1;
  const yPad = Math.max(Math.abs(yMaxRaw - yMinRaw) * 0.1, Math.abs(yMaxRaw) * 0.02, 1e-12);
  const yMin = yMinRaw - yPad;
  const yMax = yMaxRaw + yPad;

  const plotPoints = validSeries.map((point) => ({
    point,
    x: scale(point.x, xMin, xMax, plotWidth),
    y: plotHeight - scale(Number(point.y), yMin, yMax, plotHeight),
  }));
  const pathData = plotPoints
    .map((point, index) => `${index === 0 ? "M" : "L"} ${point.x.toFixed(2)} ${point.y.toFixed(2)}`)
    .join(" ");
  const xTicks = [xMin, (xMin + xMax) / 2, xMax].filter((value, index, values) => (
    index === 0 || Math.abs(value - values[index - 1]) > 1e-12
  ));
  const yTicks = [yMinRaw, (yMinRaw + yMaxRaw) / 2, yMaxRaw].filter((value, index, values) => (
    index === 0 || Math.abs(value - values[index - 1]) > 1e-12
  ));

  if (validSeries.length === 0) {
    return (
      <div className="epw-plot-shell epw-empty-plot">
        <span>{title}</span>
        <p>No finite series values were parsed for this view.</p>
      </div>
    );
  }

  return (
    <div className="epw-plot-shell">
      <svg viewBox={`0 0 ${width} ${height}`} role="img" aria-label={title}>
        <rect x={0} y={0} width={width} height={height} className="epw-plot-bg" />
        <text x={margin.left} y={18} className="epw-plot-title">{title}</text>
        <g transform={`translate(${margin.left}, ${margin.top})`}>
          {yTicks.map((tick) => {
            const y = plotHeight - scale(tick, yMin, yMax, plotHeight);
            return (
              <g key={`y-${tick}`}>
                <line x1={0} x2={plotWidth} y1={y} y2={y} className="epw-plot-grid" />
                <text x={-10} y={y + 4} textAnchor="end" className="epw-plot-tick">
                  {formatEpwNumber(tick, 3)}
                </text>
              </g>
            );
          })}
          {xTicks.map((tick) => {
            const x = scale(tick, xMin, xMax, plotWidth);
            return (
              <g key={`x-${tick}`}>
                <line x1={x} x2={x} y1={0} y2={plotHeight} className="epw-plot-grid epw-plot-grid-soft" />
                <text x={x} y={plotHeight + 20} textAnchor="middle" className="epw-plot-tick">
                  {formatEpwNumber(tick, 3)}
                </text>
              </g>
            );
          })}
          <line x1={0} x2={plotWidth} y1={plotHeight} y2={plotHeight} className="epw-plot-axis" />
          <line x1={0} x2={0} y1={0} y2={plotHeight} className="epw-plot-axis" />
          <path d={pathData} fill="none" stroke={color} strokeWidth={2.4} strokeLinecap="round" strokeLinejoin="round" />
          {plotPoints.map((point) => (
            <circle key={point.point.index} cx={point.x} cy={point.y} r={3.5} fill={color}>
              <title>
                {`${xLabel}: ${formatEpwNumber(point.point.x, 4)}; ${yLabel}: ${formatEpwNumber(point.point.y, 4)}`}
              </title>
            </circle>
          ))}
          <text x={plotWidth / 2} y={plotHeight + 40} textAnchor="middle" className="epw-plot-label">{xLabel}</text>
          <text
            x={-plotHeight / 2}
            y={-54}
            transform="rotate(-90)"
            textAnchor="middle"
            className="epw-plot-label"
          >
            {yLabel}
          </text>
        </g>
      </svg>
    </div>
  );
}

export function EpwCompactHeatmap({
  title,
  columnLabels,
  rowLabels,
  values,
  formatter = (value) => formatEpwNumber(value, 3),
}: EpwCompactHeatmapProps) {
  const finite = values.flat().filter((value): value is number => typeof value === "number" && Number.isFinite(value));
  const min = finite.length ? Math.min(...finite) : 0;
  const max = finite.length ? Math.max(...finite) : 1;
  const span = Math.abs(max - min) < 1e-12 ? 1 : max - min;
  return (
    <div className="epw-heatmap">
      <div className="epw-heatmap-title">{title}</div>
      <div className="epw-heatmap-grid" style={{ gridTemplateColumns: `minmax(48px, auto) repeat(${columnLabels.length}, minmax(42px, 1fr))` }}>
        <span className="epw-heatmap-head" />
        {columnLabels.map((label) => <span key={label} className="epw-heatmap-head">{label}</span>)}
        {values.map((row, rowIndex) => (
          <div className="epw-heatmap-row" key={rowLabels[rowIndex] || rowIndex} style={{ display: "contents" }}>
            <span className="epw-heatmap-head epw-heatmap-row-label">{rowLabels[rowIndex] || `r${rowIndex + 1}`}</span>
            {row.map((value, columnIndex) => {
              const ratio = value == null ? 0 : (value - min) / span;
              const opacity = value == null ? 0.08 : 0.18 + ratio * 0.7;
              return (
                <span
                  key={`${rowIndex}-${columnIndex}`}
                  className="epw-heatmap-cell"
                  style={{ backgroundColor: `rgba(37, 99, 235, ${opacity})` }}
                  title={value == null ? "No value" : formatter(value)}
                >
                  {value == null ? "-" : formatter(value)}
                </span>
              );
            })}
          </div>
        ))}
      </div>
    </div>
  );
}
