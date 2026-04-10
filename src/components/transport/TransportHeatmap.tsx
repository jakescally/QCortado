import { useLayoutEffect, useMemo, useRef, useState } from "react";
import { useTheme } from "../../lib/ThemeContext";
import {
  TransportAxisMode,
  TransportAxisSettings,
  TransportHeatmapCell,
  TransportMetricDefinition,
  formatTransportNumber,
} from "../../lib/transport";

function clamp(value: number, min: number, max: number): number {
  return Math.min(Math.max(value, min), max);
}

interface TransportHeatmapProps {
  metricDefinition: TransportMetricDefinition;
  axisMode: TransportAxisMode;
  axisSettings: TransportAxisSettings;
  grid: TransportHeatmapCell[][];
  selectedTemperatureIndex: number;
  selectedMuIndex: number;
  onSelectCell: (temperatureIndex: number, muIndex: number) => void;
}

function hexToRgb(hex: string): [number, number, number] {
  const normalized = hex.replace("#", "");
  const expanded = normalized.length === 3
    ? normalized.split("").map((digit) => `${digit}${digit}`).join("")
    : normalized;
  const safe = expanded.padEnd(6, "0").slice(0, 6);
  return [
    Number.parseInt(safe.slice(0, 2), 16),
    Number.parseInt(safe.slice(2, 4), 16),
    Number.parseInt(safe.slice(4, 6), 16),
  ];
}

function interpolateColor(
  start: [number, number, number],
  end: [number, number, number],
  ratio: number,
): string {
  const t = clamp(ratio, 0, 1);
  const [sr, sg, sb] = start;
  const [er, eg, eb] = end;
  const r = Math.round(sr + (er - sr) * t);
  const g = Math.round(sg + (eg - sg) * t);
  const b = Math.round(sb + (eb - sb) * t);
  return `rgb(${r}, ${g}, ${b})`;
}

function getSampleIndices(length: number): number[] {
  if (length <= 1) {
    return [0];
  }
  if (length <= 4) {
    return Array.from({ length }, (_, index) => index);
  }
  const indices = new Set<number>([
    0,
    Math.floor((length - 1) / 3),
    Math.floor((2 * (length - 1)) / 3),
    length - 1,
  ]);
  return Array.from(indices).sort((left, right) => left - right);
}

export function TransportHeatmap({
  metricDefinition,
  axisMode,
  axisSettings,
  grid,
  selectedTemperatureIndex,
  selectedMuIndex,
  onSelectCell,
}: TransportHeatmapProps) {
  const { isDark } = useTheme();
  const canvasRef = useRef<HTMLDivElement>(null);
  const [containerWidth, setContainerWidth] = useState<number | null>(null);
  const [hoveredCell, setHoveredCell] = useState<{ temperatureIndex: number; muIndex: number } | null>(null);

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
      setContainerWidth((previous) => (previous === width ? previous : width));
    });
    observer.observe(node);
    return () => observer.disconnect();
  }, []);

  const width = containerWidth ?? 332;
  const height = 320;
  const margin = { top: 20, right: 18, bottom: 56, left: 66 };
  const plotWidth = Math.max(1, width - margin.left - margin.right);
  const plotHeight = Math.max(1, height - margin.top - margin.bottom);
  const rowCount = grid.length;
  const columnCount = grid[0]?.length ?? 0;
  const cellWidth = columnCount > 0 ? plotWidth / columnCount : plotWidth;
  const cellHeight = rowCount > 0 ? plotHeight / rowCount : plotHeight;
  const hoverGrow = Math.max(2, Math.min(5, Math.min(cellWidth, cellHeight) * 0.15));
  const hoverGrowHalf = hoverGrow / 2;
  const firstTemperature = grid[0]?.[0]?.temperatureK ?? 0;
  const lastTemperature = grid[rowCount - 1]?.[0]?.temperatureK ?? firstTemperature;
  const reverseTemperatureAxis = rowCount > 1 && firstTemperature <= lastTemperature;

  const colors = useMemo(() => (
    isDark
      ? {
          bg: "#1e1e2e",
          axis: "#718096",
          text: "#e2e8f0",
          empty: "#252738",
        }
      : {
          bg: "#ffffff",
          axis: "#333333",
          text: "#000000",
          empty: "#edf2f7",
        }
  ), [isDark]);

  const validValues = useMemo(() => (
    grid.flatMap((row) => row.map((cell) => cell.value).filter((value): value is number => value != null))
  ), [grid]);
  const minValue = validValues.length > 0 ? Math.min(...validValues) : 0;
  const maxValue = validValues.length > 0 ? Math.max(...validValues) : 1;
  const hasDivergingScale = minValue < 0 && maxValue > 0;
  const accentRgb = hexToRgb(metricDefinition.accentColor);
  const negativeRgb: [number, number, number] = isDark ? [77, 171, 247] : [37, 99, 235];
  const neutralRgb: [number, number, number] = isDark ? [39, 43, 60] : [247, 250, 252];

  const xTickIndices = getSampleIndices(columnCount);
  const yTickIndices = getSampleIndices(rowCount);

  const displayRowIndex = (temperatureIndex: number): number => (
    reverseTemperatureAxis ? rowCount - 1 - temperatureIndex : temperatureIndex
  );
  const rowY = (temperatureIndex: number): number => displayRowIndex(temperatureIndex) * cellHeight;
  const selectedBandTemperatureIndex = axisMode === "mu" ? selectedTemperatureIndex : null;
  const selectedBandMuIndex = axisMode === "temperature" ? selectedMuIndex : null;
  const hoveredBandTemperatureIndex = axisMode === "mu" ? hoveredCell?.temperatureIndex ?? null : null;
  const hoveredBandMuIndex = axisMode === "temperature" ? hoveredCell?.muIndex ?? null : null;

  const getCellFill = (value: number | null): string => {
    if (value == null) {
      return colors.empty;
    }
    if (hasDivergingScale) {
      if (value >= 0) {
        const ratio = maxValue === 0 ? 1 : value / maxValue;
        return interpolateColor(neutralRgb, accentRgb, ratio);
      }
      const ratio = minValue === 0 ? 1 : value / minValue;
      return interpolateColor(neutralRgb, negativeRgb, ratio);
    }
    const span = Math.abs(maxValue - minValue) < 1e-12 ? 1 : maxValue - minValue;
    return interpolateColor(neutralRgb, accentRgb, (value - minValue) / span);
  };

  return (
    <div className="transport-heatmap-shell">
      <div className="transport-heatmap-canvas" ref={canvasRef}>
        <svg width={width} height={height} viewBox={`0 0 ${width} ${height}`} style={{ display: "block" }}>
          <rect x={0} y={0} width={width} height={height} fill={colors.bg} />

          <g transform={`translate(${margin.left}, ${margin.top})`}>
            {selectedBandTemperatureIndex != null && rowCount > 0 && (
              <>
                <rect
                  x={0}
                  y={rowY(selectedBandTemperatureIndex)}
                  width={plotWidth}
                  height={cellHeight}
                  fill={metricDefinition.accentColor}
                  opacity={0.1}
                />
              </>
            )}
            {selectedBandMuIndex != null && columnCount > 0 && (
              <>
                <rect
                  x={selectedBandMuIndex * cellWidth}
                  y={0}
                  width={cellWidth}
                  height={plotHeight}
                  fill={metricDefinition.accentColor}
                  opacity={0.1}
                />
              </>
            )}

            {grid.map((row, rowIndex) => row.map((cell, columnIndex) => {
              const x = columnIndex * cellWidth;
              const y = rowY(rowIndex);
              const isInHoveredSelectableBand = hoveredCell == null
                ? false
                : axisMode === "mu"
                  ? hoveredCell.temperatureIndex === rowIndex
                  : hoveredCell.muIndex === columnIndex;
              let cellX = x;
              let cellY = y;
              let resolvedCellWidth = Math.max(1, cellWidth);
              let resolvedCellHeight = Math.max(1, cellHeight);

              if (isInHoveredSelectableBand) {
                if (axisMode === "mu") {
                  const scaleX = (plotWidth + hoverGrow) / Math.max(plotWidth, 1);
                  const scaleY = (cellHeight + hoverGrow) / Math.max(cellHeight, 1);
                  const bandCenterX = plotWidth / 2;
                  const bandCenterY = y + cellHeight / 2;
                  cellX = bandCenterX + (x - bandCenterX) * scaleX;
                  cellY = bandCenterY + (y - bandCenterY) * scaleY;
                  resolvedCellWidth = Math.max(1, cellWidth * scaleX);
                  resolvedCellHeight = Math.max(1, cellHeight * scaleY);
                } else {
                  const scaleX = (cellWidth + hoverGrow) / Math.max(cellWidth, 1);
                  const scaleY = (plotHeight + hoverGrow) / Math.max(plotHeight, 1);
                  const bandCenterX = x + cellWidth / 2;
                  const bandCenterY = plotHeight / 2;
                  cellX = bandCenterX + (x - bandCenterX) * scaleX;
                  cellY = bandCenterY + (y - bandCenterY) * scaleY;
                  resolvedCellWidth = Math.max(1, cellWidth * scaleX);
                  resolvedCellHeight = Math.max(1, cellHeight * scaleY);
                }
              }
              return (
                <g key={`${rowIndex}-${columnIndex}`}>
                  <rect
                    className={`transport-heatmap-cell ${isInHoveredSelectableBand ? "is-hover-band" : ""}`}
                    x={cellX}
                    y={cellY}
                    width={resolvedCellWidth}
                    height={resolvedCellHeight}
                    fill={getCellFill(cell.value)}
                    stroke={colors.bg}
                    strokeWidth={1}
                    onMouseEnter={() => setHoveredCell({ temperatureIndex: rowIndex, muIndex: columnIndex })}
                    onMouseLeave={() => setHoveredCell(null)}
                    onClick={() => onSelectCell(rowIndex, columnIndex)}
                    style={{ cursor: "pointer" }}
                  />
                </g>
              );
            }))}

            {selectedBandTemperatureIndex != null && (
              <rect
                className="transport-heatmap-selected-band"
                x={0}
                y={rowY(selectedBandTemperatureIndex)}
                width={plotWidth}
                height={cellHeight}
                fill="none"
                stroke={metricDefinition.accentColor}
                strokeWidth={2.4}
                pointerEvents="none"
              />
            )}
            {selectedBandMuIndex != null && (
              <rect
                className="transport-heatmap-selected-band"
                x={selectedBandMuIndex * cellWidth}
                y={0}
                width={cellWidth}
                height={plotHeight}
                fill="none"
                stroke={metricDefinition.accentColor}
                strokeWidth={2.4}
                pointerEvents="none"
              />
            )}

            {hoveredBandTemperatureIndex != null && (
              <rect
                className="transport-heatmap-hover-band"
                x={-hoverGrowHalf}
                y={rowY(hoveredBandTemperatureIndex) - hoverGrowHalf}
                width={plotWidth + hoverGrow}
                height={cellHeight + hoverGrow}
                fill={metricDefinition.accentColor}
                fillOpacity={0.06}
                pointerEvents="none"
              />
            )}
            {hoveredBandMuIndex != null && (
              <rect
                className="transport-heatmap-hover-band"
                x={(hoveredBandMuIndex * cellWidth) - hoverGrowHalf}
                y={-hoverGrowHalf}
                width={cellWidth + hoverGrow}
                height={plotHeight + hoverGrow}
                fill={metricDefinition.accentColor}
                fillOpacity={0.06}
                pointerEvents="none"
              />
            )}

            <line x1={0} x2={plotWidth} y1={plotHeight} y2={plotHeight} stroke={colors.axis} />
            <line x1={0} x2={0} y1={0} y2={plotHeight} stroke={colors.axis} />

            {xTickIndices.map((index) => {
              const cell = grid[0]?.[index];
              if (!cell) {
                return null;
              }
              const x = index * cellWidth + cellWidth / 2;
              const tickValue = axisSettings.showAbsoluteMu ? cell.muValueEv : cell.muOffsetEv;
              return (
                <g key={`x-${index}`}>
                  <line
                    x1={x}
                    x2={x}
                    y1={plotHeight}
                    y2={plotHeight + 6}
                    stroke={colors.axis}
                  />
                  <text
                    x={x}
                    y={plotHeight + 22}
                    textAnchor="middle"
                    fill={colors.text}
                    fontSize={11}
                  >
                    {formatTransportNumber(
                      tickValue,
                      axisSettings.precisionMode,
                      axisSettings.precision,
                    )}
                  </text>
                </g>
              );
            })}

            {yTickIndices.map((index) => {
              const cell = grid[index]?.[0];
              if (!cell) {
                return null;
              }
              const y = rowY(index) + cellHeight / 2;
              return (
                <g key={`y-${index}`}>
                  <line
                    x1={-6}
                    x2={0}
                    y1={y}
                    y2={y}
                    stroke={colors.axis}
                  />
                  <text
                    x={-10}
                    y={y + 4}
                    textAnchor="end"
                    fill={colors.text}
                    fontSize={11}
                  >
                    {formatTransportNumber(
                      cell.temperatureK,
                      axisSettings.precisionMode,
                      axisSettings.precision,
                    )}
                  </text>
                </g>
              );
            })}

            <text
              x={plotWidth / 2}
              y={plotHeight + 46}
              textAnchor="middle"
              fill={colors.text}
              fontSize={13}
              fontFamily="serif"
              fontStyle="italic"
            >
              {axisSettings.showAbsoluteMu ? "μ (eV)" : "Δμ (eV)"}
            </text>
            <text
              x={-50}
              y={plotHeight / 2}
              textAnchor="middle"
              fill={colors.text}
              fontSize={13}
              fontFamily="serif"
              fontStyle="italic"
              transform={`rotate(-90 -50 ${plotHeight / 2})`}
            >
              T (K)
            </text>
          </g>
        </svg>
      </div>
    </div>
  );
}
