import { useEffect, useMemo, useState } from "react";

interface TransportPlotProps {
  data: any;
}

type TransportMetricKey = "conductivity" | "sigma_s" | "seebeck" | "kappa";
type TransportAxisMode = "mu" | "temperature";
interface TransportSeriesPoint {
  x: number;
  y: number | null;
}

interface TransportPlotPoint {
  x: number;
  y: number;
  rawX: number;
  rawY: number;
}

const METRIC_OPTIONS: Array<{ key: TransportMetricKey; label: string }> = [
  { key: "conductivity", label: "Conductivity" },
  { key: "sigma_s", label: "Sigma·S" },
  { key: "seebeck", label: "Seebeck" },
  { key: "kappa", label: "K (BoltzWann)" },
];

function findComponentIndex(labels: string[], target: string): number {
  return labels.findIndex((label) => label.trim().toLowerCase() === target);
}

function getComponentOptions(dataset: any): string[] {
  const labels = Array.isArray(dataset?.component_labels) ? dataset.component_labels.map((label: unknown) => String(label)) : [];
  const options: string[] = [];
  const xx = findComponentIndex(labels, "xx");
  const yy = findComponentIndex(labels, "yy");
  const zz = findComponentIndex(labels, "zz");
  if (xx >= 0 || yy >= 0 || zz >= 0) {
    options.push("avg");
  }
  for (const preferred of ["xx", "yy", "zz"]) {
    if (findComponentIndex(labels, preferred) >= 0) {
      options.push(preferred);
    }
  }
  for (const label of labels) {
    const normalized = label.trim().toLowerCase();
    if (!options.includes(normalized)) {
      options.push(normalized);
    }
  }
  return options.length > 0 ? options : ["value"];
}

function getDatasetComponentValue(
  dataset: any,
  component: string,
  temperatureIndex: number,
  muIndex: number,
): number | null {
  const labels = Array.isArray(dataset?.component_labels) ? dataset.component_labels.map((label: unknown) => String(label).trim().toLowerCase()) : [];
  const values = Array.isArray(dataset?.values) ? dataset.values : [];

  const readValue = (componentIndex: number): number | null => {
    const raw = values?.[componentIndex]?.[temperatureIndex]?.[muIndex];
    const numeric = Number(raw);
    return Number.isFinite(numeric) ? numeric : null;
  };

  if (component === "avg") {
    const diagonalIndices = ["xx", "yy", "zz"]
      .map((label) => findComponentIndex(labels, label))
      .filter((index) => index >= 0);
    if (diagonalIndices.length === 0) {
      return null;
    }
    const diagonalValues = diagonalIndices
      .map(readValue)
      .filter((value): value is number => value != null);
    if (diagonalValues.length === 0) {
      return null;
    }
    return diagonalValues.reduce((sum, value) => sum + value, 0) / diagonalValues.length;
  }

  const index = findComponentIndex(labels, component);
  if (index < 0) {
    return null;
  }
  return readValue(index);
}

function formatNumber(value: number | null | undefined, digits = 4): string {
  if (!Number.isFinite(Number(value))) {
    return "N/A";
  }
  return Number(value).toFixed(digits);
}

export function TransportPlot({ data }: TransportPlotProps) {
  const [metric, setMetric] = useState<TransportMetricKey>("conductivity");
  const [component, setComponent] = useState("avg");
  const [axisMode, setAxisMode] = useState<TransportAxisMode>("mu");
  const [selectedTemperatureIndex, setSelectedTemperatureIndex] = useState(0);
  const [selectedMuIndex, setSelectedMuIndex] = useState(0);

  const dataset = data?.[metric] ?? null;
  const componentOptions = useMemo(() => getComponentOptions(dataset), [dataset]);
  const muOffsets = Array.isArray(data?.mu_offsets_ev) ? data.mu_offsets_ev.map((value: unknown) => Number(value)) : [];
  const temperatureValues = Array.isArray(data?.temperature_values_k)
    ? data.temperature_values_k.map((value: unknown) => Number(value))
    : [];

  useEffect(() => {
    if (!componentOptions.includes(component)) {
      setComponent(componentOptions[0] ?? "value");
    }
  }, [component, componentOptions]);

  useEffect(() => {
    if (!data?.[metric]) {
      const fallbackMetric = METRIC_OPTIONS.find((option) => data?.[option.key] != null)?.key ?? "conductivity";
      setMetric(fallbackMetric);
    }
  }, [data, metric]);

  const safeTemperatureIndex = Math.min(
    Math.max(0, selectedTemperatureIndex),
    Math.max(0, temperatureValues.length - 1),
  );
  const safeMuIndex = Math.min(
    Math.max(0, selectedMuIndex),
    Math.max(0, muOffsets.length - 1),
  );

  const series = useMemo<TransportSeriesPoint[]>(() => {
    if (!dataset) return [];

    if (axisMode === "mu") {
      return muOffsets.map((muOffset: number, index: number) => ({
        x: muOffset,
        y: getDatasetComponentValue(dataset, component, safeTemperatureIndex, index),
      }));
    }

    return temperatureValues.map((temperature: number, index: number) => ({
      x: temperature,
      y: getDatasetComponentValue(dataset, component, index, safeMuIndex),
    }));
  }, [axisMode, component, dataset, muOffsets, safeMuIndex, safeTemperatureIndex, temperatureValues]);

  const validSeries = series.filter((point: TransportSeriesPoint) => point.y != null && Number.isFinite(point.x));
  const xValues = validSeries.map((point: TransportSeriesPoint) => point.x);
  const yValues = validSeries.map((point: TransportSeriesPoint) => Number(point.y));
  const xMin = xValues.length > 0 ? Math.min(...xValues) : 0;
  const xMax = xValues.length > 0 ? Math.max(...xValues) : 1;
  const yMinRaw = yValues.length > 0 ? Math.min(...yValues) : 0;
  const yMaxRaw = yValues.length > 0 ? Math.max(...yValues) : 1;
  const ySpan = Math.max(1e-12, yMaxRaw - yMinRaw);
  const yPad = ySpan * 0.1;
  const yMin = yMinRaw - yPad;
  const yMax = yMaxRaw + yPad;
  const width = 860;
  const height = 360;
  const paddingLeft = 68;
  const paddingRight = 24;
  const paddingTop = 28;
  const paddingBottom = 48;

  const plotWidth = width - paddingLeft - paddingRight;
  const plotHeight = height - paddingTop - paddingBottom;
  const xDenominator = Math.max(1e-12, xMax - xMin);
  const yDenominator = Math.max(1e-12, yMax - yMin);

  const points: TransportPlotPoint[] = validSeries.map((point: TransportSeriesPoint) => {
    const x = paddingLeft + ((point.x - xMin) / xDenominator) * plotWidth;
    const y = paddingTop + (1 - ((Number(point.y) - yMin) / yDenominator)) * plotHeight;
    return { x, y, rawX: point.x, rawY: Number(point.y) };
  });

  const polylinePoints = points.map((point: TransportPlotPoint) => `${point.x},${point.y}`).join(" ");
  const axisLabel = axisMode === "mu" ? "Δμ (eV)" : "Temperature (K)";
  const fixedLabel = axisMode === "mu"
    ? `T = ${formatNumber(temperatureValues[safeTemperatureIndex], 1)} K`
    : `Δμ = ${formatNumber(muOffsets[safeMuIndex], 3)} eV`;

  if (!dataset) {
    return <p>No transport data is available for plotting.</p>;
  }

  return (
    <div className="bands-viewer-content" style={{ display: "block" }}>
      <div className="calc-action-grid" style={{ marginBottom: "1rem" }}>
        {METRIC_OPTIONS.map((option) => (
          <button
            key={option.key}
            type="button"
            className="calc-action-btn"
            onClick={() => setMetric(option.key)}
            style={metric === option.key ? { borderColor: "var(--accent-color, #0f766e)" } : undefined}
          >
            <span className="calc-action-label">{option.label}</span>
            <span className="calc-action-hint">{option.key === "seebeck" ? "tau-independent" : "tau-dependent"}</span>
          </button>
        ))}
      </div>

      <div className="details-grid" style={{ marginBottom: "1rem" }}>
        <div className="detail-item">
          <label>Component</label>
          <select value={component} onChange={(event) => setComponent(event.target.value)}>
            {componentOptions.map((option) => (
              <option key={option} value={option}>{option}</option>
            ))}
          </select>
        </div>
        <div className="detail-item">
          <label>Axis</label>
          <select value={axisMode} onChange={(event) => setAxisMode(event.target.value as TransportAxisMode)}>
            <option value="mu">vs Δμ</option>
            <option value="temperature">vs T</option>
          </select>
        </div>
        {axisMode === "mu" && (
          <div className="detail-item">
            <label>Temperature Slice</label>
            <select
              value={safeTemperatureIndex}
              onChange={(event) => setSelectedTemperatureIndex(Number(event.target.value))}
            >
              {temperatureValues.map((value: number, index: number) => (
                <option key={`${value}-${index}`} value={index}>
                  {formatNumber(value, 1)} K
                </option>
              ))}
            </select>
          </div>
        )}
        {axisMode === "temperature" && (
          <div className="detail-item">
            <label>Δμ Slice</label>
            <select
              value={safeMuIndex}
              onChange={(event) => setSelectedMuIndex(Number(event.target.value))}
            >
              {muOffsets.map((value: number, index: number) => (
                <option key={`${value}-${index}`} value={index}>
                  {formatNumber(value, 3)} eV
                </option>
              ))}
            </select>
          </div>
        )}
        <div className="detail-item">
          <label>Slice</label>
          <span>{fixedLabel}</span>
        </div>
      </div>

      <svg viewBox={`0 0 ${width} ${height}`} style={{ width: "100%", maxWidth: "100%", border: "1px solid rgba(15, 23, 42, 0.08)", borderRadius: "14px", background: "rgba(248, 250, 252, 0.72)" }}>
        <rect x="0" y="0" width={width} height={height} fill="transparent" />
        <line x1={paddingLeft} y1={height - paddingBottom} x2={width - paddingRight} y2={height - paddingBottom} stroke="rgba(15,23,42,0.35)" />
        <line x1={paddingLeft} y1={paddingTop} x2={paddingLeft} y2={height - paddingBottom} stroke="rgba(15,23,42,0.35)" />
        <text x={width / 2} y={height - 12} textAnchor="middle" fontSize="13" fill="rgba(15,23,42,0.75)">
          {axisLabel}
        </text>
        <text x="18" y={height / 2} textAnchor="middle" fontSize="13" fill="rgba(15,23,42,0.75)" transform={`rotate(-90 18 ${height / 2})`}>
          {METRIC_OPTIONS.find((option) => option.key === metric)?.label || metric}
        </text>

        {points.length > 1 && (
          <polyline
            fill="none"
            stroke="var(--accent-color, #0f766e)"
            strokeWidth="3"
            strokeLinejoin="round"
            strokeLinecap="round"
            points={polylinePoints}
          />
        )}
        {points.map((point: TransportPlotPoint, index: number) => (
          <g key={`${point.rawX}-${index}`}>
            <circle cx={point.x} cy={point.y} r="4" fill="var(--accent-color, #0f766e)" />
            <title>{`${axisLabel}: ${formatNumber(point.rawX, axisMode === "mu" ? 3 : 1)}, value: ${formatNumber(point.rawY, 6)}`}</title>
          </g>
        ))}
      </svg>

      <div className="details-grid" style={{ marginTop: "1rem" }}>
        <div className="detail-item">
          <label>Engine</label>
          <span>{data?.engine || "boltzwann"}</span>
        </div>
        <div className="detail-item">
          <label>Reference Fermi Energy</label>
          <span>{formatNumber(data?.reference_fermi_energy_ev, 4)} eV</span>
        </div>
        <div className="detail-item">
          <label>Relaxation Time</label>
          <span>{formatNumber(data?.relaxation_time_fs, 2)} fs</span>
        </div>
        <div className="detail-item">
          <label>Provenance</label>
          <span>{metric === "seebeck" ? "Seebeck is tau-independent" : "Conductivity, Sigma·S, and K are tau-dependent"}</span>
        </div>
      </div>

      <p className="step-description" style={{ marginTop: "1rem" }}>
        `K (BoltzWann)` is displayed with its upstream BoltzWann meaning and is not silently relabeled as electronic thermal conductivity.
      </p>
    </div>
  );
}
