import { useCallback, useMemo } from "react";
import { InfoTooltip } from "./InfoTooltip";
import { TransportHeatmap } from "./transport/TransportHeatmap";
import { TransportLinePlot } from "./transport/TransportLinePlot";
import { TransportTdfPlot } from "./transport/TransportTdfPlot";
import { useTransportViewerState } from "./transport/useTransportViewerState";
import {
  TRANSPORT_METRIC_KEYS,
  TransportResult,
  buildTransportHeatmapGrid,
  buildTransportSeries,
  formatTransportNumber,
  getTransportMetricDefinition,
  getTransportTauBadge,
  resolveTransportUnitSuffix,
  resolveTransportXAxisLabel,
  resolveTransportYAxisLabel,
} from "../lib/transport";

interface TransportPlotProps {
  data: TransportResult;
}

function formatBytes(bytes: number): string {
  if (!Number.isFinite(bytes) || bytes <= 0) {
    return "0 B";
  }
  const units = ["B", "KB", "MB", "GB", "TB"];
  let value = bytes;
  let unitIndex = 0;
  while (value >= 1024 && unitIndex < units.length - 1) {
    value /= 1024;
    unitIndex += 1;
  }
  return `${value.toFixed(value >= 10 || unitIndex === 0 ? 0 : 1)} ${units[unitIndex]}`;
}

function joinTransportTooltip(base: string, extra?: string): string {
  return extra ? `${base} ${extra}` : base;
}

export function TransportPlot({ data }: TransportPlotProps) {
  const {
    state,
    dataset,
    componentOptions,
    tdfOptions,
    safeTemperatureIndex,
    safeMuIndex,
    rangeError,
    setMetric,
    setComponent,
    setAxisMode,
    setSelectedTemperatureIndex,
    setSelectedMuIndex,
    applyHeatmapSelection,
    setCustomYRange,
    updateAxisSettings,
    applyManualYRange,
    resetView,
    resetAxisSettings,
    resetAll,
    toggleHeatmapOpen,
    toggleSettingsOpen,
    toggleRunContextOpen,
    toggleTdfOpen,
    setTdfComponent,
  } = useTransportViewerState(data);

  const metricDefinition = getTransportMetricDefinition(state.metric);
  const tdfDefinition = getTransportMetricDefinition("tdf");
  const yAxisLabel = resolveTransportYAxisLabel(metricDefinition, state.axisSettings);
  const xAxisLabel = resolveTransportXAxisLabel(state.axisMode, state.axisSettings);
  const valueUnit = resolveTransportUnitSuffix(metricDefinition, state.axisSettings);

  const series = useMemo(() => buildTransportSeries({
    result: data,
    dataset,
    component: state.component,
    axisMode: state.axisMode,
    selectedTemperatureIndex: safeTemperatureIndex,
    selectedMuIndex: safeMuIndex,
    showAbsoluteMu: state.axisSettings.showAbsoluteMu,
  }), [
    data,
    dataset,
    safeMuIndex,
    safeTemperatureIndex,
    state.axisMode,
    state.axisSettings.showAbsoluteMu,
    state.component,
  ]);

  const heatmapGrid = useMemo(
    () => buildTransportHeatmapGrid(data, dataset, state.component),
    [data, dataset, state.component],
  );

  const selectedTemperature = data.temperature_values_k[safeTemperatureIndex] ?? 0;
  const selectedMuOffset = data.mu_offsets_ev[safeMuIndex] ?? 0;
  const selectedMuValue = data.mu_values_ev[safeMuIndex] ?? selectedMuOffset;
  const fixedContextLabel = state.axisMode === "mu"
    ? `T = ${formatTransportNumber(
        selectedTemperature,
        state.axisSettings.precisionMode,
        state.axisSettings.precision,
      )} K`
    : `${state.axisSettings.showAbsoluteMu ? "μ" : "Δμ"} = ${formatTransportNumber(
        state.axisSettings.showAbsoluteMu ? selectedMuValue : selectedMuOffset,
        state.axisSettings.precisionMode,
        state.axisSettings.precision,
      )} eV`;
  const activeComponent = componentOptions.find((option) => option.value === state.component);

  const artifactBytes = Array.isArray(data.artifact_manifest)
    ? data.artifact_manifest.reduce((sum, item) => sum + (Number(item?.size_bytes) || 0), 0)
    : 0;

  const manualRangeKeyDown = useCallback((event: React.KeyboardEvent<HTMLInputElement>) => {
    if (event.key !== "Enter") {
      return;
    }
    event.preventDefault();
    applyManualYRange();
  }, [applyManualYRange]);

  const metricTooltip = joinTransportTooltip(
    metricDefinition.tooltip,
    metricDefinition.tauBehavior === "dependent"
      ? "Displayed values scale with the transport relaxation time."
      : metricDefinition.tauBehavior === "independent"
        ? "Displayed values are treated as relaxation-time independent in this workflow."
        : undefined,
  );

  return (
    <div className="band-plot-layout transport-plot-layout">
      <div className="band-plot-main transport-plot-main">
        <div className="transport-plot-header">
          <div className="transport-plot-header-main">
            <div className="transport-plot-header-copy">
              <span className="transport-plot-kicker">BoltzWann Transport</span>
              <h3>{metricDefinition.symbol} · {activeComponent?.label ?? state.component}</h3>
            </div>
          </div>

          <div className="transport-plot-header-side">
            <div className="transport-plot-header-controls">
              <div className="band-control-row transport-plot-control">
                <label>
                  Component
                  <InfoTooltip
                    text={activeComponent?.tooltip ?? "Choose the tensor component or the diagonal average."}
                  />
                </label>
                <select
                  value={state.component}
                  onChange={(event) => setComponent(event.target.value)}
                >
                  {componentOptions.map((option) => (
                    <option key={option.value} value={option.value}>{option.label}</option>
                  ))}
                </select>
              </div>

              <div className="band-control-row transport-plot-control">
                <label>
                  Plot Against
                  <InfoTooltip text="Plot the selected metric against chemical potential at fixed temperature, or against temperature at fixed chemical potential." />
                </label>
                <select
                  value={state.axisMode}
                  onChange={(event) => setAxisMode(event.target.value as typeof state.axisMode)}
                >
                  <option value="mu">Chemical potential</option>
                  <option value="temperature">Temperature</option>
                </select>
              </div>
            </div>
          </div>
        </div>

        <div className="transport-lineplot-stage">
          <TransportLinePlot
            lockScopeKey={`${data.source_wannier_calc_id || data.seedname}:${state.metric}:${state.component}:${state.axisMode}`}
            metricDefinition={metricDefinition}
            axisMode={state.axisMode}
            axisSettings={state.axisSettings}
            component={state.component}
            series={series}
            xAxisLabel={xAxisLabel}
            yAxisLabel={yAxisLabel}
            valueUnit={valueUnit}
            fixedContextLabel={fixedContextLabel}
            onCustomYRangeChange={setCustomYRange}
          />
        </div>

        <div className="band-plot-info transport-plot-info">
          <span>seedname {data.seedname}</span>
          <span>E_F = {formatTransportNumber(data.reference_fermi_energy_ev, "manual", 4)} eV</span>
          <span>τ = {formatTransportNumber(data.relaxation_time_fs, "manual", 2)} fs</span>
          <span>{data.mu_values_ev.length} μ points</span>
          <span>{data.temperature_values_k.length} T points</span>
          {data.is_2d && (
            <span className="band-plot-projection-pill">
              2D · {data.boltz_2d_dir?.trim() || "enabled"}
            </span>
          )}
          {data.warnings.length > 0 && (
            <span className="transport-plot-status transport-plot-status-warning">
              {data.warnings.length} warning{data.warnings.length === 1 ? "" : "s"}
            </span>
          )}
        </div>
      </div>

      <aside className="band-plot-sidebar transport-plot-sidebar">
        <div className="band-plot-controls">
          <button type="button" onClick={resetView} className="band-plot-reset">
            Reset View
          </button>
          <button type="button" onClick={resetAll} className="band-plot-export">
            Reset Viewer
          </button>
          <span className="band-plot-hint">Scroll: zoom Y | Shift+Scroll: pan range</span>
        </div>

        <div className="band-plot-control-panel">
          <div className="band-control-section">
            <div className="band-control-header transport-static-header">
              <span>Transport View</span>
              <InfoTooltip text={metricTooltip} />
            </div>
            <div className="band-control-grid transport-plot-meta-panel">
              <div className="transport-plot-metric-tabs" role="tablist" aria-label="Transport metric">
                {TRANSPORT_METRIC_KEYS.map((metricKey) => {
                  const definition = getTransportMetricDefinition(metricKey);
                  const isActive = metricKey === state.metric;
                  return (
                    <button
                      key={metricKey}
                      type="button"
                      className={`transport-plot-metric-button ${isActive ? "is-active" : ""}`}
                      onClick={() => setMetric(metricKey)}
                      role="tab"
                      aria-selected={isActive}
                      title={joinTransportTooltip(
                        definition.tooltip,
                        definition.tauBehavior === "dependent"
                          ? "Displayed values scale with the transport relaxation time."
                          : definition.tauBehavior === "independent"
                            ? "Displayed values are treated as relaxation-time independent in this workflow."
                            : undefined,
                      )}
                    >
                      <span className="transport-plot-metric-symbol">{definition.symbol}</span>
                      <span className="transport-plot-metric-label">{definition.shortLabel}</span>
                      <span className="transport-plot-metric-badge">{getTransportTauBadge(definition)}</span>
                    </button>
                  );
                })}
              </div>
            </div>
          </div>

          <div className="band-control-section">
            <button type="button" className="band-control-header" onClick={toggleHeatmapOpen}>
              <span className={`collapse-icon ${state.heatmapOpen ? "expanded" : ""}`}>▶</span>
              μ-T Map
            </button>
            {state.heatmapOpen && (
              <div className="band-control-grid transport-heatmap-panel">
                <div className="band-control-row">
                  <label>{state.axisMode === "mu" ? "Temperature Slice" : "Chemical-Potential Slice"}</label>
                  {state.axisMode === "mu" ? (
                    <select
                      value={safeTemperatureIndex}
                      onChange={(event) => setSelectedTemperatureIndex(Number(event.target.value))}
                    >
                      {data.temperature_values_k.map((temperature, index) => (
                        <option key={`${temperature}-${index}`} value={index}>
                          {formatTransportNumber(
                            temperature,
                            state.axisSettings.precisionMode,
                            state.axisSettings.precision,
                          )} K
                        </option>
                      ))}
                    </select>
                  ) : (
                    <select
                      value={safeMuIndex}
                      onChange={(event) => setSelectedMuIndex(Number(event.target.value))}
                    >
                      {data.mu_offsets_ev.map((muOffset, index) => (
                        <option key={`${muOffset}-${index}`} value={index}>
                          {formatTransportNumber(
                            state.axisSettings.showAbsoluteMu ? (data.mu_values_ev[index] ?? muOffset) : muOffset,
                            state.axisSettings.precisionMode,
                            state.axisSettings.precision,
                          )} eV
                        </option>
                      ))}
                    </select>
                  )}
                </div>

                <div className="band-control-note">
                  Click the heatmap to move through the chemical-potential and temperature grid while keeping the main plot focused on a single slice.
                </div>

                <TransportHeatmap
                  metricDefinition={metricDefinition}
                  axisMode={state.axisMode}
                  axisSettings={state.axisSettings}
                  grid={heatmapGrid}
                  selectedTemperatureIndex={safeTemperatureIndex}
                  selectedMuIndex={safeMuIndex}
                  onSelectCell={applyHeatmapSelection}
                />
              </div>
            )}
          </div>

          <div className="band-control-section">
            <button type="button" className="band-control-header" onClick={toggleSettingsOpen}>
              <span className={`collapse-icon ${state.settingsOpen ? "expanded" : ""}`}>▶</span>
              Plot Settings
            </button>
            {state.settingsOpen && (
              <div className="band-control-grid">
                <div className="band-control-row">
                  <label>μ Display</label>
                  <select
                    value={state.axisSettings.showAbsoluteMu ? "absolute" : "relative"}
                    onChange={(event) => updateAxisSettings({
                      showAbsoluteMu: event.target.value === "absolute",
                    })}
                  >
                    <option value="relative">Δμ relative to E_F</option>
                    <option value="absolute">Absolute μ</option>
                  </select>
                </div>

                <div className="band-control-row">
                  <label>Y Range</label>
                  <div className="band-control-range-inputs">
                    <input
                      type="number"
                      step="any"
                      placeholder="min"
                      value={state.axisSettings.manualYMinInput}
                      onChange={(event) => updateAxisSettings({ manualYMinInput: event.target.value })}
                      onKeyDown={manualRangeKeyDown}
                    />
                    <span className="band-control-range-separator">to</span>
                    <input
                      type="number"
                      step="any"
                      placeholder="max"
                      value={state.axisSettings.manualYMaxInput}
                      onChange={(event) => updateAxisSettings({ manualYMaxInput: event.target.value })}
                      onKeyDown={manualRangeKeyDown}
                    />
                  </div>
                  <button type="button" className="band-control-apply" onClick={() => applyManualYRange()}>
                    Apply
                  </button>
                  {rangeError && <span className="band-control-range-error">{rangeError}</span>}
                </div>

                <div className="band-control-row">
                  <label>
                    <input
                      type="checkbox"
                      checked={state.axisSettings.keepYRangeAcrossSlices}
                      onChange={(event) => updateAxisSettings({
                        keepYRangeAcrossSlices: event.target.checked,
                      })}
                    />
                    Keep Y range across slices
                    <InfoTooltip text="Freeze the current auto-scaled Y range while selecting different heatmap rows or columns. Changing the transport metric, component, or plot axis re-bases the auto range. Manual Y min and max still take precedence." />
                  </label>
                  <div className="band-control-note">
                    Useful for stepping through neighboring heatmap slices without the plot re-scaling each time.
                  </div>
                </div>

                <div className="band-control-row">
                  <label>Tick Density</label>
                  <select
                    value={state.axisSettings.tickDensity}
                    onChange={(event) => updateAxisSettings({
                      tickDensity: event.target.value as typeof state.axisSettings.tickDensity,
                    })}
                  >
                    <option value="sparse">Sparse</option>
                    <option value="normal">Normal</option>
                    <option value="dense">Dense</option>
                  </select>
                </div>

                <div className="band-control-row">
                  <label>Precision</label>
                  <div className="transport-inline-settings">
                    <select
                      value={state.axisSettings.precisionMode}
                      onChange={(event) => updateAxisSettings({
                        precisionMode: event.target.value as typeof state.axisSettings.precisionMode,
                      })}
                    >
                      <option value="auto">Auto</option>
                      <option value="manual">Manual</option>
                    </select>
                    {state.axisSettings.precisionMode === "manual" && (
                      <input
                        type="number"
                        min={0}
                        max={8}
                        step={1}
                        value={state.axisSettings.precision}
                        onChange={(event) => updateAxisSettings({
                          precision: Math.min(8, Math.max(0, Number(event.target.value))),
                        })}
                      />
                    )}
                  </div>
                </div>

                <div className="band-control-row">
                  <label>X Label Override</label>
                  <input
                    type="text"
                    value={state.axisSettings.xLabelOverride}
                    onChange={(event) => updateAxisSettings({ xLabelOverride: event.target.value })}
                    placeholder={xAxisLabel}
                  />
                </div>

                <div className="band-control-row">
                  <label>Y Label Override</label>
                  <input
                    type="text"
                    value={state.axisSettings.yLabelOverride}
                    onChange={(event) => updateAxisSettings({ yLabelOverride: event.target.value })}
                    placeholder={metricDefinition.symbol}
                  />
                </div>

                <div className="band-control-row">
                  <label>Unit Suffix Override</label>
                  <input
                    type="text"
                    value={state.axisSettings.unitSuffixOverride}
                    onChange={(event) => updateAxisSettings({ unitSuffixOverride: event.target.value })}
                    placeholder={metricDefinition.defaultUnit}
                  />
                </div>

                <div className="transport-checkbox-row">
                  <label>
                    <input
                      type="checkbox"
                      checked={state.axisSettings.showGrid}
                      onChange={(event) => updateAxisSettings({ showGrid: event.target.checked })}
                    />
                    Grid
                  </label>
                  <label>
                    <input
                      type="checkbox"
                      checked={state.axisSettings.showZeroLine}
                      onChange={(event) => updateAxisSettings({ showZeroLine: event.target.checked })}
                    />
                    Zero line
                  </label>
                </div>

                <button
                  type="button"
                  className="band-control-apply"
                  onClick={resetAxisSettings}
                >
                  Reset to Defaults
                </button>
              </div>
            )}
          </div>

          {data.tdf && (
            <div className="band-control-section">
              <button type="button" className="band-control-header" onClick={toggleTdfOpen}>
                <span className={`collapse-icon ${state.tdfOpen ? "expanded" : ""}`}>▶</span>
                <span>
                  TDF
                  <InfoTooltip text={tdfDefinition.tooltip} />
                </span>
              </button>
              {state.tdfOpen && (
                <div className="band-control-grid">
                  <div className="band-control-row">
                    <label>Component</label>
                    <select
                      value={state.tdfComponent}
                      onChange={(event) => setTdfComponent(event.target.value)}
                    >
                      {tdfOptions.map((option) => (
                        <option key={option.value} value={option.value}>{option.label}</option>
                      ))}
                    </select>
                  </div>
                  <TransportTdfPlot
                    data={data.tdf}
                    component={state.tdfComponent}
                    axisSettings={state.axisSettings}
                    metricDefinition={tdfDefinition}
                  />
                </div>
              )}
            </div>
          )}

          <div className="band-control-section">
            <button type="button" className="band-control-header" onClick={toggleRunContextOpen}>
              <span className={`collapse-icon ${state.runContextOpen ? "expanded" : ""}`}>▶</span>
              Run Context
            </button>
            {state.runContextOpen && (
              <div className="band-control-grid">
                <div className="transport-context-grid">
                  <div className="band-control-readout transport-context-card">
                    <span>Engine</span>
                    <strong>{data.engine || "boltzwann"}</strong>
                  </div>
                  <div className="band-control-readout transport-context-card">
                    <span>Seedname</span>
                    <strong>{data.seedname}</strong>
                  </div>
                  <div className="band-control-readout transport-context-card">
                    <span>μ points</span>
                    <strong>{data.mu_values_ev.length}</strong>
                  </div>
                  <div className="band-control-readout transport-context-card">
                    <span>T points</span>
                    <strong>{data.temperature_values_k.length}</strong>
                  </div>
                  <div className="band-control-readout transport-context-card">
                    <span>Artifacts</span>
                    <strong>{formatBytes(artifactBytes)}</strong>
                  </div>
                  <div className="band-control-readout transport-context-card">
                    <span>Tensor mode</span>
                    <strong>{data.is_2d ? "2D" : "3D"}</strong>
                  </div>
                </div>

                {data.warnings.length > 0 && (
                  <div className="band-control-warning transport-context-stack">
                    {data.warnings.map((warning, index) => (
                      <p key={`${warning}-${index}`}>{warning}</p>
                    ))}
                  </div>
                )}

                {data.notes.length > 0 && (
                  <div className="band-control-note transport-context-stack">
                    {data.notes.map((note, index) => (
                      <p key={`${note}-${index}`}>{note}</p>
                    ))}
                  </div>
                )}
              </div>
            )}
          </div>
        </div>
      </aside>
    </div>
  );
}
