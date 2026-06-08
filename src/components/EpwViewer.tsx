import { useMemo, useState } from "react";
import { AppHeaderPortal } from "./AppHeaderPortal";
import { InfoTooltip } from "./InfoTooltip";
import { EpwCompactHeatmap, EpwLinePlot } from "./epw/EpwLinePlot";
import {
  EpwMobilityDataset,
  EpwMobilityMetric,
  EpwParsedTable,
  EpwSeriesPoint,
  EpwViewerPayload,
  buildEpwMobilitySeries,
  buildEpwTableSeries,
  collectEpwWarnings,
  formatEpwBytes,
  formatEpwNumber,
  getDefaultEpwTab,
  getEpwComponentOptions,
  hasEpwSuperconductivityData,
  hasEpwTransportData,
  mergeEpwArtifacts,
} from "../lib/engines/qe/epw";

export type {
  EpwArtifactManifestEntry,
  EpwErrorRecord,
  EpwNamedArtifact,
  EpwResultSummary,
  EpwSourceRef,
  EpwViewerData,
  EpwViewerPayload,
} from "../lib/engines/qe/epw";

interface EpwViewerProps {
  payload: EpwViewerPayload;
  onBack: () => void;
}

interface MetricCardProps {
  label: string;
  value: unknown;
  unit?: string;
  tooltip?: string;
}

const MOBILITY_METRICS: Array<{ key: EpwMobilityMetric; label: string; unit: string; tooltip: string }> = [
  {
    key: "mobility",
    label: "Mobility",
    unit: "cm^2/Vs",
    tooltip: "Carrier mobility parsed from EPW transport output. SERTA values use self-energy relaxation time; IBTE values are iterative Boltzmann transport results when available.",
  },
  {
    key: "fermi",
    label: "Fermi",
    unit: "eV",
    tooltip: "Fine-grid Fermi level reported by EPW for the carrier density and temperature point.",
  },
  {
    key: "density",
    label: "Density",
    unit: "cm^-3",
    tooltip: "Carrier density used in the EPW mobility calculation.",
  },
  {
    key: "population",
    label: "Population",
    unit: "e/cell",
    tooltip: "Population term reported by EPW in tensor-style transport tables. This is mainly a diagnostic quantity.",
  },
];

function sourceIdLabel(source: { calc_id?: string; calc_type?: string } | null | undefined): string {
  const raw = String(source?.calc_id || "").trim();
  return raw ? raw.slice(0, 8) : "N/A";
}

function joinMesh(mesh: [number, number, number] | undefined): string {
  return Array.isArray(mesh) ? mesh.join(" x ") : "N/A";
}

function metricValue(data: EpwViewerPayload["data"], key: string): unknown {
  return data.result_summary?.core_metrics?.[key] ?? data.result_summary?.unknown_metrics?.[key] ?? null;
}

function MetricCard({ label, value, unit, tooltip }: MetricCardProps) {
  return (
    <div className="epw-metric-card">
      <span>
        {label}
        {tooltip && <InfoTooltip text={tooltip} />}
      </span>
      <strong>{formatEpwNumber(value, 4)}{unit ? ` ${unit}` : ""}</strong>
    </div>
  );
}

function EmptyState({ title, body }: { title: string; body: string }) {
  return (
    <div className="epw-empty-state">
      <h3>{title}</h3>
      <p>{body}</p>
    </div>
  );
}

function TablePreview({ table, maxRows = 12 }: { table: EpwParsedTable; maxRows?: number }) {
  if (table.skipped) {
    return (
      <div className="epw-table-skip">
        {table.skip_reason || "This file was not previewed."}
      </div>
    );
  }
  if (!table.rows?.length) {
    return <div className="epw-table-skip">No numeric rows were parsed.</div>;
  }
  return (
    <div className="epw-table-wrap">
      <table className="epw-table">
        <thead>
          <tr>
            {(table.column_labels || []).map((label, index) => (
              <th key={`${label}-${index}`}>{label || `c${index + 1}`}</th>
            ))}
          </tr>
        </thead>
        <tbody>
          {table.rows.slice(0, maxRows).map((row, rowIndex) => (
            <tr key={rowIndex}>
              {row.map((value, colIndex) => (
                <td key={`${rowIndex}-${colIndex}`}>{value == null ? "-" : formatEpwNumber(value, 5)}</td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
      {table.rows.length > maxRows && (
        <p className="epw-table-note">Showing {maxRows} of {table.rows.length} parsed rows.</p>
      )}
    </div>
  );
}

function OverviewTab({ payload }: { payload: EpwViewerPayload }) {
  const data = payload.data || {};
  const summary = data.result_summary || {};
  const warnings = collectEpwWarnings(data);
  const parsedTables = data.parsed_tables || [];
  const artifacts = mergeEpwArtifacts(data);
  const manifests = Array.isArray(data.sources?.manifests) ? data.sources!.manifests! : [];
  const parsedFamilies = Array.from(new Set(parsedTables.map((table) => table.family))).filter(Boolean);

  return (
    <div className="epw-tab-panel">
      <section className="epw-section">
        <div className="epw-section-header">
          <h3>Run Summary</h3>
          <span className={`epw-status-pill ${summary.completed ? "is-complete" : "is-partial"}`}>
            {summary.completed ? "Complete" : "Incomplete"}
          </span>
        </div>
        <div className="epw-metric-grid">
          <MetricCard label="Schema" value={data.schema_version ?? "v1"} />
          <MetricCard label="Elapsed" value={summary.elapsed_seconds} unit="s" />
          <MetricCard label="Parsed Results" value={summary.parse_partial ? "Limited" : "Complete"} />
          <MetricCard label="Artifacts" value={artifacts.length} />
          <MetricCard label="Parsed Tables" value={parsedTables.length} />
          <MetricCard label="Warnings" value={warnings.length} />
        </div>
      </section>

      <section className="epw-section">
        <h3>Top Metrics</h3>
        <div className="epw-metric-grid">
          <MetricCard label="lambda" value={metricValue(data, "lambda")} tooltip="Dimensionless electron-phonon coupling strength. Use mesh and smearing convergence checks before treating this as final." />
          <MetricCard label="lambda_tr" value={metricValue(data, "lambda_tr")} tooltip="Transport electron-phonon coupling, when EPW reports it separately from lambda." />
          <MetricCard label="Tc" value={metricValue(data, "tc")} unit="K" tooltip="Critical temperature estimate parsed from EPW output. The exact model is shown in the superconductivity tab when available." />
          <MetricCard label="Mobility Sets" value={metricValue(data, "epw_mobility_datasets")} />
          <MetricCard label="Self-Energy Modes" value={metricValue(data, "self_energy_modes")} />
          <MetricCard label="Coverage" value={(summary.parse_coverage || []).join(", ") || "N/A"} />
        </div>
      </section>

      <section className="epw-section">
        <h3>Provenance</h3>
        <div className="epw-context-grid">
          <MetricCard label="Source Phonon" value={sourceIdLabel(data.sources?.phonon)} />
          <MetricCard label="Source Wannier" value={sourceIdLabel(data.sources?.wannier)} />
          <MetricCard label="Source SCF" value={sourceIdLabel(data.sources?.scf)} />
          <MetricCard label="Manifest Entries" value={manifests.length} />
          <MetricCard label="K Mesh" value={joinMesh(data.input?.k_mesh)} />
          <MetricCard label="Q Mesh" value={joinMesh(data.input?.q_mesh)} />
          <MetricCard label="fsthick" value={data.input?.fsthick_ev} unit="eV" />
          <MetricCard label="degaussw" value={data.input?.degaussw_ev} unit="eV" />
          <MetricCard label="nbndsub" value={data.input?.nbndsub ?? "N/A"} />
          <MetricCard label="MPI Pools" value={data.runtime?.pools ?? "N/A"} />
          <MetricCard label="Sync Mode" value={data.runtime?.artifact_sync_mode || "N/A"} />
          <MetricCard label="Table Families" value={parsedFamilies.join(", ") || "N/A"} />
        </div>
      </section>

      {warnings.length > 0 && (
        <section className="epw-section">
          <h3>Warnings</h3>
          <ul className="epw-warning-list">
            {warnings.map((warning, index) => <li key={`${warning}-${index}`}>{warning}</li>)}
          </ul>
        </section>
      )}
    </div>
  );
}

function getDatasetLabel(dataset: EpwMobilityDataset, index: number): string {
  const pieces = [
    dataset.method || "EPW",
    dataset.carrier_type || "carrier",
    dataset.iteration != null ? `iter ${dataset.iteration}` : null,
  ].filter(Boolean);
  return pieces.join(" / ") || `dataset ${index + 1}`;
}

function yLabelForMobility(metric: EpwMobilityMetric, component: string): string {
  if (metric === "mobility") return `Mobility ${component} (cm^2/Vs)`;
  if (metric === "fermi") return "Fermi level (eV)";
  if (metric === "density") return "Carrier density (cm^-3)";
  return "Population (e/cell)";
}

function TransportTab({ payload }: { payload: EpwViewerPayload }) {
  const datasets = payload.data.transport?.mobility || [];
  const [datasetIndex, setDatasetIndex] = useState(0);
  const [metric, setMetric] = useState<EpwMobilityMetric>("mobility");
  const selectedDataset = datasets[Math.min(datasetIndex, Math.max(datasets.length - 1, 0))];
  const componentOptions = getEpwComponentOptions(selectedDataset?.component_labels || []);
  const [component, setComponent] = useState(componentOptions[0]?.value || "avg");
  const safeComponent = componentOptions.some((option) => option.value === component)
    ? component
    : componentOptions[0]?.value || "avg";
  const series = selectedDataset ? buildEpwMobilitySeries(selectedDataset, metric, safeComponent) : [];
  const activeMetric = MOBILITY_METRICS.find((entry) => entry.key === metric) || MOBILITY_METRICS[0];

  if (!datasets.length) {
    return (
      <EmptyState
        title="No EPW transport tables parsed"
        body="The run did not expose EPW mobility/scattering stdout sections in the supported SERTA or IBTE formats. Check that transport keywords such as scattering, iterative_bte, carrier, ltrans_crta, or related EPW options were enabled, then inspect Files + Log for raw output."
      />
    );
  }

  return (
    <div className="epw-tab-panel epw-transport-panel">
      <section className="epw-section">
        <div className="epw-section-header epw-controls-header">
          <h3>EPW Transport</h3>
          <div className="epw-control-row">
            <label>
              Dataset
              <InfoTooltip text="EPW may print multiple transport datasets: SERTA, each IBTE iteration, electron/hole carriers, or intermediate/final mobility blocks." />
            </label>
            <select value={datasetIndex} onChange={(event) => setDatasetIndex(Number(event.target.value))}>
              {datasets.map((dataset, index) => (
                <option key={`${dataset.method}-${dataset.carrier_type}-${index}`} value={index}>
                  {getDatasetLabel(dataset, index)}
                </option>
              ))}
            </select>
          </div>
          {metric === "mobility" && (
            <div className="epw-control-row">
              <label>
                Component
                <InfoTooltip text={componentOptions.find((option) => option.value === safeComponent)?.tooltip || "Tensor component parsed from EPW output."} />
              </label>
              <select value={safeComponent} onChange={(event) => setComponent(event.target.value)}>
                {componentOptions.map((option) => <option key={option.value} value={option.value}>{option.label}</option>)}
              </select>
            </div>
          )}
        </div>

        <div className="epw-metric-tabs" role="tablist" aria-label="EPW transport metric">
          {MOBILITY_METRICS.map((entry) => (
            <button
              key={entry.key}
              type="button"
              className={metric === entry.key ? "is-active" : ""}
              onClick={() => setMetric(entry.key)}
              title={entry.tooltip}
            >
              <span>{entry.label}</span>
              <small>{entry.unit}</small>
            </button>
          ))}
        </div>

        <EpwLinePlot
          title={`${activeMetric.label} / ${getDatasetLabel(selectedDataset, datasetIndex)}`}
          xLabel="Temperature (K)"
          yLabel={yLabelForMobility(metric, safeComponent)}
          series={series}
          color="#0f766e"
        />

        {metric === "mobility" && selectedDataset.component_labels.length > 1 && (
          <EpwCompactHeatmap
            title="Component by Temperature"
            columnLabels={(selectedDataset.temperature_values_k || []).map((temp) => `${formatEpwNumber(temp, 0)} K`)}
            rowLabels={selectedDataset.component_labels || []}
            values={selectedDataset.mobility_values || []}
            formatter={(value) => formatEpwNumber(value, 2)}
          />
        )}
      </section>

      <section className="epw-section">
        <h3>Transport Table</h3>
        <div className="epw-table-wrap">
          <table className="epw-table">
            <thead>
              <tr>
                <th>T (K)</th>
                <th>Fermi (eV)</th>
                <th>Density (cm^-3)</th>
                <th>Population</th>
                {selectedDataset.component_labels.map((label) => <th key={label}>{label}</th>)}
              </tr>
            </thead>
            <tbody>
              {(selectedDataset.temperature_values_k || []).map((temperature, rowIndex) => (
                <tr key={`${temperature}-${rowIndex}`}>
                  <td>{formatEpwNumber(temperature, 2)}</td>
                  <td>{formatEpwNumber(selectedDataset.fermi_values_ev?.[rowIndex], 4)}</td>
                  <td>{formatEpwNumber(selectedDataset.density_values_cm3?.[rowIndex], 4)}</td>
                  <td>{formatEpwNumber(selectedDataset.population_values?.[rowIndex], 4)}</td>
                  {selectedDataset.component_labels.map((_, componentIndex) => (
                    <td key={`${rowIndex}-${componentIndex}`}>
                      {formatEpwNumber(selectedDataset.mobility_values?.[componentIndex]?.[rowIndex], 4)}
                    </td>
                  ))}
                </tr>
              ))}
            </tbody>
          </table>
        </div>
        <div className="epw-note-row">
          <span>Method: {selectedDataset.method || "N/A"}</span>
          <span>Carrier: {selectedDataset.carrier_type || "N/A"}</span>
          <span>Converged: {selectedDataset.converged == null ? "N/A" : selectedDataset.converged ? "Yes" : "No"}</span>
          <span>Max error: {formatEpwNumber(selectedDataset.max_error, 4)}</span>
        </div>
      </section>
    </div>
  );
}

function SuperconductivityTab({ payload }: { payload: EpwViewerPayload }) {
  const data = payload.data.superconductivity;
  const [iterationMetric, setIterationMetric] = useState<"delta" | "ethr" | "znorm">("delta");

  if (!data || !hasEpwSuperconductivityData(payload.data)) {
    return (
      <EmptyState
        title="No structured superconductivity data"
        body="The run did not expose supported Eliashberg, a2F, lambda, Tc, or gap output. Use Files + Log to inspect raw EPW output and generated artifacts."
      />
    );
  }

  const iterations = data.eliashberg_iterations || [];
  const iterationSeries: EpwSeriesPoint[] = iterations.map((row, index) => ({
    index,
    x: row.iteration,
    y: iterationMetric === "delta" ? row.delta_mev : iterationMetric === "ethr" ? row.ethr : row.znorm,
    label: row.temperature_k == null ? undefined : `${formatEpwNumber(row.temperature_k, 2)} K`,
  }));
  const iterationLabel = iterationMetric === "delta" ? "Delta (meV)" : iterationMetric === "ethr" ? "ethr" : "znorm";
  const spectralTables = data.spectral_tables || [];

  return (
    <div className="epw-tab-panel">
      <section className="epw-section">
        <div className="epw-section-header">
          <h3>Superconductivity Summary</h3>
          <span className={`epw-status-pill ${data.eliashberg_converged ? "is-complete" : "is-partial"}`}>
            {data.eliashberg_converged ? "Eliashberg converged" : "Convergence N/A"}
          </span>
        </div>
        <div className="epw-metric-grid">
          <MetricCard label="lambda" value={data.lambda ?? data.electron_phonon_coupling} tooltip="Dimensionless electron-phonon coupling strength. Converge this with respect to k/q meshes, Fermi-surface window, and smearing." />
          <MetricCard label="lambda_tr" value={data.lambda_tr} tooltip="Transport electron-phonon coupling, when EPW prints it separately." />
          <MetricCard label="Tc McMillan" value={data.tc_mcmillan_k} unit="K" tooltip="McMillan critical-temperature estimate using the parsed Coulomb pseudopotential mu* when present." />
          <MetricCard label="Tc Allen-Dynes" value={data.tc_allen_dynes_k} unit="K" tooltip="Allen-Dynes modified McMillan critical-temperature estimate." />
          <MetricCard label="Tc SISSO" value={data.tc_sisso_k} unit="K" tooltip="Machine-learning Tc estimate printed by EPW when enabled." />
          <MetricCard label="w_log" value={data.w_log_mev} unit="meV" />
          <MetricCard label="BCS Gap" value={data.bcs_gap_mev} unit="meV" />
          <MetricCard label="mu*" value={data.muc} />
          <MetricCard label="Frequency Points" value={data.frequency_points} />
          <MetricCard label="wscut" value={data.frequency_cutoff_ev} unit="eV" />
        </div>
      </section>

      {iterations.length > 0 && (
        <section className="epw-section">
          <div className="epw-section-header epw-controls-header">
            <h3>Eliashberg Convergence</h3>
            <div className="epw-control-row">
              <label>Metric</label>
              <select value={iterationMetric} onChange={(event) => setIterationMetric(event.target.value as typeof iterationMetric)}>
                <option value="delta">delta</option>
                <option value="ethr">ethr</option>
                <option value="znorm">znorm</option>
              </select>
            </div>
          </div>
          <EpwLinePlot
            title="Eliashberg Iterations"
            xLabel="Iteration"
            yLabel={iterationLabel}
            series={iterationSeries}
            color="#be185d"
          />
        </section>
      )}

      {(data.gap_summaries || []).length > 0 && (
        <section className="epw-section">
          <h3>Gap Summary</h3>
          <div className="epw-table-wrap">
            <table className="epw-table">
              <thead>
                <tr>
                  <th>T (K)</th>
                  <th>Free Energy (meV)</th>
                  <th>Gap Min (meV)</th>
                  <th>Gap Max (meV)</th>
                </tr>
              </thead>
              <tbody>
                {(data.gap_summaries || []).map((row, index) => (
                  <tr key={index}>
                    <td>{formatEpwNumber(row.temperature_k, 3)}</td>
                    <td>{formatEpwNumber(row.free_energy_mev, 5)}</td>
                    <td>{formatEpwNumber(row.gap_min_mev, 5)}</td>
                    <td>{formatEpwNumber(row.gap_max_mev, 5)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </section>
      )}

      {spectralTables.length > 0 && (
        <section className="epw-section">
          <h3>a2F / DOS Tables</h3>
          <div className="epw-table-card-grid">
            {spectralTables.slice(0, 3).map((table) => (
              <div className="epw-table-card" key={table.file_name}>
                <EpwLinePlot
                  title={table.title || table.file_name}
                  xLabel={table.column_labels?.[0] || "x"}
                  yLabel={table.column_labels?.[1] || "value"}
                  series={buildEpwTableSeries(table, 0, 1)}
                  color="#b45309"
                />
                <TablePreview table={table} maxRows={6} />
              </div>
            ))}
          </div>
        </section>
      )}
    </div>
  );
}

function SpectralTab({ payload }: { payload: EpwViewerPayload }) {
  const tables = useMemo(() => {
    const fromSpectral = payload.data.spectral?.tables || [];
    const allParsed = payload.data.parsed_tables || [];
    const merged = new Map<string, EpwParsedTable>();
    for (const table of [...fromSpectral, ...allParsed]) {
      merged.set(table.file_name, table);
    }
    return Array.from(merged.values());
  }, [payload.data.parsed_tables, payload.data.spectral?.tables]);
  const modes = payload.data.spectral?.self_energy_modes || [];
  const [tableIndex, setTableIndex] = useState(0);
  const [xColumn, setXColumn] = useState(0);
  const [yColumn, setYColumn] = useState(1);
  const selectedTable = tables[Math.min(tableIndex, Math.max(tables.length - 1, 0))];
  const maxColumns = selectedTable?.column_labels?.length || selectedTable?.rows?.[0]?.length || 0;
  const safeXColumn = Math.min(xColumn, Math.max(maxColumns - 1, 0));
  const safeYColumn = Math.min(yColumn, Math.max(maxColumns - 1, 0));

  return (
    <div className="epw-tab-panel">
      {modes.length > 0 && (
        <section className="epw-section">
          <h3>Lambda / Self-Energy Modes</h3>
          <div className="epw-table-wrap">
            <table className="epw-table">
              <thead>
                <tr>
                  <th>Mode</th>
                  <th>lambda</th>
                  <th>lambda_tr</th>
                  <th>gamma (meV)</th>
                  <th>gamma_tr (meV)</th>
                  <th>omega (meV)</th>
                </tr>
              </thead>
              <tbody>
                {modes.slice(0, 80).map((mode, index) => (
                  <tr key={`${mode.mode_label}-${index}`}>
                    <td>{mode.mode_label}</td>
                    <td>{formatEpwNumber(mode.lambda, 5)}</td>
                    <td>{formatEpwNumber(mode.lambda_tr, 5)}</td>
                    <td>{formatEpwNumber(mode.gamma_mev, 5)}</td>
                    <td>{formatEpwNumber(mode.gamma_tr_mev, 5)}</td>
                    <td>{formatEpwNumber(mode.omega_mev, 5)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </section>
      )}

      {tables.length === 0 ? (
        <EmptyState
          title="No spectral or generic numeric tables"
          body="No chartable EPW artifacts were parsed. The raw files are still listed under Files + Log."
        />
      ) : (
        <section className="epw-section">
          <div className="epw-section-header epw-controls-header">
            <h3>Parsed Numeric Artifacts</h3>
            <div className="epw-control-row">
              <label>Table</label>
              <select value={tableIndex} onChange={(event) => {
                setTableIndex(Number(event.target.value));
                setXColumn(0);
                setYColumn(1);
              }}>
                {tables.map((table, index) => (
                  <option key={table.file_name} value={index}>
                    {table.family} / {table.file_name}
                  </option>
                ))}
              </select>
            </div>
            {maxColumns > 1 && (
              <>
                <div className="epw-control-row">
                  <label>X</label>
                  <select value={safeXColumn} onChange={(event) => setXColumn(Number(event.target.value))}>
                    {selectedTable.column_labels.map((label, index) => <option key={`${label}-${index}`} value={index}>{label}</option>)}
                  </select>
                </div>
                <div className="epw-control-row">
                  <label>Y</label>
                  <select value={safeYColumn} onChange={(event) => setYColumn(Number(event.target.value))}>
                    {selectedTable.column_labels.map((label, index) => <option key={`${label}-${index}`} value={index}>{label}</option>)}
                  </select>
                </div>
              </>
            )}
          </div>

          {selectedTable && !selectedTable.skipped && maxColumns > 1 && (
            <EpwLinePlot
              title={selectedTable.title || selectedTable.file_name}
              xLabel={selectedTable.column_labels?.[safeXColumn] || `c${safeXColumn + 1}`}
              yLabel={selectedTable.column_labels?.[safeYColumn] || `c${safeYColumn + 1}`}
              series={buildEpwTableSeries(selectedTable, safeXColumn, safeYColumn)}
              color="#7c3aed"
            />
          )}
          {selectedTable && selectedTable.rows?.length > 0 && selectedTable.rows.length <= 60 && maxColumns > 2 && (
            <EpwCompactHeatmap
              title="Table Heatmap"
              columnLabels={selectedTable.column_labels.slice(1)}
              rowLabels={selectedTable.rows.map((row, index) => formatEpwNumber(row[0] ?? index, 3))}
              values={selectedTable.rows.map((row) => row.slice(1))}
            />
          )}
          {selectedTable && <TablePreview table={selectedTable} />}
        </section>
      )}
    </div>
  );
}

function FilesTab({ payload }: { payload: EpwViewerPayload }) {
  const data = payload.data || {};
  const rawOutput = typeof payload.rawOutput === "string" ? payload.rawOutput : "";
  const [query, setQuery] = useState("");
  const artifacts = mergeEpwArtifacts(data);
  const parsedTables = data.parsed_tables || [];
  const notes = [
    ...(data.result_summary?.notes || []),
    ...(data.transport?.notes || []),
    ...(data.superconductivity?.notes || []),
    ...(data.spectral?.notes || []),
  ];
  const filteredLog = query.trim()
    ? rawOutput
        .split("\n")
        .filter((line) => line.toLowerCase().includes(query.trim().toLowerCase()))
        .join("\n")
    : rawOutput;

  return (
    <div className="epw-tab-panel">
      <section className="epw-section">
        <h3>Output Files</h3>
        {artifacts.length > 0 ? (
          <div className="epw-artifact-list">
            {artifacts.map((artifact) => (
              <div key={artifact.file_name} className="epw-artifact-row">
                <span>{artifact.file_name}</span>
                <strong>{formatEpwBytes(Number(artifact.size_bytes) || 0)}</strong>
              </div>
            ))}
          </div>
        ) : (
          <p className="param-hint">No output artifacts recorded.</p>
        )}
      </section>

      <section className="epw-section">
        <h3>Parsed Table Inventory</h3>
        {parsedTables.length > 0 ? (
          <div className="epw-artifact-list">
            {parsedTables.map((table) => (
              <div key={table.file_name} className="epw-artifact-row">
                <span>{table.family} / {table.file_name}</span>
                <strong>{table.skipped ? "Skipped" : `${table.rows?.length || 0} rows`}</strong>
              </div>
            ))}
          </div>
        ) : (
          <p className="param-hint">No chartable numeric tables were parsed.</p>
        )}
      </section>

      {notes.length > 0 && (
        <section className="epw-section">
          <h3>Parser Notes</h3>
          <ul className="epw-warning-list">
            {Array.from(new Set(notes)).map((note, index) => <li key={`${note}-${index}`}>{note}</li>)}
          </ul>
        </section>
      )}

      <section className="epw-section epw-log-section">
        <div className="epw-section-header epw-controls-header">
          <h3>Run Log</h3>
          <div className="epw-control-row">
            <label>Filter</label>
            <input value={query} onChange={(event) => setQuery(event.target.value)} placeholder="warning, lambda, mobility..." />
          </div>
        </div>
        <pre className="output-text epw-log-output">{filteredLog || "No saved log output available for this calculation."}</pre>
      </section>
    </div>
  );
}

export function EpwViewer({ payload, onBack }: EpwViewerProps) {
  const data = payload.data || {};
  const [activeTab, setActiveTab] = useState(() => getDefaultEpwTab(data));
  const warnings = collectEpwWarnings(data);
  const tabs = [
    { key: "overview", label: "Overview", badge: data.schema_version ? `v${data.schema_version}` : "v1" },
    { key: "transport", label: "Transport", badge: hasEpwTransportData(data) ? String(data.transport?.mobility?.length || 0) : "" },
    { key: "superconductivity", label: "Superconductivity", badge: hasEpwSuperconductivityData(data) ? "data" : "" },
    { key: "spectral", label: "Spectral + Self-Energy", badge: data.parsed_tables?.length ? String(data.parsed_tables.length) : "" },
    { key: "files", label: "Files + Log", badge: warnings.length ? String(warnings.length) : "" },
  ] as const;

  return (
    <div className="bands-viewer-container epw-viewer-container">
      <AppHeaderPortal className="bands-viewer-header epw-viewer-header">
        <button className="back-button" onClick={onBack}>
          ← Back to Dashboard
        </button>
        <div className="epw-viewer-title">
          <span>EPW Results</span>
          <h2>{hasEpwTransportData(data) ? "Transport" : hasEpwSuperconductivityData(data) ? "Superconductivity" : "Results"}</h2>
        </div>
      </AppHeaderPortal>

      <div className="epw-tab-bar" role="tablist" aria-label="EPW result views">
        {tabs.map((tab) => (
          <button
            key={tab.key}
            type="button"
            className={activeTab === tab.key ? "is-active" : ""}
            onClick={() => setActiveTab(tab.key)}
            role="tab"
            aria-selected={activeTab === tab.key}
          >
            <span>{tab.label}</span>
            {tab.badge && <small>{tab.badge}</small>}
          </button>
        ))}
      </div>

      <div className="bands-viewer-content bands-viewer-content-stacked epw-viewer-content">
        {activeTab === "overview" && <OverviewTab payload={payload} />}
        {activeTab === "transport" && <TransportTab payload={payload} />}
        {activeTab === "superconductivity" && <SuperconductivityTab payload={payload} />}
        {activeTab === "spectral" && <SpectralTab payload={payload} />}
        {activeTab === "files" && <FilesTab payload={payload} />}
      </div>
    </div>
  );
}
