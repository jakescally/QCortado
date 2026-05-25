import { useEffect, useMemo, useRef, useState } from "react";
import { listen, type UnlistenFn } from "@tauri-apps/api/event";
import type { CrystalData } from "../../lib/types";
import {
  DEFAULT_WIEN2K_STRUCTURE_CONTROLS,
  discardWien2kStructureSession,
  normalizeWien2kCaseName,
  prepareWien2kStructureDraft,
  runWien2kStructureStage,
  saveWien2kStructureSource,
  startWien2kStructureSession,
  validateWien2kStructureControls,
} from "../../lib/engines/wien2k";
import type {
  Wien2kStructureControls,
  Wien2kStructureDraft,
  Wien2kStructureSession,
  Wien2kStructureStage,
  Wien2kStructureStageResult,
  Wien2kStructureSite,
} from "../../lib/engines/wien2k";
import { LiveOutputPanel } from "../LiveOutputPanel";
import { UnitCellViewer } from "../UnitCellViewer";

interface Wien2kStructureWizardProps {
  projectId: string;
  cifId: string;
  crystalData: CrystalData;
  onBack: () => void;
  onSaved: () => void;
}

function defaultCaseName(crystalData: CrystalData): string {
  const formula = crystalData.formula_sum ?? crystalData.formula_structural ?? "case";
  const sanitized = formula.replace(/[^A-Za-z0-9_.-]/g, "");
  return normalizeWien2kCaseName(sanitized) ?? "case";
}

function displaySiteValue(site: Wien2kStructureSite, controls: Wien2kStructureControls, field: "npt" | "r0" | "rmt") {
  const override = controls.siteOverrides.find((entry) => entry.siteIndex === site.siteIndex);
  return override?.[field] ?? site[field];
}

export function Wien2kStructureWizard({
  projectId,
  cifId,
  crystalData,
  onBack,
  onSaved,
}: Wien2kStructureWizardProps) {
  const [caseName, setCaseName] = useState(() => defaultCaseName(crystalData));
  const [controls, setControls] = useState<Wien2kStructureControls>(DEFAULT_WIEN2K_STRUCTURE_CONTROLS);
  const [draft, setDraft] = useState<Wien2kStructureDraft | null>(null);
  const [session, setSession] = useState<Wien2kStructureSession | null>(null);
  const [lastResult, setLastResult] = useState<Wien2kStructureStageResult | null>(null);
  const [outputLines, setOutputLines] = useState<string[]>([]);
  const [draftStale, setDraftStale] = useState(false);
  const [isPreparing, setIsPreparing] = useState(false);
  const [isRunning, setIsRunning] = useState(false);
  const [isSaving, setIsSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const outputUnlistenRef = useRef<UnlistenFn | null>(null);
  const draftRequestRef = useRef(0);
  const caseNameValid = normalizeWien2kCaseName(caseName) !== null;
  const controlsError = validateWien2kStructureControls(controls);
  const draftSiteSettingsKey = JSON.stringify(controls.siteOverrides.map(({ siteIndex, npt, r0 }) => ({
    siteIndex,
    npt,
    r0,
  })));
  const sites = lastResult?.sites ?? draft?.sites ?? [];
  const canOperate = Boolean(draft && !draftStale && !controlsError && caseNameValid);
  const displayedCrystalData = useMemo<CrystalData>(() => {
    if (!draft) return crystalData;
    let atomIndex = 0;
    return {
      ...crystalData,
      cell_length_a: { value: draft.cellParameters[0] },
      cell_length_b: { value: draft.cellParameters[1] },
      cell_length_c: { value: draft.cellParameters[2] },
      cell_angle_alpha: { value: draft.cellParameters[3] },
      cell_angle_beta: { value: draft.cellParameters[4] },
      cell_angle_gamma: { value: draft.cellParameters[5] },
      space_group_HM: draft.internationalSymbol,
      space_group_IT_number: draft.spacegroupNumber,
      atom_sites: draft.sites.flatMap((site) => site.positions.map((position) => {
        atomIndex += 1;
        return {
          label: `${site.symbol}${atomIndex}`,
          type_symbol: site.symbol,
          fract_x: position[0],
          fract_y: position[1],
          fract_z: position[2],
          occupancy: 1,
        };
      })),
      symmetry_operations: [],
      anisotropic_params: [],
    };
  }, [crystalData, draft]);

  useEffect(() => {
    return () => {
      outputUnlistenRef.current?.();
    };
  }, []);

  useEffect(() => {
    if (session || !caseNameValid || controlsError) return;

    const normalized = normalizeWien2kCaseName(caseName);
    if (!normalized) return;

    const requestId = draftRequestRef.current + 1;
    draftRequestRef.current = requestId;
    setDraftStale(true);
    setIsPreparing(true);
    setError(null);

    const timer = window.setTimeout(() => {
      void prepareWien2kStructureDraft(projectId, cifId, normalized, controls)
        .then((nextDraft) => {
          if (draftRequestRef.current !== requestId) return;
          setDraft(nextDraft);
          setLastResult(null);
          setDraftStale(false);
          setOutputLines(["Generated local WIEN2k structure draft from the selected project CIF."]);
        })
        .catch((reason) => {
          if (draftRequestRef.current !== requestId) return;
          setError(String(reason));
        })
        .finally(() => {
          if (draftRequestRef.current === requestId) {
            setIsPreparing(false);
          }
        });
    }, draft ? 150 : 0);

    return () => {
      window.clearTimeout(timer);
      if (draftRequestRef.current === requestId) {
        draftRequestRef.current += 1;
        setIsPreparing(false);
      }
    };
    // Stage-only controls (RMT reduction, RMT overrides and SGROUP tolerance)
    // do not alter the local draft file and are applied during native refinement.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [caseName, caseNameValid, cifId, controlsError, draftSiteSettingsKey, projectId, session]);

  async function attachOutputListener(sessionId: string) {
    outputUnlistenRef.current?.();
    outputUnlistenRef.current = await listen<string>(`wien2k-structure-output:${sessionId}`, (event) => {
      setOutputLines((current) => [...current, event.payload].slice(-1000));
    });
  }

  function updateSiteOverride(siteIndex: number, field: "npt" | "r0" | "rmt", rawValue: string) {
    const value = rawValue.trim() === "" ? null : Number(rawValue);
    setControls((current) => {
      const existing = current.siteOverrides.find((entry) => entry.siteIndex === siteIndex) ?? { siteIndex };
      const replacement = { ...existing, [field]: value };
      return {
        ...current,
        siteOverrides: [
          ...current.siteOverrides.filter((entry) => entry.siteIndex !== siteIndex),
          replacement,
        ],
      };
    });
    if (field !== "rmt") {
      setDraftStale(true);
    }
  }

  async function ensureSession(): Promise<Wien2kStructureSession | null> {
    if (session) return session;
    if (!draft) return null;
    const created = await startWien2kStructureSession(draft);
    setSession(created);
    setOutputLines((current) => [...current, `Remote refinement staged in ${created.remoteCaseDir}.`]);
    await attachOutputListener(created.sessionId);
    return created;
  }

  async function runStage(stage: Wien2kStructureStage) {
    if (!canOperate) return;
    setIsRunning(true);
    setError(null);
    try {
      const activeSession = await ensureSession();
      if (!activeSession) {
        throw new Error("A local structure draft is required before remote refinement.");
      }
      const result = await runWien2kStructureStage(activeSession.sessionId, stage, controls);
      setSession((current) => current ? { ...current, phase: result.phase, currentStruct: result.candidateStruct } : current);
      setLastResult(result);
      if (result.diagnostics.length > 0) {
        setOutputLines((current) => [...current, ...result.diagnostics.map((item) => `[diagnostic] ${item}`)]);
      }
    } catch (reason) {
      setError(String(reason));
    } finally {
      setIsRunning(false);
    }
  }

  async function discardSessionAndReturn() {
    setError(null);
    try {
      if (session) {
        await discardWien2kStructureSession(session.sessionId);
      }
      outputUnlistenRef.current?.();
      outputUnlistenRef.current = null;
      onBack();
    } catch (reason) {
      setError(String(reason));
    }
  }

  async function saveSource() {
    if (!session || !lastResult?.saveAllowed) return;
    setIsSaving(true);
    setError(null);
    try {
      await saveWien2kStructureSource(session.sessionId);
      outputUnlistenRef.current?.();
      outputUnlistenRef.current = null;
      onSaved();
    } catch (reason) {
      setError(String(reason));
    } finally {
      setIsSaving(false);
    }
  }

  const stageLabel = useMemo(() => {
    switch (session?.phase) {
      case "rmt_ready": return "RMT candidate ready";
      case "sgroup_ready": return "Space-group candidate ready";
      case "symmetry_ready": return "Symmetry candidate ready";
      default: return isPreparing ? "Generating draft" : draft ? "Draft ready" : "Preparing draft";
    }
  }, [draft, isPreparing, session?.phase]);

  return (
    <div className="wizard-container wien2k-structure-wizard">
      <div className="wizard-header">
        <button className="back-btn" type="button" disabled={isRunning || isSaving} onClick={() => void discardSessionAndReturn()}>
          Back
        </button>
        <h2>WIEN2k Structure</h2>
        <span className="wien2k-stage-badge">{stageLabel}</span>
      </div>
      {error && <div className="error-banner">{error}</div>}
      <div className="wien2k-structure-content">
        <section className="wien2k-structure-controls">
          <div className="wien2k-control-grid">
            <label>
              Case name
              <input
                value={caseName}
                disabled={Boolean(session)}
                onChange={(event) => {
                  setCaseName(event.target.value);
                  setDraftStale(true);
                }}
              />
            </label>
            <label>
              RMT reduction (%)
              <input
                type="number"
                min="0"
                max="99.999"
                step="0.1"
                value={controls.rmtReductionPercent}
                onChange={(event) => setControls((current) => ({
                  ...current,
                  rmtReductionPercent: Number(event.target.value),
                }))}
              />
            </label>
          </div>
          {!caseNameValid && (
            <p className="wien2k-validation">Case name may contain only letters, numbers, underscore, hyphen, and period.</p>
          )}
          {controlsError && <p className="wien2k-validation">{controlsError}</p>}
          {!draft && isPreparing && <p className="wien2k-preparing">Generating initial structure draft...</p>}

          {draft && (
            <>
              <div className="wien2k-summary">
                <h3>Standardized conventional cell</h3>
                <p>
                  {draft.internationalSymbol} (#{draft.spacegroupNumber}), lattice {draft.latticeType};
                  {" "}a={draft.cellParameters[0].toFixed(4)} A, b={draft.cellParameters[1].toFixed(4)} A,
                  {" "}c={draft.cellParameters[2].toFixed(4)} A
                </p>
              </div>
              <h3 className="wien2k-site-heading">Muffin-tin radii</h3>
              <table className="wien2k-site-table">
                <thead>
                  <tr><th>Site</th><th>Multiplicity</th><th>RMT</th></tr>
                </thead>
                <tbody>
                  {sites.map((site) => (
                    <tr key={site.siteIndex}>
                      <td>{site.symbol}</td>
                      <td>{site.positions.length}</td>
                      <td>
                        <input
                          type="number"
                          min="0.000001"
                          step="0.01"
                          disabled={Boolean(session && session.phase !== "rmt_ready")}
                          value={displaySiteValue(site, controls, "rmt")}
                          onChange={(event) => updateSiteOverride(site.siteIndex, "rmt", event.target.value)}
                        />
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
              <details className="wien2k-advanced">
                <summary>Advanced structure settings</summary>
                <label>
                  SGROUP tolerance
                  <input
                    type="number"
                    min="0.0000001"
                    max="0.001"
                    step="0.000001"
                    value={controls.sgroupTolerance}
                    onChange={(event) => setControls((current) => ({
                      ...current,
                      sgroupTolerance: Number(event.target.value),
                    }))}
                  />
                </label>
                <table className="wien2k-site-table">
                  <thead>
                    <tr><th>Site</th><th>NPT</th><th>R0</th></tr>
                  </thead>
                  <tbody>
                    {sites.map((site) => (
                      <tr key={site.siteIndex}>
                        <td>{site.symbol}</td>
                        <td>
                          <input
                            type="number"
                            min="1"
                            step="2"
                            disabled={Boolean(session)}
                            value={displaySiteValue(site, controls, "npt")}
                            onChange={(event) => updateSiteOverride(site.siteIndex, "npt", event.target.value)}
                          />
                        </td>
                        <td>
                          <input
                            type="number"
                            min="0.000001"
                            step="0.000001"
                            disabled={Boolean(session)}
                            value={displaySiteValue(site, controls, "r0")}
                            onChange={(event) => updateSiteOverride(site.siteIndex, "r0", event.target.value)}
                          />
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
                <h4>Structure preview</h4>
                <pre className="wien2k-struct-preview">{lastResult?.candidateStruct ?? draft.structContent}</pre>
                {lastResult?.nativeOutput && (
                  <>
                    <h4>Native diagnostic output</h4>
                    <pre className="wien2k-struct-preview">{lastResult.nativeOutput}</pre>
                  </>
                )}
              </details>
            </>
          )}
        </section>
        <section className="wien2k-preview-column">
          <div className="wien2k-viewer-frame">
            <UnitCellViewer crystalData={displayedCrystalData} />
          </div>
          <LiveOutputPanel
            title="WIEN2k structure refinement"
            output={outputLines.join("\n")}
            totalLineCount={outputLines.length}
            visibleLineCount={outputLines.length}
            panelClassName="output-panel wien2k-inline-output"
            outputClassName="output-text wien2k-inline-output-text"
          />
        </section>
      </div>
      <div className="step-actions">
        {session && (
          <button type="button" disabled={isRunning || isSaving} onClick={() => void discardSessionAndReturn()}>
            Discard Session
          </button>
        )}
        {draft && !session && (
          <button type="button" className="run-btn" disabled={!canOperate || isRunning} onClick={() => void runStage("rmt")}>
            Start RMT Refinement
          </button>
        )}
        {session?.phase === "rmt_ready" && (
          <>
            <button type="button" disabled={isRunning} onClick={() => void runStage("rmt")}>Re-run RMT</button>
            <button type="button" className="run-btn" disabled={isRunning || Boolean(controlsError)} onClick={() => void runStage("sgroup")}>
              Accept RMT and Run SGROUP
            </button>
          </>
        )}
        {session?.phase === "sgroup_ready" && (
          <>
            <button type="button" disabled={isRunning} onClick={() => void runStage("sgroup")}>Re-run SGROUP</button>
            <button type="button" className="run-btn" disabled={isRunning || Boolean(controlsError)} onClick={() => void runStage("symmetry")}>
              Accept SGROUP and Run SYMMETRY
            </button>
          </>
        )}
        {session?.phase === "symmetry_ready" && (
          <>
            {!lastResult?.saveAllowed && (
              <button type="button" disabled={isRunning || isSaving} onClick={() => void runStage("symmetry")}>
                Retry SYMMETRY
              </button>
            )}
            <button type="button" className="run-btn" disabled={isSaving || !lastResult?.saveAllowed} onClick={() => void saveSource()}>
              {isSaving ? "Saving..." : "Save Structure Source"}
            </button>
          </>
        )}
      </div>
    </div>
  );
}
