import { useMemo, useState } from "react";
import { invoke } from "@tauri-apps/api/core";
import { open, save } from "@tauri-apps/plugin-dialog";
import { readTextFile, writeTextFile } from "@tauri-apps/plugin-fs";
import { parseCIF } from "../lib/cifParser";
import {
  CifAtomElementGroup,
  CifSubstitutionMapping,
  CifSubstitutionResult,
  inspectCifAtomGroups,
  substituteCifElements,
} from "../lib/cifSubstitution";
import { isElementSymbol, normalizeElementSymbol } from "../lib/elements";
import { CrystalData } from "../lib/types";
import { InfoTooltip } from "./InfoTooltip";

interface CifVariant {
  id: string;
  filename: string;
  formula: string;
  added_at: string;
  calculations: unknown[];
}

interface CifSubstitutionDialogProps {
  isOpen: boolean;
  projectId: string;
  onClose: () => void;
  onSaved: (variant: CifVariant, result: CifSubstitutionResult) => Promise<void> | void;
}

interface MappingDraft {
  id: string;
  from: string;
  to: string;
}

function basename(path: string): string {
  return path.split("/").pop() || "structure.cif";
}

function formulaLabel(data: CrystalData | null): string {
  return data?.formula_sum || data?.formula_structural || "Unknown";
}

function createInitialDraft(groups: CifAtomElementGroup[]): MappingDraft {
  return {
    id: `mapping_${Date.now()}_${Math.random().toString(16).slice(2)}`,
    from: groups[0]?.element ?? "",
    to: "",
  };
}

function normalizeDrafts(drafts: MappingDraft[]): CifSubstitutionMapping[] {
  return drafts
    .map((draft) => ({
      from: normalizeElementSymbol(draft.from),
      to: normalizeElementSymbol(draft.to),
    }))
    .filter((mapping) => mapping.from && mapping.to);
}

export function CifSubstitutionDialog({
  isOpen,
  projectId,
  onClose,
  onSaved,
}: CifSubstitutionDialogProps) {
  const [sourcePath, setSourcePath] = useState("");
  const [sourceFilename, setSourceFilename] = useState("");
  const [sourceContent, setSourceContent] = useState("");
  const [sourceCrystalData, setSourceCrystalData] = useState<CrystalData | null>(null);
  const [groups, setGroups] = useState<CifAtomElementGroup[]>([]);
  const [drafts, setDrafts] = useState<MappingDraft[]>([]);
  const [error, setError] = useState<string | null>(null);
  const [status, setStatus] = useState<string | null>(null);
  const [isSelecting, setIsSelecting] = useState(false);
  const [isSaving, setIsSaving] = useState(false);
  const [isExporting, setIsExporting] = useState(false);
  const [savedResult, setSavedResult] = useState<CifSubstitutionResult | null>(null);

  const substitutionComputation = useMemo((): { result: CifSubstitutionResult | null; error: string | null } => {
    if (!sourceContent) return { result: null, error: null };
    const mappings = normalizeDrafts(drafts);
    if (mappings.length === 0) return { result: null, error: null };
    try {
      const result = substituteCifElements(sourceContent, mappings, { sourceFilename });
      return { result, error: null };
    } catch (e) {
      return { result: null, error: String(e instanceof Error ? e.message : e) };
    }
  }, [drafts, sourceContent, sourceFilename]);

  const substitutionResult = substitutionComputation.result;

  const substitutionError = useMemo(() => {
    if (!sourceContent) return null;
    const hasReplacement = drafts.some((draft) => draft.to.trim().length > 0);
    if (!hasReplacement) return null;
    try {
      normalizeDrafts(drafts).forEach((mapping) => {
        if (!isElementSymbol(mapping.to)) {
          throw new Error(`"${mapping.to}" is not a valid element symbol.`);
        }
      });
      return substitutionComputation.error;
    } catch (e) {
      return String(e instanceof Error ? e.message : e);
    }
  }, [drafts, sourceContent, substitutionComputation.error]);

  if (!isOpen) return null;

  async function handleSelectSource() {
    setIsSelecting(true);
    setError(null);
    setStatus(null);
    setSavedResult(null);
    try {
      const selected = await open({
        multiple: false,
        filters: [{ name: "CIF Files", extensions: ["cif"] }],
        title: "Select Source CIF File",
      });

      if (!selected || typeof selected !== "string") {
        return;
      }

      const content = await readTextFile(selected);
      const parsed = parseCIF(content);
      const inspectedGroups = inspectCifAtomGroups(content);
      if (inspectedGroups.length === 0) {
        throw new Error("No atom-site records were found in this CIF.");
      }

      setSourcePath(selected);
      setSourceFilename(basename(selected));
      setSourceContent(content);
      setSourceCrystalData(parsed);
      setGroups(inspectedGroups);
      setDrafts([createInitialDraft(inspectedGroups)]);
    } catch (e) {
      console.error("Failed to select source CIF:", e);
      setError(`Failed to load source CIF: ${e}`);
    } finally {
      setIsSelecting(false);
    }
  }

  function updateDraft(id: string, patch: Partial<MappingDraft>) {
    setSavedResult(null);
    setDrafts((current) => current.map((draft) => (
      draft.id === id ? { ...draft, ...patch } : draft
    )));
  }

  function addDraft() {
    setSavedResult(null);
    const used = new Set(drafts.map((draft) => draft.from));
    const nextGroup = groups.find((group) => !used.has(group.element)) ?? groups[0];
    setDrafts((current) => [...current, createInitialDraft(nextGroup ? [nextGroup] : groups)]);
  }

  function removeDraft(id: string) {
    setSavedResult(null);
    setDrafts((current) => current.length > 1 ? current.filter((draft) => draft.id !== id) : current);
  }

  async function handleSaveToProject() {
    if (!substitutionResult || substitutionError) return;
    setIsSaving(true);
    setError(null);
    setStatus(null);
    try {
      const parsedData = parseCIF(substitutionResult.content);
      const formula = parsedData.formula_sum || parsedData.formula_structural || substitutionResult.newFormula || "Unknown";
      const variant = await invoke<CifVariant>("add_cif_to_project", {
        projectId,
        cifData: {
          filename: substitutionResult.suggestedFilename,
          formula,
          content: substitutionResult.content,
          crystal_data: parsedData,
        },
      });
      setSavedResult(substitutionResult);
      setStatus(`Saved ${formula} to project.`);
      await onSaved(variant, substitutionResult);
    } catch (e) {
      console.error("Failed to save substituted CIF:", e);
      setError(`Failed to save modified CIF: ${e}`);
    } finally {
      setIsSaving(false);
    }
  }

  async function handleExport() {
    const result = savedResult ?? substitutionResult;
    if (!result) return;
    setIsExporting(true);
    setError(null);
    try {
      const destinationPath = await save({
        title: "Export Modified CIF",
        defaultPath: result.suggestedFilename,
        filters: [{ name: "CIF Files", extensions: ["cif"] }],
      });
      if (!destinationPath) return;
      await writeTextFile(destinationPath, result.content);
      setStatus(`Exported modified CIF to ${destinationPath}`);
    } catch (e) {
      console.error("Failed to export modified CIF:", e);
      setError(`Failed to export modified CIF: ${e}`);
    } finally {
      setIsExporting(false);
    }
  }

  const sourceFormula = formulaLabel(sourceCrystalData);
  const canAddDraft = groups.length > drafts.length;
  const canSave = Boolean(substitutionResult) && !substitutionError && !isSaving;

  return (
    <div className="dialog-overlay" onClick={onClose}>
      <div className="dialog-content cif-substitution-dialog" onClick={(e) => e.stopPropagation()}>
        <div className="dialog-header">
          <div className="cif-substitution-title">
            <h2>Modify Existing CIF</h2>
            <InfoTooltip text="Element substitutions are structural guesses. Not all substitutions produce a physically valid structure, so run structure relaxation before using results.">
              <span className="cif-substitution-help" aria-label="Substitution warning">?</span>
            </InfoTooltip>
          </div>
          <button className="dialog-close" onClick={onClose} aria-label="Close">
            &times;
          </button>
        </div>

        <div className="dialog-body cif-substitution-body">
          {error && <div className="dialog-error">{error}</div>}
          {status && <div className="info-banner">{status}</div>}

          <section className="cif-substitution-section">
            <div className="cif-substitution-section-header">
              <h3>Source CIF</h3>
              <button
                type="button"
                className="add-structure-btn"
                onClick={handleSelectSource}
                disabled={isSelecting || isSaving}
              >
                {sourceContent ? "Choose Different CIF" : (isSelecting ? "Selecting..." : "Select Source CIF")}
              </button>
            </div>
            {sourceContent ? (
              <div className="cif-substitution-source">
                <span>{sourceFilename}</span>
                <span>{sourceFormula}</span>
                <code>{sourcePath}</code>
              </div>
            ) : (
              <p className="cif-substitution-empty">Select a CIF to inspect its atom-site records.</p>
            )}
          </section>

          {groups.length > 0 && (
            <>
              <section className="cif-substitution-section">
                <h3>Atomic Sites</h3>
                <div className="cif-substitution-groups">
                  {groups.map((group) => (
                    <div key={group.element} className="cif-substitution-group">
                      <strong>{group.element}</strong>
                      <span>{group.count} site{group.count === 1 ? "" : "s"}</span>
                      <small>
                        {group.sites.slice(0, 3).map((site) => site.label || site.typeSymbol).join(", ")}
                        {group.sites.length > 3 ? "..." : ""}
                      </small>
                    </div>
                  ))}
                </div>
              </section>

              <section className="cif-substitution-section">
                <div className="cif-substitution-section-header">
                  <h3>Substitutions</h3>
                  <button
                    type="button"
                    className="add-structure-btn"
                    onClick={addDraft}
                    disabled={!canAddDraft || isSaving}
                  >
                    Add Substitution
                  </button>
                </div>
                <div className="cif-substitution-mappings">
                  {drafts.map((draft) => {
                    const normalizedTo = normalizeElementSymbol(draft.to);
                    const invalid = draft.to.trim().length > 0 && !isElementSymbol(normalizedTo);
                    return (
                      <div key={draft.id} className="cif-substitution-mapping-row">
                        <label>
                          <span>Replace</span>
                          <select
                            value={draft.from}
                            onChange={(e) => updateDraft(draft.id, { from: e.target.value })}
                            disabled={isSaving}
                          >
                            {groups.map((group) => (
                              <option key={group.element} value={group.element}>
                                {group.element}
                              </option>
                            ))}
                          </select>
                        </label>
                        <label>
                          <span>With</span>
                          <input
                            type="text"
                            value={draft.to}
                            onChange={(e) => updateDraft(draft.id, { to: e.target.value })}
                            onBlur={() => updateDraft(draft.id, { to: normalizedTo })}
                            className={invalid ? "invalid" : ""}
                            placeholder="Tl"
                            disabled={isSaving}
                          />
                        </label>
                        <button
                          type="button"
                          className="dialog-btn cancel cif-substitution-remove"
                          onClick={() => removeDraft(draft.id)}
                          disabled={drafts.length <= 1 || isSaving}
                          aria-label="Remove substitution"
                        >
                          Remove
                        </button>
                      </div>
                    );
                  })}
                </div>
                {substitutionError && <div className="dialog-error">{substitutionError}</div>}
              </section>

              {substitutionResult && !substitutionError && (
                <section className="cif-substitution-section">
                  <h3>Preview</h3>
                  <div className="cif-substitution-preview-summary">
                    <span>{substitutionResult.originalFormula || sourceFormula}</span>
                    <span aria-hidden="true">-&gt;</span>
                    <strong>{substitutionResult.newFormula || "Modified structure"}</strong>
                  </div>
                  {substitutionResult.warnings.length > 0 && (
                    <div className="cif-substitution-warnings">
                      {substitutionResult.warnings.map((warning) => (
                        <span key={warning}>{warning}</span>
                      ))}
                    </div>
                  )}
                  <div className="cif-substitution-preview-table">
                    <table>
                      <thead>
                        <tr>
                          <th>Site</th>
                          <th>Type</th>
                          <th>Fractional Position</th>
                          <th>Occupancy</th>
                        </tr>
                      </thead>
                      <tbody>
                        {substitutionResult.changedSites.map((site) => (
                          <tr key={`${site.label}-${site.newLabel}-${site.fractX}-${site.fractY}-${site.fractZ}`}>
                            <td>{site.label} -&gt; {site.newLabel}</td>
                            <td>{site.typeSymbol || "N/A"} -&gt; {site.newTypeSymbol || "N/A"}</td>
                            <td>{site.fractX}, {site.fractY}, {site.fractZ}</td>
                            <td>{site.occupancy || "N/A"}</td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                </section>
              )}
            </>
          )}
        </div>

        <div className="dialog-footer">
          <button className="dialog-btn cancel" onClick={onClose} disabled={isSaving || isExporting}>
            {savedResult ? "Done" : "Cancel"}
          </button>
          {savedResult && (
            <button
              className="dialog-btn cancel"
              onClick={handleExport}
              disabled={isSaving || isExporting}
            >
              {isExporting ? "Exporting..." : "Export CIF"}
            </button>
          )}
          <button
            className="dialog-btn save"
            onClick={handleSaveToProject}
            disabled={!canSave}
          >
            {isSaving ? "Saving..." : "Save to Project"}
          </button>
        </div>
      </div>
    </div>
  );
}
