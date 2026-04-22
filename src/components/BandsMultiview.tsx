import {
  CSSProperties,
  KeyboardEvent,
  PointerEvent as ReactPointerEvent,
  useCallback,
  useEffect,
  useMemo,
  useRef,
  useState,
} from "react";
import { invoke } from "@tauri-apps/api/core";
import { listen } from "@tauri-apps/api/event";
import {
  BandData,
  BandPlot,
  buildBandPlotProjectionOptions,
  BandPlotSharedSettings,
  FermiReferenceMode,
  getDefaultBandPlotEnergyRange,
} from "./BandPlot";
import { InfoTooltip } from "./InfoTooltip";

export interface BandsMultiviewCalculation {
  folder_id?: string | null;
  folder_name?: string | null;
  project_id: string;
  project_name: string;
  cif_id: string;
  cif_filename: string;
  cif_formula: string;
  calc_id: string;
  parameters: Record<string, unknown>;
  tags?: string[];
  started_at: string;
  completed_at: string;
  storage_bytes?: number | null;
  band_data: BandData;
  scf_fermi_energy?: number | null;
}

interface BandsMultiviewProps {
  onBack: () => void;
  initialCalculations?: BandsMultiviewCalculation[] | null;
}

export interface BandsMultiviewScanProgress {
  progress_event_id: string;
  phase: string;
  found_count: number;
  scanned_projects: number;
  total_projects: number;
}

interface MultiviewProjectGroup {
  projectId: string;
  projectName: string;
  calculations: BandsMultiviewCalculation[];
}

interface MultiviewFolderGroup {
  folderId: string | null;
  folderName: string;
  projects: MultiviewProjectGroup[];
}

type CalcTagType = "info" | "feature";

interface TileDragSession {
  calcId: string;
  width: number;
  height: number;
  offsetX: number;
  offsetY: number;
  startClientX: number;
  startClientY: number;
}

const DEFAULT_SHARED_SETTINGS: BandPlotSharedSettings = {
  fermiReferenceMode: "scf",
  lineWidth: 1.5,
  lineOpacity: 0.85,
  plotTextScale: 1,
  colorMode: "rainbow",
  singleBandColor: "#1565c0",
  rainbowPalette: "jet",
  plotBgWhite: true,
  showBandGapOverlay: true,
};
const DEFAULT_MULTIVIEW_PLOT_HEIGHT = 430;
const DEFAULT_MULTIVIEW_PLOT_WIDTH = 520;
const DEFAULT_PROJECTION_OPTIONS = [{ value: "none", label: "none" }];

function createProgressEventId(): string {
  if (typeof crypto !== "undefined" && typeof crypto.randomUUID === "function") {
    return crypto.randomUUID();
  }
  return `bands-multiview-${Date.now()}-${Math.random().toString(36).slice(2)}`;
}

function formatDate(isoString: string): string {
  try {
    return new Date(isoString).toLocaleDateString(undefined, {
      year: "numeric",
      month: "short",
      day: "numeric",
      hour: "2-digit",
      minute: "2-digit",
    });
  } catch {
    return isoString;
  }
}

function formatAxisInputValue(value: number): string {
  if (!Number.isFinite(value)) return "";
  return Number.parseFloat(value.toFixed(6)).toString();
}

function getParamText(params: Record<string, unknown>, key: string): string | null {
  const value = params[key];
  if (value == null) return null;
  const text = String(value).trim();
  return text.length > 0 ? text : null;
}

function isHpcCalculation(calc: Pick<BandsMultiviewCalculation, "parameters">): boolean {
  const params = calc.parameters || {};
  const backend = String(params.execution_backend || "").trim().toLowerCase();
  if (backend === "hpc") {
    return true;
  }
  if (params.remote_job_id || params.remote_workdir || params.remote_project_path) {
    return true;
  }
  return false;
}

function getCalcTagClass(tag: { label: string; type: CalcTagType }): string {
  const normalizedLabel = tag.label.trim().toUpperCase();
  const isHpcTag = normalizedLabel === "HPC";
  return `calc-tag calc-tag-${tag.type}${isHpcTag ? " calc-tag-hpc" : ""}`;
}

function getBandsTags(calc: BandsMultiviewCalculation): { label: string; type: CalcTagType }[] {
  const tags: { label: string; type: CalcTagType }[] = [];
  const params = calc.parameters || {};
  const pushTag = (label: string, type: CalcTagType) => {
    if (!tags.some((tag) => tag.label === label)) {
      tags.push({ label, type });
    }
  };

  if (params.total_k_points) {
    pushTag(`${params.total_k_points} k-pts`, "info");
  }
  if (params.lspinorb) {
    pushTag("SOC", "feature");
  }
  if (params.nspin === 2) {
    pushTag("Magnetic", "feature");
  }
  if (params.lda_plus_u) {
    pushTag("DFT+U", "feature");
  }
  if (params.vdw_corr && params.vdw_corr !== "none") {
    pushTag("vdW", "feature");
  }

  const hasProjectionTag = Array.isArray(calc.tags)
    && calc.tags.some((tag) => {
      const normalized = String(tag).trim().toLowerCase();
      return normalized === "proj" || normalized === "orb";
    });
  if (hasProjectionTag || params.fat_bands_requested) {
    pushTag("Proj", "feature");
  }

  if (isHpcCalculation(calc)) {
    pushTag("HPC", "feature");
  }

  return tags;
}

function buildFolderGroups(calculations: BandsMultiviewCalculation[]): MultiviewFolderGroup[] {
  const folders = new Map<string, MultiviewFolderGroup>();

  for (const calculation of calculations) {
    const folderId = calculation.folder_id ?? null;
    const folderKey = folderId ?? "__root__";
    const folderName = calculation.folder_name?.trim() || "Root";
    let folderGroup = folders.get(folderKey);
    if (!folderGroup) {
      folderGroup = {
        folderId,
        folderName,
        projects: [],
      };
      folders.set(folderKey, folderGroup);
    }

    let projectGroup = folderGroup.projects.find(
      (project) => project.projectId === calculation.project_id,
    );
    if (!projectGroup) {
      projectGroup = {
        projectId: calculation.project_id,
        projectName: calculation.project_name,
        calculations: [],
      };
      folderGroup.projects.push(projectGroup);
    }
    projectGroup.calculations.push(calculation);
  }

  const groups = Array.from(folders.values());
  for (const folder of groups) {
    folder.projects.sort((a, b) => a.projectName.localeCompare(b.projectName));
    for (const project of folder.projects) {
      project.calculations.sort((a, b) => b.completed_at.localeCompare(a.completed_at));
    }
  }

  groups.sort((a, b) => {
    if (a.folderId === null && b.folderId !== null) return -1;
    if (a.folderId !== null && b.folderId === null) return 1;
    return a.folderName.localeCompare(b.folderName);
  });
  return groups;
}

export function BandsMultiview({
  onBack,
  initialCalculations,
}: BandsMultiviewProps) {
  const [calculations, setCalculations] = useState<BandsMultiviewCalculation[]>(
    initialCalculations ?? [],
  );
  const [isLoading, setIsLoading] = useState(initialCalculations === undefined);
  const [error, setError] = useState<string | null>(null);
  const [loadingProgress, setLoadingProgress] = useState<BandsMultiviewScanProgress>({
    progress_event_id: "",
    phase: "loading",
    found_count: 0,
    scanned_projects: 0,
    total_projects: 0,
  });
  const [showSelector, setShowSelector] = useState(true);
  const [selectedIds, setSelectedIds] = useState<string[]>([]);
  const [sharedSettings, setSharedSettings] = useState<BandPlotSharedSettings>(
    DEFAULT_SHARED_SETTINGS,
  );
  const [manualRange, setManualRange] = useState<[number, number] | null>(null);
  const [manualRangeMinInput, setManualRangeMinInput] = useState("");
  const [manualRangeMaxInput, setManualRangeMaxInput] = useState("");
  const [manualRangeError, setManualRangeError] = useState<string | null>(null);
  const [appearanceExpanded, setAppearanceExpanded] = useState(true);
  const [bandGapExpanded, setBandGapExpanded] = useState(true);
  const [plotHeight, setPlotHeight] = useState(DEFAULT_MULTIVIEW_PLOT_HEIGHT);
  const [plotWidth, setPlotWidth] = useState(DEFAULT_MULTIVIEW_PLOT_WIDTH);
  const [projectionSelectionsByCalcId, setProjectionSelectionsByCalcId] = useState<
    Record<string, string>
  >({});
  const [tileDrag, setTileDrag] = useState<TileDragSession | null>(null);
  const [dragTargetCalcId, setDragTargetCalcId] = useState<string | null>(null);
  const cardRefs = useRef(new Map<string, HTMLElement>());
  const floatingCardRef = useRef<HTMLDivElement | null>(null);
  const selectedIdsRef = useRef<string[]>([]);
  const dragPointerRef = useRef<{ x: number; y: number } | null>(null);
  const dragTargetCalcIdRef = useRef<string | null>(null);

  useEffect(() => {
    if (initialCalculations === undefined) {
      return;
    }
    setCalculations(initialCalculations ?? []);
    setIsLoading(false);
    setError(null);
  }, [initialCalculations]);

  useEffect(() => {
    if (initialCalculations !== undefined) {
      return;
    }
    let cancelled = false;
    let listenerRemoved = false;
    const progressEventId = createProgressEventId();
    const unlistenPromise = listen<BandsMultiviewScanProgress>(
      "multiview-bands-progress",
      (event) => {
        const payload = event.payload;
        if (cancelled || payload.progress_event_id !== progressEventId) {
          return;
        }
        setLoadingProgress(payload);
      },
    );

    const cleanupListener = () => {
      if (listenerRemoved) {
        return;
      }
      listenerRemoved = true;
      void unlistenPromise.then((unlisten) => unlisten());
    };

    async function loadCalculations() {
      setIsLoading(true);
      setError(null);
      setLoadingProgress({
        progress_event_id: progressEventId,
        phase: "loading",
        found_count: 0,
        scanned_projects: 0,
        total_projects: 0,
      });
      try {
        const result = await invoke<BandsMultiviewCalculation[]>(
          "list_multiview_band_calculations",
          { progressEventId },
        );
        if (cancelled) return;
        setCalculations(result);
      } catch (e) {
        if (cancelled) return;
        console.error("Failed to load multiview band calculations:", e);
        setError(`Failed to load saved band calculations: ${e}`);
      } finally {
        if (!cancelled) {
          setIsLoading(false);
        }
        cleanupListener();
      }
    }

    void loadCalculations();

    return () => {
      cancelled = true;
      cleanupListener();
    };
  }, [initialCalculations]);

  const folderGroups = useMemo(() => buildFolderGroups(calculations), [calculations]);
  const calculationById = useMemo(
    () => new Map(calculations.map((calc) => [calc.calc_id, calc])),
    [calculations],
  );
  const selectedIdSet = useMemo(() => new Set(selectedIds), [selectedIds]);
  const selectedCalculations = useMemo(
    () => selectedIds
      .map((calcId) => calculationById.get(calcId))
      .filter((calc): calc is BandsMultiviewCalculation => Boolean(calc)),
    [calculationById, selectedIds],
  );
  const projectionOptionsByCalcId = useMemo(() => {
    const options = new Map<string, { value: string; label: string }[]>();
    for (const calc of selectedCalculations) {
      options.set(calc.calc_id, buildBandPlotProjectionOptions(calc.band_data));
    }
    return options;
  }, [selectedCalculations]);
  const selectedCalculationCards = useMemo(() => selectedCalculations.map((calc) => {
    const tags = getBandsTags(calc);
    const kPath = getParamText(calc.parameters, "k_path");
    const projectionOptions =
      projectionOptionsByCalcId.get(calc.calc_id) ?? DEFAULT_PROJECTION_OPTIONS;
    const requestedProjectionSelection =
      projectionSelectionsByCalcId[calc.calc_id] ?? "none";
    const projectionSelection = projectionOptions.some(
      (option) => option.value === requestedProjectionSelection,
    )
      ? requestedProjectionSelection
      : "none";

    return {
      calc,
      tags,
      kPath,
      projectionOptions,
      projectionSelection,
    };
  }), [projectionOptionsByCalcId, projectionSelectionsByCalcId, selectedCalculations]);
  const draggingCalcId = tileDrag?.calcId ?? null;
  const draggedCalculationCard = useMemo(
    () => selectedCalculationCards.find((card) => card.calc.calc_id === draggingCalcId) ?? null,
    [draggingCalcId, selectedCalculationCards],
  );
  const previewGridStyle = useMemo(
    () => ({ "--bands-multiview-plot-width": `${plotWidth}px` } as CSSProperties),
    [plotWidth],
  );
  const loadingProjectCount = loadingProgress.total_projects > 0
    ? Math.min(loadingProgress.scanned_projects, loadingProgress.total_projects)
    : 0;

  useEffect(() => {
    setSelectedIds((current) => current.filter((calcId) => calculationById.has(calcId)));
  }, [calculationById]);

  useEffect(() => {
    selectedIdsRef.current = selectedIds;
  }, [selectedIds]);

  useEffect(() => {
    dragTargetCalcIdRef.current = dragTargetCalcId;
  }, [dragTargetCalcId]);

  useEffect(() => {
    if (!tileDrag) {
      return;
    }
    if (showSelector || !selectedIdSet.has(tileDrag.calcId)) {
      setTileDrag(null);
      setDragTargetCalcId(null);
      dragTargetCalcIdRef.current = null;
      dragPointerRef.current = null;
      document.body.classList.remove("bands-multiview-tile-dragging");
    }
  }, [selectedIdSet, showSelector, tileDrag]);

  const autoEnergyRange = useMemo<[number, number]>(() => {
    if (selectedCalculations.length === 0) {
      return [-5, 5];
    }

    let globalMin = Number.POSITIVE_INFINITY;
    let globalMax = Number.NEGATIVE_INFINITY;
    for (const calc of selectedCalculations) {
      const [calcMin, calcMax] = getDefaultBandPlotEnergyRange(
        calc.band_data,
        calc.scf_fermi_energy ?? undefined,
        sharedSettings.fermiReferenceMode,
      );
      globalMin = Math.min(globalMin, calcMin);
      globalMax = Math.max(globalMax, calcMax);
    }

    if (!Number.isFinite(globalMin) || !Number.isFinite(globalMax) || globalMax <= globalMin) {
      return [-5, 5];
    }
    return [globalMin, globalMax];
  }, [selectedCalculations, sharedSettings.fermiReferenceMode]);

  const resolvedEnergyRange = manualRange ?? autoEnergyRange;

  useEffect(() => {
    setManualRangeMinInput(formatAxisInputValue(resolvedEnergyRange[0]));
    setManualRangeMaxInput(formatAxisInputValue(resolvedEnergyRange[1]));
  }, [resolvedEnergyRange]);

  function toggleSelectedCalculation(calcId: string) {
    setSelectedIds((current) => (
      current.includes(calcId)
        ? current.filter((id) => id !== calcId)
        : [...current, calcId]
    ));
  }

  function handleSelectAll() {
    setSelectedIds(calculations.map((calc) => calc.calc_id));
  }

  function handleClearSelection() {
    setSelectedIds([]);
  }

  const setCardRef = useCallback((calcId: string, node: HTMLElement | null) => {
    if (node) {
      cardRefs.current.set(calcId, node);
    } else {
      cardRefs.current.delete(calcId);
    }
  }, []);

  const reorderSelectedCalculations = useCallback((movedCalcId: string, targetCalcId: string) => {
    if (movedCalcId === targetCalcId) {
      return;
    }

    setSelectedIds((current) => {
      const movedIndex = current.indexOf(movedCalcId);
      const targetIndex = current.indexOf(targetCalcId);
      if (movedIndex === -1 || targetIndex === -1) {
        return current;
      }

      const next = [...current];
      const [movedId] = next.splice(movedIndex, 1);
      if (!movedId) {
        return current;
      }

      const insertionIndex = movedIndex < targetIndex ? targetIndex : targetIndex;
      next.splice(insertionIndex, 0, movedId);
      return next;
    });
  }, []);

  const updateFloatingCardPosition = useCallback((
    clientX: number,
    clientY: number,
    session: TileDragSession,
  ) => {
    const floatingCard = floatingCardRef.current;
    if (!floatingCard) {
      return;
    }
    const left = clientX - session.offsetX;
    const top = clientY - session.offsetY;
    floatingCard.style.transform = `translate3d(${left}px, ${top}px, 0)`;
  }, []);

  const findNearestCalculationCard = useCallback((
    clientX: number,
    clientY: number,
    movedCalcId: string,
  ): string | null => {
    let bestCalcId: string | null = null;
    let bestScore = Number.POSITIVE_INFINITY;

    for (const calcId of selectedIdsRef.current) {
      if (calcId === movedCalcId) {
        continue;
      }
      const card = cardRefs.current.get(calcId);
      if (!card) {
        continue;
      }
      const rect = card.getBoundingClientRect();
      const inside =
        clientX >= rect.left
        && clientX <= rect.right
        && clientY >= rect.top
        && clientY <= rect.bottom;
      const dx = clientX - (rect.left + rect.width / 2);
      const dy = clientY - (rect.top + rect.height / 2);
      const score = inside ? 0 : (dx * dx) + (dy * dy);
      if (score < bestScore) {
        bestScore = score;
        bestCalcId = calcId;
      }
    }

    return bestCalcId;
  }, []);

  function applyManualRange() {
    const parsedMin = Number.parseFloat(manualRangeMinInput.trim());
    const parsedMax = Number.parseFloat(manualRangeMaxInput.trim());
    if (!Number.isFinite(parsedMin) || !Number.isFinite(parsedMax)) {
      setManualRangeError("Enter valid numeric values for Y min and Y max.");
      return;
    }
    if (parsedMax <= parsedMin) {
      setManualRangeError("Y max must be greater than Y min.");
      return;
    }
    setManualRange([parsedMin, parsedMax]);
    setManualRangeError(null);
  }

  function resetEnergyRange() {
    setManualRange(null);
    setManualRangeError(null);
  }

  function handleManualRangeKeyDown(event: KeyboardEvent<HTMLInputElement>) {
    if (event.key !== "Enter") return;
    event.preventDefault();
    applyManualRange();
  }

  function openPreview() {
    if (selectedIds.length === 0) return;
    setShowSelector(false);
  }

  function handleExportPdf() {
    if (typeof window === "undefined" || selectedCalculations.length === 0) {
      return;
    }
    window.print();
  }

  const finishTileDrag = useCallback(() => {
    setTileDrag(null);
    setDragTargetCalcId(null);
    dragTargetCalcIdRef.current = null;
    dragPointerRef.current = null;
    document.body.classList.remove("bands-multiview-tile-dragging");
  }, []);

  function handleTileHandlePointerDown(
    event: ReactPointerEvent<HTMLButtonElement>,
    calcId: string,
  ) {
    if (selectedCalculations.length < 2 || event.button !== 0) {
      return;
    }
    const card = cardRefs.current.get(calcId);
    if (!card) {
      return;
    }
    const rect = card.getBoundingClientRect();
    dragPointerRef.current = { x: event.clientX, y: event.clientY };
    setTileDrag({
      calcId,
      width: rect.width,
      height: rect.height,
      offsetX: event.clientX - rect.left,
      offsetY: event.clientY - rect.top,
      startClientX: event.clientX,
      startClientY: event.clientY,
    });
    setDragTargetCalcId(calcId);
    dragTargetCalcIdRef.current = calcId;
    event.preventDefault();
  }

  useEffect(() => {
    if (!tileDrag) {
      return;
    }

    document.body.classList.add("bands-multiview-tile-dragging");
    const initialPointer = dragPointerRef.current;
    if (initialPointer) {
      requestAnimationFrame(() => {
        updateFloatingCardPosition(initialPointer.x, initialPointer.y, tileDrag);
      });
    }

    const handlePointerMove = (event: PointerEvent) => {
      dragPointerRef.current = { x: event.clientX, y: event.clientY };
      updateFloatingCardPosition(event.clientX, event.clientY, tileDrag);

      const dragDistance = Math.hypot(
        event.clientX - tileDrag.startClientX,
        event.clientY - tileDrag.startClientY,
      );
      if (dragDistance < 10) {
        return;
      }

      const nearestCalcId = findNearestCalculationCard(
        event.clientX,
        event.clientY,
        tileDrag.calcId,
      );
      if (!nearestCalcId || nearestCalcId === dragTargetCalcIdRef.current) {
        return;
      }

      dragTargetCalcIdRef.current = nearestCalcId;
      setDragTargetCalcId(nearestCalcId);
      reorderSelectedCalculations(tileDrag.calcId, nearestCalcId);
    };

    const handlePointerUp = () => {
      finishTileDrag();
    };

    window.addEventListener("pointermove", handlePointerMove);
    window.addEventListener("pointerup", handlePointerUp);
    window.addEventListener("pointercancel", handlePointerUp);

    return () => {
      window.removeEventListener("pointermove", handlePointerMove);
      window.removeEventListener("pointerup", handlePointerUp);
      window.removeEventListener("pointercancel", handlePointerUp);
      document.body.classList.remove("bands-multiview-tile-dragging");
    };
  }, [
    findNearestCalculationCard,
    finishTileDrag,
    reorderSelectedCalculations,
    tileDrag,
    updateFloatingCardPosition,
  ]);

  return (
    <div className="bands-viewer-container bands-multiview-screen">
      <div className="bands-viewer-header bands-multiview-header">
        <button className="back-button" onClick={onBack}>
          ← Back to Projects
        </button>
        <h2>Band Multiview</h2>
        <div className="bands-multiview-header-actions">
          {!showSelector && (
            <>
              <button
                type="button"
                className="secondary-project-btn"
                onClick={() => setShowSelector(true)}
              >
                Change Selection
              </button>
              <button
                type="button"
                className="view-bands-btn"
                onClick={handleExportPdf}
                disabled={selectedCalculations.length === 0}
              >
                Export PDF
              </button>
            </>
          )}
        </div>
      </div>

      {error && <div className="error-banner">{error}</div>}

      {isLoading ? (
        <div className="browser-content bands-multiview-loading-shell">
          <div className="bands-multiview-loading-state">
            <div className="bands-multiview-loading-spinner" aria-hidden="true" />
            <div className="bands-multiview-loading-copy">
              <h3>
                Found {loadingProgress.found_count} band calculation
                {loadingProgress.found_count !== 1 ? "s" : ""}
                ...
              </h3>
              <p>
                {loadingProgress.total_projects > 0
                  ? `Scanned ${loadingProjectCount} of ${loadingProgress.total_projects} saved projects.`
                  : "Scanning saved projects for completed band calculations."}
              </p>
            </div>
          </div>
        </div>
      ) : calculations.length === 0 ? (
        <div className="browser-content">
          <div className="empty-state">
            <h3>No Saved Bands Calculations</h3>
            <p>Run and save a band structure calculation before opening multiview.</p>
          </div>
        </div>
      ) : showSelector ? (
        <div className="bands-multiview-selector-screen">
          <div className="bands-multiview-selector-toolbar">
            <div>
              <h3>Select Saved Band Calculations</h3>
              <p>
                Completed band runs are grouped by folder and project. Select as many as you want
                for the shared two-column preview.
              </p>
            </div>
            <div className="bands-multiview-selector-actions">
              <button
                type="button"
                className="secondary-project-btn"
                onClick={handleSelectAll}
              >
                Select All
              </button>
              <button
                type="button"
                className="secondary-project-btn"
                onClick={handleClearSelection}
                disabled={selectedIds.length === 0}
              >
                Clear
              </button>
              <button
                type="button"
                className="view-bands-btn"
                onClick={openPreview}
                disabled={selectedIds.length === 0}
              >
                Open Preview ({selectedIds.length})
              </button>
            </div>
          </div>

          <div className="bands-multiview-selection-list">
            {folderGroups.map((folder) => (
              <section key={folder.folderId ?? "root"} className="bands-multiview-folder-group">
                <div className="bands-multiview-folder-header">
                  <h3>{folder.folderName}</h3>
                  <span>
                    {folder.projects.reduce((count, project) => count + project.calculations.length, 0)}
                    {" "}
                    calculations
                  </span>
                </div>

                {folder.projects.map((project) => (
                  <div key={project.projectId} className="bands-multiview-project-group">
                    <div className="bands-multiview-project-header">
                      <h4>{project.projectName}</h4>
                      <span>
                        {project.calculations.length} saved run
                        {project.calculations.length !== 1 ? "s" : ""}
                      </span>
                    </div>

                    <div className="bands-multiview-project-calculations">
                      {project.calculations.map((calc) => {
                        const isSelected = selectedIdSet.has(calc.calc_id);
                        const tags = getBandsTags(calc);
                        const kPath = getParamText(calc.parameters, "k_path");
                        const bandCount = getParamText(calc.parameters, "n_bands");
                        return (
                          <label
                            key={calc.calc_id}
                            className={`bands-multiview-select-card ${isSelected ? "selected" : ""}`}
                          >
                            <div className="bands-multiview-select-checkbox">
                              <input
                                type="checkbox"
                                checked={isSelected}
                                onChange={() => toggleSelectedCalculation(calc.calc_id)}
                              />
                            </div>
                            <div className="bands-multiview-select-content">
                              <div className="bands-multiview-select-header">
                                <div>
                                  <h5>{calc.cif_formula || calc.cif_filename}</h5>
                                  <p>{calc.cif_filename}</p>
                                </div>
                                <span>{formatDate(calc.completed_at)}</span>
                              </div>
                              <div className="bands-multiview-select-meta">
                                {kPath && (
                                  <span className="calc-kpath">{kPath}</span>
                                )}
                                {bandCount && (
                                  <span>{bandCount} bands</span>
                                )}
                              </div>
                              {tags.length > 0 && (
                                <div className="calc-tags">
                                  {tags.map((tag) => (
                                    <span key={tag.label} className={getCalcTagClass(tag)}>
                                      {tag.label}
                                    </span>
                                  ))}
                                </div>
                              )}
                            </div>
                          </label>
                        );
                      })}
                    </div>
                  </div>
                ))}
              </section>
            ))}
          </div>
        </div>
      ) : (
        <div className="bands-multiview-content">
          <div className="bands-multiview-preview-pane">
            <div className="bands-multiview-scroll-window">
              <div className="bands-multiview-grid" style={previewGridStyle}>
                {selectedCalculationCards.map((card) => {
                  const {
                    calc,
                    tags,
                    kPath,
                    projectionOptions,
                    projectionSelection,
                  } = card;
                  return (
                    <article
                      key={calc.calc_id}
                      ref={(node) => setCardRef(calc.calc_id, node)}
                      className={[
                        "bands-multiview-card",
                        draggingCalcId === calc.calc_id ? "is-dragging" : "",
                        draggingCalcId === calc.calc_id ? "is-drag-placeholder" : "",
                        dragTargetCalcId === calc.calc_id && draggingCalcId !== calc.calc_id
                          ? "is-drop-target"
                          : "",
                      ].filter(Boolean).join(" ")}
                      data-multiview-card-id={calc.calc_id}
                    >
                      <div className="bands-multiview-card-header">
                        <div>
                          <h3>{calc.project_name}</h3>
                          <p>
                            {calc.folder_name?.trim() || "Root"} · {calc.cif_formula || calc.cif_filename}
                          </p>
                        </div>
                        <div className="bands-multiview-card-header-side">
                          <span>{formatDate(calc.completed_at)}</span>
                          <InfoTooltip text="Drag to reorder">
                            <button
                              type="button"
                              className="bands-multiview-reorder-handle"
                              disabled={selectedCalculations.length < 2}
                              onPointerDown={(event) => handleTileHandlePointerDown(event, calc.calc_id)}
                              aria-label={`Drag to reorder ${calc.cif_formula || calc.cif_filename}`}
                            >
                              ⋮⋮
                            </button>
                          </InfoTooltip>
                        </div>
                      </div>

                      {(kPath || tags.length > 0) && (
                        <div className="bands-multiview-card-meta">
                          {kPath && (
                            <span className="calc-kpath">{kPath}</span>
                          )}
                          {tags.length > 0 && (
                            <div className="calc-tags">
                              {tags.map((tag) => (
                                <span key={tag.label} className={getCalcTagClass(tag)}>
                                  {tag.label}
                                </span>
                              ))}
                            </div>
                          )}
                        </div>
                      )}

                      <div className="bands-multiview-card-controls">
                        <label className="bands-multiview-projection-picker">
                          <span>Projection:</span>
                          <select
                            value={projectionSelection}
                            onChange={(event) =>
                              setProjectionSelectionsByCalcId((current) => ({
                                ...current,
                                [calc.calc_id]: event.target.value,
                              }))
                            }
                            disabled={projectionOptions.length <= 1}
                            aria-label={`Projection for ${calc.project_name}`}
                          >
                            {projectionOptions.map((option) => (
                              <option key={option.value} value={option.value}>
                                {option.label}
                              </option>
                            ))}
                          </select>
                        </label>
                      </div>

                      <div
                        className="bands-multiview-plot-frame"
                        style={{ height: `${plotHeight}px`, minHeight: `${plotHeight}px` }}
                      >
                        <BandPlot
                          data={calc.band_data}
                          width={plotWidth}
                          height={plotHeight}
                          energyRange={resolvedEnergyRange}
                          scfFermiEnergy={calc.scf_fermi_energy ?? undefined}
                          viewerType="electronic"
                          sharedSettings={sharedSettings}
                          showSidebar={false}
                          projectionSelection={projectionSelection}
                          enableWheelRangeControl={false}
                          enableHoverScrollLock={false}
                        />
                      </div>
                    </article>
                  );
                })}
              </div>
            </div>
          </div>

          <aside className="band-plot-sidebar bands-multiview-settings-sidebar">
            <div className="band-plot-controls">
              <button onClick={resetEnergyRange} className="band-plot-reset">
                Reset View
              </button>
              <button
                type="button"
                className="band-plot-export"
                onClick={handleExportPdf}
                disabled={selectedCalculations.length === 0}
              >
                Export PDF
              </button>
              <span className="band-plot-hint">
                Plot wheel zoom is disabled here so page scrolling always moves through the grid.
              </span>
              <span className="band-plot-hint">
                Drag a tile handle to reorder the preview stack.
              </span>
            </div>

            <div className="band-plot-control-panel">
              <div className="band-control-section">
                <button
                  type="button"
                  className="band-control-header"
                  onClick={() => setAppearanceExpanded((prev) => !prev)}
                >
                  <span className={`collapse-icon ${appearanceExpanded ? "expanded" : ""}`}>▶</span>
                  Appearance
                </button>
                {appearanceExpanded && (
                  <div className="band-control-grid">
                    <div className="band-control-row">
                      <label>Fermi Reference</label>
                      <select
                        value={sharedSettings.fermiReferenceMode}
                        onChange={(event) => {
                          setSharedSettings((current) => ({
                            ...current,
                            fermiReferenceMode: event.target.value as FermiReferenceMode,
                          }));
                          setManualRange(null);
                          setManualRangeError(null);
                        }}
                      >
                        <option value="scf">SCF</option>
                        <option value="bands">Bands run</option>
                      </select>
                    </div>

                    <div className="band-control-row">
                      <label>Y Range (eV)</label>
                      <div className="band-control-range-inputs">
                        <input
                          type="number"
                          step="any"
                          value={manualRangeMinInput}
                          onChange={(event) => {
                            setManualRangeMinInput(event.target.value);
                            setManualRangeError(null);
                          }}
                          onKeyDown={handleManualRangeKeyDown}
                          aria-label="Y minimum"
                        />
                        <span className="band-control-range-separator">to</span>
                        <input
                          type="number"
                          step="any"
                          value={manualRangeMaxInput}
                          onChange={(event) => {
                            setManualRangeMaxInput(event.target.value);
                            setManualRangeError(null);
                          }}
                          onKeyDown={handleManualRangeKeyDown}
                          aria-label="Y maximum"
                        />
                      </div>
                      <button
                        type="button"
                        className="band-control-apply"
                        onClick={applyManualRange}
                      >
                        Apply
                      </button>
                      {manualRangeError && (
                        <span className="band-control-range-error">{manualRangeError}</span>
                      )}
                    </div>

                    <div className="band-control-row">
                      <label>Line Thickness</label>
                      <input
                        type="range"
                        min={0.5}
                        max={5}
                        step={0.1}
                        value={sharedSettings.lineWidth}
                        onChange={(event) =>
                          setSharedSettings((current) => ({
                            ...current,
                            lineWidth: Number.parseFloat(event.target.value),
                          }))
                        }
                      />
                      <span className="band-control-value">{sharedSettings.lineWidth.toFixed(1)} px</span>
                    </div>

                    <div className="band-control-row">
                      <label>Line Opacity</label>
                      <input
                        type="range"
                        min={0.1}
                        max={1}
                        step={0.05}
                        value={sharedSettings.lineOpacity}
                        onChange={(event) =>
                          setSharedSettings((current) => ({
                            ...current,
                            lineOpacity: Number.parseFloat(event.target.value),
                          }))
                        }
                      />
                      <span className="band-control-value">{sharedSettings.lineOpacity.toFixed(2)}</span>
                    </div>

                    <div className="band-control-row">
                      <label>Plot Text Size</label>
                      <input
                        type="range"
                        min={0.7}
                        max={2}
                        step={0.05}
                        value={sharedSettings.plotTextScale}
                        onChange={(event) =>
                          setSharedSettings((current) => ({
                            ...current,
                            plotTextScale: Number.parseFloat(event.target.value),
                          }))
                        }
                      />
                      <span className="band-control-value">{sharedSettings.plotTextScale.toFixed(2)}x</span>
                    </div>

                    <div className="band-control-row">
                      <label>Plot Width</label>
                      <input
                        type="range"
                        min={320}
                        max={860}
                        step={10}
                        value={plotWidth}
                        onChange={(event) => setPlotWidth(Number.parseFloat(event.target.value))}
                      />
                      <span className="band-control-value">{plotWidth} u</span>
                    </div>

                    <div className="band-control-row">
                      <label>Plot Height</label>
                      <input
                        type="range"
                        min={280}
                        max={1200}
                        step={10}
                        value={plotHeight}
                        onChange={(event) => setPlotHeight(Number.parseFloat(event.target.value))}
                      />
                      <span className="band-control-value">{plotHeight} u</span>
                    </div>

                    <div className="band-control-row">
                      <label>Band Color Mode</label>
                      <select
                        value={sharedSettings.colorMode}
                        onChange={(event) =>
                          setSharedSettings((current) => ({
                            ...current,
                            colorMode: event.target.value as BandPlotSharedSettings["colorMode"],
                          }))
                        }
                      >
                        <option value="single">Single color</option>
                        <option value="rainbow">Rainbow by band</option>
                      </select>
                    </div>

                    {sharedSettings.colorMode === "single" && (
                      <div className="band-control-row">
                        <label>Band Color</label>
                        <input
                          type="color"
                          value={sharedSettings.singleBandColor}
                          onChange={(event) =>
                            setSharedSettings((current) => ({
                              ...current,
                              singleBandColor: event.target.value,
                            }))
                          }
                        />
                      </div>
                    )}

                    {sharedSettings.colorMode === "rainbow" && (
                      <div className="band-control-row">
                        <label>Rainbow Palette</label>
                        <select
                          value={sharedSettings.rainbowPalette}
                          onChange={(event) =>
                            setSharedSettings((current) => ({
                              ...current,
                              rainbowPalette: event.target.value as BandPlotSharedSettings["rainbowPalette"],
                            }))
                          }
                        >
                          <option value="jet">Jet-like</option>
                          <option value="sinebow">Sinebow</option>
                        </select>
                      </div>
                    )}

                    <div className="band-control-row">
                      <label>Plot Background</label>
                      <select
                        value={sharedSettings.plotBgWhite ? "white" : "theme"}
                        onChange={(event) =>
                          setSharedSettings((current) => ({
                            ...current,
                            plotBgWhite: event.target.value === "white",
                          }))
                        }
                      >
                        <option value="white">Always White</option>
                        <option value="theme">Match Theme</option>
                      </select>
                    </div>
                  </div>
                )}
              </div>

              <div className="band-control-section">
                <button
                  type="button"
                  className="band-control-header"
                  onClick={() => setBandGapExpanded((prev) => !prev)}
                >
                  <span className={`collapse-icon ${bandGapExpanded ? "expanded" : ""}`}>▶</span>
                  Band Gap
                </button>
                {bandGapExpanded && (
                  <div className="band-control-grid">
                    <div className="band-control-row">
                      <label>Show Gap Overlay</label>
                      <input
                        type="checkbox"
                        checked={sharedSettings.showBandGapOverlay}
                        onChange={(event) =>
                          setSharedSettings((current) => ({
                            ...current,
                            showBandGapOverlay: event.target.checked,
                          }))
                        }
                      />
                    </div>

                    <div className="band-control-note">
                      Gap detection and overlay styling stay synchronized across all selected plots.
                    </div>
                  </div>
                )}
              </div>
            </div>
          </aside>
        </div>
      )}

      {tileDrag && draggedCalculationCard && !showSelector && (
        <div
          ref={floatingCardRef}
          className="bands-multiview-floating-card"
          style={{
            width: `${tileDrag.width}px`,
            minWidth: `${tileDrag.width}px`,
          }}
          aria-hidden="true"
        >
          <article className="bands-multiview-card is-floating">
            <div className="bands-multiview-card-header">
              <div>
                <h3>{draggedCalculationCard.calc.project_name}</h3>
                <p>
                  {draggedCalculationCard.calc.folder_name?.trim() || "Root"}
                  {" · "}
                  {draggedCalculationCard.calc.cif_formula || draggedCalculationCard.calc.cif_filename}
                </p>
              </div>
              <div className="bands-multiview-card-header-side">
                <span>{formatDate(draggedCalculationCard.calc.completed_at)}</span>
                <button
                  type="button"
                  className="bands-multiview-reorder-handle"
                  disabled
                  tabIndex={-1}
                >
                  ⋮⋮
                </button>
              </div>
            </div>

            {(draggedCalculationCard.kPath || draggedCalculationCard.tags.length > 0) && (
              <div className="bands-multiview-card-meta">
                {draggedCalculationCard.kPath && (
                  <span className="calc-kpath">{draggedCalculationCard.kPath}</span>
                )}
                {draggedCalculationCard.tags.length > 0 && (
                  <div className="calc-tags">
                    {draggedCalculationCard.tags.map((tag) => (
                      <span key={tag.label} className={getCalcTagClass(tag)}>
                        {tag.label}
                      </span>
                    ))}
                  </div>
                )}
              </div>
            )}

            <div className="bands-multiview-card-controls">
              <label className="bands-multiview-projection-picker">
                <span>Projection:</span>
                <select value={draggedCalculationCard.projectionSelection} disabled tabIndex={-1}>
                  {draggedCalculationCard.projectionOptions.map((option) => (
                    <option key={option.value} value={option.value}>
                      {option.label}
                    </option>
                  ))}
                </select>
              </label>
            </div>

            <div
              className="bands-multiview-plot-frame"
              style={{ height: `${plotHeight}px`, minHeight: `${plotHeight}px` }}
            >
              <BandPlot
                data={draggedCalculationCard.calc.band_data}
                width={plotWidth}
                height={plotHeight}
                energyRange={resolvedEnergyRange}
                scfFermiEnergy={draggedCalculationCard.calc.scf_fermi_energy ?? undefined}
                viewerType="electronic"
                sharedSettings={sharedSettings}
                showSidebar={false}
                projectionSelection={draggedCalculationCard.projectionSelection}
                enableWheelRangeControl={false}
                enableHoverScrollLock={false}
              />
            </div>
          </article>
        </div>
      )}
    </div>
  );
}
