import {
  BandStructureWizard,
  ElectronicDOSWizard,
  EpwWizard,
  FermiSurfaceWizard,
  HubbardLrtWizard,
  PhononWizard,
  SCFWizard,
  TransportWizard,
  WannierWizard,
} from "../qe";
import type { BandData } from "../BandPlot";
import type { ElectronicDOSData } from "../ElectronicDOSPlot";
import type {
  CalculationRun,
  SavedBandsCalculationContext,
  WannierBandOverlayOption,
} from "../ProjectDashboard";
import type { TransportResult } from "../../lib/transport";
import type {
  CrystalData,
  ExecutionMode,
  HpcProfile,
  OptimizedStructureOption,
  QeSmearingType,
  SCFPreset,
} from "../../lib/types";
import type { EngineWorkflowHostRoute } from "../../lib/engines";
import type { EngineWorkflowView } from "../../lib/engines/plugin";
import type { EpwViewerPayload } from "../qe";
import { Wien2kStructureWizard } from "../wien2k/Wien2kStructureWizard";
import { Wien2kScfWizard } from "../wien2k/Wien2kScfWizard";
import { Wien2kBandsWizard } from "../wien2k/Wien2kBandsWizard";

export interface ScfWorkflowContext {
  cifId: string;
  crystalData: CrystalData;
  cifContent: string;
  filename: string;
  projectId: string;
  initialPreset?: SCFPreset;
  presetLock?: boolean;
  optimizedStructures?: OptimizedStructureOption[];
  calculations?: CalculationRun[];
}

export interface SourceScfWorkflowContext {
  cifId: string;
  crystalData: CrystalData;
  projectId: string;
  scfCalculations: CalculationRun[];
}

export interface TransportWorkflowContext {
  cifId: string;
  crystalData: CrystalData;
  projectId: string;
  wannierCalculations: CalculationRun[];
}

export interface EpwWorkflowContext {
  cifId: string;
  crystalData: CrystalData;
  projectId: string;
  calculations: CalculationRun[];
}

export interface Wien2kStructureWorkflowContext {
  cifId: string;
  crystalData: CrystalData;
  projectId: string;
}

export interface EngineWorkflowHostContexts {
  scf: ScfWorkflowContext | null;
  bands: SourceScfWorkflowContext | null;
  dos: SourceScfWorkflowContext | null;
  wannier: SourceScfWorkflowContext | null;
  fermiSurface: SourceScfWorkflowContext | null;
  hubbardLrt: SourceScfWorkflowContext | null;
  phonons: SourceScfWorkflowContext | null;
  transport: TransportWorkflowContext | null;
  epw: EpwWorkflowContext | null;
  structureSetup: Wien2kStructureWorkflowContext | null;
  wien2kScf: ScfWorkflowContext | null;
}

export interface PhononWorkflowViewerData {
  dos_data: unknown | null;
  dispersion_data: unknown | null;
}

export type EngineWorkflowBackDestination = "project-browser" | "project-dashboard";

interface EngineWorkflowHostRuntime {
  qePath: string | null;
  defaultSmearing: QeSmearingType;
  executionMode: ExecutionMode;
  onExecutionModeChange: (mode: ExecutionMode) => Promise<void> | void;
  activeHpcProfile: HpcProfile | null;
}

export interface EngineWorkflowHostProps {
  route: EngineWorkflowHostRoute;
  runtime: EngineWorkflowHostRuntime;
  contexts: EngineWorkflowHostContexts;
  reconnectTaskId: string | null;
  onBack: (view: EngineWorkflowView, destination: EngineWorkflowBackDestination) => void;
  onViewBands: (
    bandData: BandData,
    fermiEnergy: number | null,
    calculationParameters?: Record<string, unknown> | null,
    calculationContext?: SavedBandsCalculationContext | null,
  ) => void;
  onViewDos: (dosData: ElectronicDOSData, fermiEnergy: number | null) => void;
  onViewWannier: (
    result: unknown,
    fermiEnergy: number | null,
    overlayOptions?: WannierBandOverlayOption[],
  ) => void;
  onViewTransport: (transportData: TransportResult) => void;
  onViewPhonons: (phononData: PhononWorkflowViewerData, viewMode: "bands" | "dos") => void;
  onViewEpw: (epwData: EpwViewerPayload["data"], rawOutput?: string | null) => void;
  onStructureSourceSaved: () => void;
  onWien2kScfSaved: () => void;
}

const EMPTY_CRYSTAL_DATA: CrystalData = {
  cell_length_a: { value: 0 },
  cell_length_b: { value: 0 },
  cell_length_c: { value: 0 },
  cell_angle_alpha: { value: 0 },
  cell_angle_beta: { value: 0 },
  cell_angle_gamma: { value: 0 },
  atom_sites: [],
  symmetry_operations: [],
  anisotropic_params: [],
};

export function canRenderEngineWorkflowHost({
  route,
  runtime,
  contexts,
  reconnectTaskId,
}: Pick<EngineWorkflowHostProps, "route" | "runtime" | "contexts" | "reconnectTaskId">): boolean {
  if (route.engineId === "wien2k") {
    if (route.view === "wien2k-structure-wizard") return Boolean(contexts.structureSetup);
    if (route.view === "wien2k-scf-wizard") return Boolean(contexts.wien2kScf);
    if (route.view === "bands-wizard") return Boolean(contexts.bands);
    return false;
  }
  if (route.engineId !== "qe") {
    return false;
  }

  const hasRuntime = Boolean(runtime.qePath) || runtime.executionMode === "hpc";
  if (!hasRuntime) {
    return false;
  }

  switch (route.view) {
    case "scf-wizard":
      return true;
    case "bands-wizard":
      return Boolean(contexts.bands || reconnectTaskId);
    case "dos-wizard":
      return Boolean(contexts.dos || reconnectTaskId);
    case "wannier-wizard":
      return Boolean(contexts.wannier || reconnectTaskId);
    case "transport-wizard":
      return Boolean(contexts.transport || reconnectTaskId);
    case "fermi-surface-wizard":
      return Boolean(contexts.fermiSurface || reconnectTaskId);
    case "hubbard-lrt-wizard":
      return Boolean(contexts.hubbardLrt || reconnectTaskId);
    case "phonon-wizard":
      return Boolean(contexts.phonons || reconnectTaskId);
    case "epw-wizard":
      return Boolean(contexts.epw || reconnectTaskId);
    default:
      return false;
  }
}

export function EngineWorkflowHost(props: EngineWorkflowHostProps) {
  if (!canRenderEngineWorkflowHost(props)) {
    return null;
  }

  const { route, runtime, contexts, reconnectTaskId, onBack } = props;
  const qePath = runtime.qePath || "";
  const reconnectTaskIdValue = reconnectTaskId ?? undefined;

  switch (route.view) {
    case "wien2k-structure-wizard": {
      const context = contexts.structureSetup;
      if (!context) return null;
      return (
        <Wien2kStructureWizard
          projectId={context.projectId}
          cifId={context.cifId}
          crystalData={context.crystalData}
          onBack={() => onBack(route.view, "project-dashboard")}
          onSaved={props.onStructureSourceSaved}
        />
      );
    }
    case "wien2k-scf-wizard": {
      const context = contexts.wien2kScf;
      if (!context) return null;
      return (
        <Wien2kScfWizard
          projectId={context.projectId}
          cifId={context.cifId}
          calculations={context.calculations ?? []}
          activeHpcProfile={runtime.activeHpcProfile}
          onBack={() => onBack(route.view, "project-dashboard")}
          onSaved={props.onWien2kScfSaved}
        />
      );
    }
    case "bands-wizard": {
      if (route.engineId === "wien2k") {
        const context = contexts.bands;
        if (!context) return null;
        return (
          <Wien2kBandsWizard
            projectId={context.projectId}
            cifId={context.cifId}
            crystalData={context.crystalData}
            scfCalculations={context.scfCalculations}
            activeHpcProfile={runtime.activeHpcProfile}
            onBack={() => onBack(route.view, "project-dashboard")}
            onViewBands={props.onViewBands}
          />
        );
      }
      const context = contexts.bands;
      return (
        <BandStructureWizard
          qePath={qePath}
          defaultSmearing={runtime.defaultSmearing}
          executionMode={runtime.executionMode}
          onExecutionModeChange={runtime.onExecutionModeChange}
          activeHpcProfile={runtime.activeHpcProfile}
          onViewBands={props.onViewBands}
          onBack={() => onBack(route.view, "project-dashboard")}
          projectId={context?.projectId ?? ""}
          cifId={context?.cifId ?? ""}
          crystalData={context?.crystalData ?? EMPTY_CRYSTAL_DATA}
          scfCalculations={context?.scfCalculations ?? []}
          reconnectTaskId={reconnectTaskIdValue}
        />
      );
    }
    case "dos-wizard": {
      const context = contexts.dos;
      return (
        <ElectronicDOSWizard
          qePath={qePath}
          defaultSmearing={runtime.defaultSmearing}
          executionMode={runtime.executionMode}
          onExecutionModeChange={runtime.onExecutionModeChange}
          activeHpcProfile={runtime.activeHpcProfile}
          onViewDos={props.onViewDos}
          onBack={() => onBack(route.view, "project-dashboard")}
          projectId={context?.projectId ?? ""}
          cifId={context?.cifId ?? ""}
          crystalData={context?.crystalData ?? EMPTY_CRYSTAL_DATA}
          scfCalculations={context?.scfCalculations ?? []}
          reconnectTaskId={reconnectTaskIdValue}
        />
      );
    }
    case "wannier-wizard": {
      const context = contexts.wannier;
      return (
        <WannierWizard
          qePath={qePath}
          defaultSmearing={runtime.defaultSmearing}
          executionMode={runtime.executionMode}
          onExecutionModeChange={runtime.onExecutionModeChange}
          activeHpcProfile={runtime.activeHpcProfile}
          onViewWannier={props.onViewWannier}
          onBack={() => onBack(route.view, "project-dashboard")}
          projectId={context?.projectId ?? ""}
          cifId={context?.cifId ?? ""}
          crystalData={context?.crystalData ?? EMPTY_CRYSTAL_DATA}
          scfCalculations={context?.scfCalculations ?? []}
          reconnectTaskId={reconnectTaskIdValue}
        />
      );
    }
    case "transport-wizard": {
      const context = contexts.transport;
      return (
        <TransportWizard
          qePath={qePath}
          executionMode={runtime.executionMode}
          onExecutionModeChange={runtime.onExecutionModeChange}
          activeHpcProfile={runtime.activeHpcProfile}
          onViewTransport={props.onViewTransport}
          onBack={() => onBack(route.view, "project-dashboard")}
          projectId={context?.projectId ?? ""}
          cifId={context?.cifId ?? ""}
          crystalData={context?.crystalData ?? EMPTY_CRYSTAL_DATA}
          wannierCalculations={context?.wannierCalculations ?? []}
          reconnectTaskId={reconnectTaskIdValue}
        />
      );
    }
    case "fermi-surface-wizard": {
      const context = contexts.fermiSurface;
      return (
        <FermiSurfaceWizard
          qePath={qePath}
          defaultSmearing={runtime.defaultSmearing}
          executionMode={runtime.executionMode}
          onExecutionModeChange={runtime.onExecutionModeChange}
          activeHpcProfile={runtime.activeHpcProfile}
          onBack={() => onBack(route.view, "project-dashboard")}
          projectId={context?.projectId ?? ""}
          cifId={context?.cifId ?? ""}
          crystalData={context?.crystalData ?? EMPTY_CRYSTAL_DATA}
          scfCalculations={context?.scfCalculations ?? []}
          reconnectTaskId={reconnectTaskIdValue}
        />
      );
    }
    case "hubbard-lrt-wizard": {
      const context = contexts.hubbardLrt;
      return (
        <HubbardLrtWizard
          qePath={qePath}
          executionMode={runtime.executionMode}
          onExecutionModeChange={runtime.onExecutionModeChange}
          activeHpcProfile={runtime.activeHpcProfile}
          onBack={() => onBack(route.view, "project-dashboard")}
          projectId={context?.projectId ?? ""}
          cifId={context?.cifId ?? ""}
          crystalData={context?.crystalData ?? EMPTY_CRYSTAL_DATA}
          scfCalculations={context?.scfCalculations ?? []}
          reconnectTaskId={reconnectTaskIdValue}
        />
      );
    }
    case "phonon-wizard": {
      const context = contexts.phonons;
      return (
        <PhononWizard
          qePath={qePath}
          executionMode={runtime.executionMode}
          onExecutionModeChange={runtime.onExecutionModeChange}
          activeHpcProfile={runtime.activeHpcProfile}
          onViewPhonons={props.onViewPhonons}
          onBack={() => onBack(route.view, "project-dashboard")}
          projectId={context?.projectId ?? ""}
          cifId={context?.cifId ?? ""}
          crystalData={context?.crystalData ?? EMPTY_CRYSTAL_DATA}
          scfCalculations={context?.scfCalculations ?? []}
          reconnectTaskId={reconnectTaskIdValue}
        />
      );
    }
    case "epw-wizard": {
      const context = contexts.epw;
      return (
        <EpwWizard
          qePath={qePath}
          executionMode={runtime.executionMode}
          activeHpcProfile={runtime.activeHpcProfile}
          onViewEPW={props.onViewEpw}
          onBack={() => onBack(route.view, "project-dashboard")}
          projectId={context?.projectId ?? ""}
          cifId={context?.cifId ?? ""}
          crystalData={context?.crystalData ?? EMPTY_CRYSTAL_DATA}
          calculations={context?.calculations ?? []}
          reconnectTaskId={reconnectTaskIdValue}
        />
      );
    }
    case "scf-wizard": {
      const context = contexts.scf;
      return (
        <SCFWizard
          qePath={qePath}
          defaultSmearing={runtime.defaultSmearing}
          executionMode={runtime.executionMode}
          activeHpcProfile={runtime.activeHpcProfile}
          onBack={() => onBack(route.view, context ? "project-dashboard" : "project-browser")}
          initialCif={context ?? undefined}
          initialPreset={context?.initialPreset}
          presetLock={context?.presetLock}
          optimizedStructures={context?.optimizedStructures}
          reconnectTaskId={reconnectTaskIdValue}
        />
      );
    }
  }
}
