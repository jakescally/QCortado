import { getEngineWorkflowView } from "./registry";
import type { EngineWorkflowView } from "./plugin";
import type { CalculationKind, EngineId } from "./types";

export interface EngineWorkflowHostRoute {
  engineId: EngineId;
  kind: CalculationKind;
  view: EngineWorkflowView;
}

export const ENGINE_WORKFLOW_VIEWS: readonly EngineWorkflowView[] = [
  "scf-wizard",
  "bands-wizard",
  "dos-wizard",
  "fermi-surface-wizard",
  "hubbard-lrt-wizard",
  "phonon-wizard",
  "epw-wizard",
  "wannier-wizard",
  "transport-wizard",
];

export function isEngineWorkflowView(view: string): view is EngineWorkflowView {
  return ENGINE_WORKFLOW_VIEWS.includes(view as EngineWorkflowView);
}

export function resolveEngineWorkflowHostRoute(
  engineId: EngineId,
  kind: CalculationKind,
): EngineWorkflowHostRoute | null {
  const view = getEngineWorkflowView(engineId, kind);
  if (!view) return null;
  return {
    engineId,
    kind,
    view,
  };
}
