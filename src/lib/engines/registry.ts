import type { EngineDescriptor, EngineId, ImplementedEngineId } from "./types";

export const DEFAULT_ENGINE_ID: ImplementedEngineId = "qe";

export const QE_ENGINE_DESCRIPTOR: EngineDescriptor = {
  id: "qe",
  label: "Quantum ESPRESSO",
  status: "implemented",
  executionModes: ["local", "hpc"],
  calculationKinds: [
    "scf",
    "structure_optimization",
    "bands",
    "dos",
    "fermi_surface",
    "phonon",
    "hubbard_lrt",
    "wannier",
    "transport",
    "epw",
  ],
};

export const FALLBACK_ENGINE_DESCRIPTORS: readonly EngineDescriptor[] = [
  QE_ENGINE_DESCRIPTOR,
];

export function isImplementedEngineId(engineId: EngineId): engineId is ImplementedEngineId {
  return engineId === DEFAULT_ENGINE_ID;
}

export function getEngineLabel(
  descriptors: readonly EngineDescriptor[],
  engineId: EngineId,
): string {
  return descriptors.find((descriptor) => descriptor.id === engineId)?.label ?? engineId;
}

export function getEngineShortLabel(engine: EngineDescriptor): string {
  if (engine.id === "qe") {
    return "QE";
  }
  if (engine.id === "wien2k") {
    return "WIEN2k";
  }
  return engine.label;
}
