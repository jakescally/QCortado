/**
 * Shared engine identity and calculation kind types.
 *
 * These types describe platform-facing metadata only. They intentionally do not
 * define engine input schemas: QE inputs stay QE-specific, and future Wien2k
 * inputs should model Wien2k case workflows directly.
 */

export type ImplementedEngineId = "qe";

/**
 * Reserved for future remote-only Wien2k work. This is a type-level name only;
 * it does not indicate that a Wien2k backend is implemented.
 */
export type ReservedEngineId = "wien2k";

export type EngineId = ImplementedEngineId | ReservedEngineId;

export type EngineImplementationStatus = "implemented" | "configured" | "reserved";

export type EngineExecutionMode = "local" | "hpc" | "remote";

/**
 * Normalized calculation categories for project/result metadata.
 *
 * Existing engine-native payloads may use different calc_type strings during
 * migration. For example, current QE structure optimizations may be saved as
 * "optimization", "relax", or "vcrelax", while platform-facing metadata can
 * group them as "structure_optimization".
 */
export type CalculationKind =
  | "scf"
  | "structure_optimization"
  | "bands"
  | "dos"
  | "fermi_surface"
  | "phonon"
  | "hubbard_lrt"
  | "wannier"
  | "transport"
  | "epw"
  | "engine_setup"
  | "other";

export type CurrentQeCalculationType =
  | "scf"
  | "bands"
  | "dos"
  | "fermi_surface"
  | "hubbard_lrt"
  | "phonon"
  | "epw"
  | "wannier"
  | "transport"
  | "optimization"
  | "relax"
  | "vcrelax";

export interface EngineDescriptor {
  id: EngineId;
  label: string;
  status: EngineImplementationStatus;
  executionModes: readonly EngineExecutionMode[];
  calculationKinds: readonly CalculationKind[];
}

export interface EngineTagged {
  /**
   * Legacy saved calculations that do not contain this field should be treated
   * as QE until the storage migration adds explicit engine metadata.
   */
  engineId: EngineId;
}
