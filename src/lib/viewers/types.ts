import type { CalculationKind, EngineId } from "../engines/types";

export type DatasetSchema =
  | "cortado.scf_summary.v1"
  | "cortado.band_path.v1"
  | "cortado.dos.v1"
  | "cortado.line_path.v1"
  | "cortado.phonon_dos.v1"
  | "cortado.transport_tensor.v1"
  | "cortado.numeric_table.v1"
  | "cortado.series.v1";

export type EnergyUnit = "eV" | "Ry" | "Ha";
export type ForceUnit = "Ry/Bohr" | "eV/angstrom";
export type StressUnit = "kbar" | "GPa";
export type FrequencyUnit = "cm-1" | "THz";

export type Vector3 = [number, number, number];
export type Matrix3x3 = [Vector3, Vector3, Vector3];

export interface ScalarQuantity<Unit extends string = string> {
  value: number;
  unit: Unit;
}

export type DiagnosticSeverity = "info" | "warning" | "error";

export type DiagnosticSource =
  | "platform"
  | "engine"
  | "parser"
  | "viewer"
  | "hpc"
  | "storage";

export interface DiagnosticHint {
  id?: string;
  severity: DiagnosticSeverity;
  source: DiagnosticSource;
  code?: string;
  title: string;
  message: string;
  suggestedAction?: string;
  relatedArtifactIds?: string[];
  metadata?: Record<string, unknown>;
}

export type ArtifactLocationKind =
  | "local_path"
  | "remote_path"
  | "project_relative_path"
  | "archive_relative_path";

export type ArtifactRole =
  | "input"
  | "output"
  | "log"
  | "scratch"
  | "viewer"
  | "provenance"
  | "unknown";

export interface ArtifactReference {
  id: string;
  label?: string;
  role: ArtifactRole;
  locationKind: ArtifactLocationKind;
  path: string;
  engineId?: EngineId;
  sizeBytes?: number | null;
  mediaType?: string | null;
  checksum?: string | null;
  metadata?: Record<string, unknown>;
}

export interface DatasetProvenance {
  engineId: EngineId;
  calculationKind: CalculationKind;
  calculationId?: string | null;
  projectId?: string | null;
  cifId?: string | null;
  sourceCalculationIds?: string[];
  generatedAt?: string | null;
}

export interface ResultDatasetBase<Schema extends DatasetSchema = DatasetSchema> {
  schema: Schema;
  provenance: DatasetProvenance;
  diagnostics?: DiagnosticHint[];
  artifacts?: ArtifactReference[];
  metadata?: Record<string, unknown>;
}

export type ScfConvergenceState =
  | "converged"
  | "not_converged"
  | "failed"
  | "cancelled"
  | "unknown";

export interface NormalizedScfSummary
  extends ResultDatasetBase<"cortado.scf_summary.v1"> {
  convergence: ScfConvergenceState;
  totalEnergy?: ScalarQuantity<EnergyUnit> | null;
  fermiEnergyEv?: number | null;
  scfSteps?: number | null;
  wallTimeSeconds?: number | null;
  totalMagnetization?: number | null;
  atomicMagneticMoments?: Vector3[] | null;
  forces?: {
    values: Vector3[];
    unit: ForceUnit;
  } | null;
  stress?: {
    values: Matrix3x3;
    unit: StressUnit;
  } | null;
}

export interface NumericSeries {
  id: string;
  label: string;
  values: Array<number | null>;
  unit?: string | null;
  metadata?: Record<string, unknown>;
}

export interface PathMarker {
  x: number;
  label: string;
}

export interface NumericTableDataset
  extends ResultDatasetBase<"cortado.numeric_table.v1"> {
  title: string;
  columns: string[];
  rows: Array<Array<number | string | null>>;
}

export interface SeriesDataset extends ResultDatasetBase<"cortado.series.v1"> {
  x: number[];
  y: Array<number | null>;
  xLabel: string;
  yLabel: string;
  xUnit?: string | null;
  yUnit?: string | null;
}

