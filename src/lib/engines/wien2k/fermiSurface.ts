import type { SlurmResourceRequest } from "../../types";
import type { Wien2kBandsSpinChannel } from "./bands";

export interface Wien2kFermiSurfaceSettings {
  kMesh: [number, number, number];
  spinChannel: Wien2kBandsSpinChannel;
  spinOrbit: boolean;
  diagnosticLog: boolean;
}

export interface Wien2kSurfaceFileData {
  fileName: string;
  sizeBytes: number;
}

export interface Wien2kFermiSurfaceResult {
  calculationId: string;
  kGrid: [number, number, number];
  fermiEnergy?: number | null;
  primaryFile: string;
  bxsfFiles: Wien2kSurfaceFileData[];
  nativeOutput: string;
  diagnostics: string[];
}

export interface StartWien2kFermiSurfaceTaskParams {
  projectId: string;
  cifId: string;
  sourceScfCalculationId: string;
  settings: Wien2kFermiSurfaceSettings;
  resources?: SlurmResourceRequest | null;
}
