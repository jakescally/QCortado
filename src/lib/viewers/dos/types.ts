import type { NumericSeries, ResultDatasetBase } from "../types";

export type DosDatasetSchema = "cortado.dos.v1";
export type DosQuantity = "electronic_dos";
export type DosEnergyUnit = "eV";
export type DosSpinChannel = "total" | "up" | "down" | "noncollinear" | "unknown";

export interface DosChannel extends NumericSeries {
  spin?: DosSpinChannel | null;
}

/**
 * Engine-neutral electronic DOS dataset for shared DOS viewers.
 *
 * The dataset describes plotted output only. Engine-specific DOS inputs and
 * file parsing stay in engine modules.
 */
export interface DosDataset extends ResultDatasetBase<DosDatasetSchema> {
  quantity: DosQuantity;
  x: number[];
  y: number[];
  xUnit: DosEnergyUnit;
  yUnit: string;
  referenceEnergyEv?: number | null;
  channels?: DosChannel[];
  energyRangeEv?: [number, number] | null;
  maxDos?: number | null;
}

