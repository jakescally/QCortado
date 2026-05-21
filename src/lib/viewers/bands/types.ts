import type { NumericSeries, PathMarker, ResultDatasetBase } from "../types";

export type BandDatasetSchema = "cortado.band_path.v1";
export type BandQuantity = "electronic_energy";
export type BandPathXAxisUnit = "path_distance";
export type BandEnergyUnit = "eV";

export interface NormalizedBandGap {
  valueEv: number;
  isDirect: boolean;
  vbmX?: number | null;
  cbmX?: number | null;
  vbmEnergyEv?: number | null;
  cbmEnergyEv?: number | null;
}

export type BandProjectionGroupKind = "atom" | "orbital" | "element_orbital" | "other";

export interface BandProjectionGroup {
  id: string;
  label: string;
  kind: BandProjectionGroupKind;
  weights: number[][];
  metadata?: Record<string, unknown>;
}

export interface BandProjectionDataset {
  source: string;
  groups: BandProjectionGroup[];
  normalization?: "raw" | "global" | "band" | null;
  metadata?: Record<string, unknown>;
}

/**
 * Engine-neutral electronic band dataset for shared band viewers.
 *
 * This is an output/viewer contract only. QE and future engines should adapt
 * their native parser results into this shape without sharing input configs.
 */
export interface BandPathDataset extends ResultDatasetBase<BandDatasetSchema> {
  quantity: BandQuantity;
  x: number[];
  series: NumericSeries[];
  xUnit: BandPathXAxisUnit;
  yUnit: BandEnergyUnit;
  referenceEnergyEv?: number | null;
  markers: PathMarker[];
  bandGap?: NormalizedBandGap | null;
  projections?: BandProjectionDataset | null;
}

export type BandDataset = BandPathDataset;
