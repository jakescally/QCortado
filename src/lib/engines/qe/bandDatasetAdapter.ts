import type {
  BandData as LegacyBandData,
  BandProjectionData as LegacyBandProjectionData,
  BandProjectionGroup as LegacyBandProjectionGroup,
} from "../../../components/BandPlot";
import type { EngineId } from "../types";
import type {
  ArtifactReference,
  DatasetProvenance,
  DiagnosticHint,
  NumericSeries,
  PathMarker,
} from "../../viewers/types";
import type {
  BandDataset,
  BandProjectionDataset,
  BandProjectionGroup,
  BandProjectionGroupKind,
  NormalizedBandGap,
} from "../../viewers/bands/types";

export interface LegacyBandDataToBandDatasetOptions {
  engineId?: EngineId;
  calculationId?: string | null;
  projectId?: string | null;
  cifId?: string | null;
  sourceCalculationIds?: string[];
  generatedAt?: string | null;
  diagnostics?: DiagnosticHint[];
  artifacts?: ArtifactReference[];
  metadata?: Record<string, unknown>;
  referenceEnergyEv?: number | null;
}

export type QeBandDataToBandDatasetOptions = Omit<
  LegacyBandDataToBandDatasetOptions,
  "engineId"
>;

function buildBandDatasetProvenance(
  options: LegacyBandDataToBandDatasetOptions,
): DatasetProvenance {
  return {
    engineId: options.engineId ?? "qe",
    calculationKind: "bands",
    calculationId: options.calculationId ?? null,
    projectId: options.projectId ?? null,
    cifId: options.cifId ?? null,
    sourceCalculationIds: options.sourceCalculationIds,
    generatedAt: options.generatedAt ?? null,
  };
}

function mapBandSeries(energies: number[][]): NumericSeries[] {
  return energies.map((values, bandIndex) => ({
    id: `band-${bandIndex + 1}`,
    label: `Band ${bandIndex + 1}`,
    values,
    unit: "eV",
    metadata: {
      legacyBandIndex: bandIndex,
    },
  }));
}

function mapPathMarkers(markers: LegacyBandData["high_symmetry_points"]): PathMarker[] {
  return markers.map((marker) => ({
    x: marker.k_distance,
    label: marker.label,
  }));
}

function mapBandGap(gap: LegacyBandData["band_gap"]): NormalizedBandGap | null {
  if (!gap) {
    return null;
  }

  return {
    valueEv: gap.value,
    isDirect: gap.is_direct,
    vbmX: gap.vbm_k,
    cbmX: gap.cbm_k,
    vbmEnergyEv: gap.vbm_energy,
    cbmEnergyEv: gap.cbm_energy,
  };
}

function normalizeProjectionKind(
  kind: LegacyBandProjectionGroup["kind"],
  fallbackKind: BandProjectionGroupKind,
): BandProjectionGroupKind {
  if (
    kind === "atom" ||
    kind === "orbital" ||
    kind === "element_orbital" ||
    kind === "other"
  ) {
    return kind;
  }

  return fallbackKind;
}

function mapProjectionGroup(
  group: LegacyBandProjectionGroup,
  fallbackKind: BandProjectionGroupKind,
): BandProjectionGroup {
  const kind = normalizeProjectionKind(group.kind, fallbackKind);
  const metadata = kind === group.kind ? undefined : { legacyKind: group.kind };

  return {
    id: group.id,
    label: group.label,
    kind,
    weights: group.weights,
    ...(metadata ? { metadata } : {}),
  };
}

function mapProjectionDataset(
  projections: LegacyBandProjectionData | null | undefined,
): BandProjectionDataset | null {
  if (!projections) {
    return null;
  }

  return {
    source: projections.source,
    groups: [
      ...projections.atom_groups.map((group) => mapProjectionGroup(group, "atom")),
      ...projections.orbital_groups.map((group) => mapProjectionGroup(group, "orbital")),
      ...(projections.element_orbital_groups ?? []).map((group) =>
        mapProjectionGroup(group, "element_orbital"),
      ),
    ],
    normalization: "raw",
  };
}

/**
 * Converts the current BandPlot/BandStructureWizard payload into the
 * engine-neutral viewer dataset without changing existing rendering paths.
 */
export function legacyBandDataToBandDataset(
  data: LegacyBandData,
  options: LegacyBandDataToBandDatasetOptions = {},
): BandDataset {
  const referenceEnergyEv =
    options.referenceEnergyEv === undefined ? data.fermi_energy : options.referenceEnergyEv;

  return {
    schema: "cortado.band_path.v1",
    provenance: buildBandDatasetProvenance(options),
    quantity: "electronic_energy",
    x: data.k_points,
    series: mapBandSeries(data.energies),
    xUnit: "path_distance",
    yUnit: "eV",
    referenceEnergyEv,
    markers: mapPathMarkers(data.high_symmetry_points),
    bandGap: mapBandGap(data.band_gap),
    projections: mapProjectionDataset(data.projections),
    diagnostics: options.diagnostics,
    artifacts: options.artifacts,
    metadata: {
      nBands: data.n_bands,
      nKpoints: data.n_kpoints,
      energyRangeEv: data.energy_range,
      legacyFermiEnergyEv: data.fermi_energy,
      ...options.metadata,
    },
  };
}

export function qeBandDataToBandDataset(
  data: LegacyBandData,
  options: QeBandDataToBandDatasetOptions = {},
): BandDataset {
  return legacyBandDataToBandDataset(data, {
    ...options,
    engineId: "qe",
  });
}
