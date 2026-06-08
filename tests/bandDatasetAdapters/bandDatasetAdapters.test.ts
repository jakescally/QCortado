import assert from "node:assert/strict";
import test from "node:test";
import { bandDatasetToBandData } from "../../src/components/BandPlot";
import type { BandData } from "../../src/components/BandPlot";
import {
  legacyBandDataToBandDataset,
  qeBandDataToBandDataset,
} from "../../src/lib/engines/qe";

const legacyBandData: BandData = {
  k_points: [0, 0.5, 1],
  energies: [
    [-1.2, -0.5, -1.0],
    [1.6, 0.7, 1.4],
  ],
  fermi_energy: 0.2,
  high_symmetry_points: [
    { k_distance: 0, label: "G" },
    { k_distance: 1, label: "X" },
  ],
  n_bands: 2,
  n_kpoints: 3,
  band_gap: {
    value: 1.2,
    is_direct: true,
    vbm_k: 0.5,
    cbm_k: 0.5,
    vbm_energy: -0.5,
    cbm_energy: 0.7,
  },
  energy_range: [-1.2, 1.6],
  projections: {
    source: "projwfc.x",
    atom_groups: [
      {
        id: "element-si",
        label: "Si total",
        kind: "atom",
        weights: [
          [0.8, 0.7, 0.6],
          [0.2, 0.3, 0.4],
        ],
      },
    ],
    orbital_groups: [
      {
        id: "orbital-p",
        label: "P orbitals",
        kind: "orbital",
        weights: [
          [0.1, 0.2, 0.3],
          [0.4, 0.5, 0.6],
        ],
      },
    ],
    element_orbital_groups: [
      {
        id: "element-orbital-si-p",
        label: "Si P orbitals",
        kind: "element_orbital",
        weights: [
          [0.05, 0.06, 0.07],
          [0.08, 0.09, 0.1],
        ],
      },
    ],
  },
};

test("legacyBandDataToBandDataset maps the current BandData shape to normalized band data", () => {
  const dataset = legacyBandDataToBandDataset(legacyBandData, {
    calculationId: "bands-1",
    projectId: "project-1",
    cifId: "cif-1",
    generatedAt: "2026-05-21T00:00:00.000Z",
    metadata: {
      sourceFormat: "legacy-band-data",
    },
  });

  assert.equal(dataset.schema, "cortado.band_path.v1");
  assert.equal(dataset.provenance.engineId, "qe");
  assert.equal(dataset.provenance.calculationKind, "bands");
  assert.equal(dataset.provenance.calculationId, "bands-1");
  assert.equal(dataset.provenance.projectId, "project-1");
  assert.equal(dataset.provenance.cifId, "cif-1");
  assert.equal(dataset.provenance.generatedAt, "2026-05-21T00:00:00.000Z");
  assert.equal(dataset.quantity, "electronic_energy");
  assert.equal(dataset.xUnit, "path_distance");
  assert.equal(dataset.yUnit, "eV");
  assert.deepEqual(dataset.x, legacyBandData.k_points);
  assert.equal(dataset.referenceEnergyEv, 0.2);
  assert.deepEqual(dataset.markers, [
    { x: 0, label: "G" },
    { x: 1, label: "X" },
  ]);
  assert.deepEqual(dataset.series, [
    {
      id: "band-1",
      label: "Band 1",
      values: [-1.2, -0.5, -1.0],
      unit: "eV",
      metadata: { legacyBandIndex: 0 },
    },
    {
      id: "band-2",
      label: "Band 2",
      values: [1.6, 0.7, 1.4],
      unit: "eV",
      metadata: { legacyBandIndex: 1 },
    },
  ]);
  assert.deepEqual(dataset.bandGap, {
    valueEv: 1.2,
    isDirect: true,
    vbmX: 0.5,
    cbmX: 0.5,
    vbmEnergyEv: -0.5,
    cbmEnergyEv: 0.7,
  });
  assert.equal(dataset.projections?.source, "projwfc.x");
  assert.equal(dataset.projections?.normalization, "raw");
  assert.deepEqual(
    dataset.projections?.groups.map((group) => group.kind),
    ["atom", "orbital", "element_orbital"],
  );
  assert.deepEqual(dataset.metadata, {
    nBands: 2,
    nKpoints: 3,
    energyRangeEv: [-1.2, 1.6],
    legacyFermiEnergyEv: 0.2,
    sourceFormat: "legacy-band-data",
  });
});

test("qeBandDataToBandDataset pins the adapter provenance to the QE engine", () => {
  const dataset = qeBandDataToBandDataset(legacyBandData, {
    referenceEnergyEv: null,
  });

  assert.equal(dataset.provenance.engineId, "qe");
  assert.equal(dataset.referenceEnergyEv, null);
});

test("bandDatasetToBandData preserves legacy BandPlot fields for existing rendering", () => {
  const dataset = legacyBandDataToBandDataset(legacyBandData);
  const restored = bandDatasetToBandData(dataset);

  assert.deepEqual(restored.k_points, legacyBandData.k_points);
  assert.deepEqual(restored.energies, legacyBandData.energies);
  assert.equal(restored.fermi_energy, legacyBandData.fermi_energy);
  assert.deepEqual(restored.high_symmetry_points, legacyBandData.high_symmetry_points);
  assert.equal(restored.n_bands, legacyBandData.n_bands);
  assert.equal(restored.n_kpoints, legacyBandData.n_kpoints);
  assert.deepEqual(restored.band_gap, legacyBandData.band_gap);
  assert.deepEqual(restored.energy_range, legacyBandData.energy_range);
  assert.deepEqual(
    restored.projections?.element_orbital_groups,
    legacyBandData.projections?.element_orbital_groups,
  );
});
