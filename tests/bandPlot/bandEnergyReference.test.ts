import assert from "node:assert/strict";
import test from "node:test";
import {
  BandData,
  calculateDisplayedBandGap,
  calculateDisplayedBandGapFromSelectedValenceBand,
  resolveBandPlotEnergyReference,
} from "../../src/components/BandPlot";

function assertApprox(actual: number, expected: number, epsilon = 1e-9) {
  assert.ok(Math.abs(actual - expected) <= epsilon, `expected ${actual} to be within ${epsilon} of ${expected}`);
}

const semiconductorBands: BandData = {
  k_points: [0, 0.5, 1],
  energies: [
    [-1.2, -0.5, -1.0],
    [1.6, 0.7, 1.4],
  ],
  fermi_energy: 0.2,
  high_symmetry_points: [],
  n_bands: 2,
  n_kpoints: 3,
  band_gap: null,
  energy_range: [-1.2, 1.6],
  projections: null,
};

const metallicBands: BandData = {
  ...semiconductorBands,
  energies: [
    [-0.4, 0.3, -0.2],
    [0.8, 1.0, 0.9],
  ],
  energy_range: [-0.4, 1.0],
};

test("VBM zero mode shifts the plot zero while keeping the Fermi line offset", () => {
  const resolved = resolveBandPlotEnergyReference(
    semiconductorBands,
    0.2,
    "scf",
    "vbm",
    -0.48,
  );

  assert.equal(resolved.fermiEnergy, 0.2);
  assert.equal(resolved.zeroEnergyReferenceMode, "vbm");
  assert.equal(resolved.valenceBandMaximum, -0.48);
  assert.equal(resolved.zeroEnergy, -0.48);
});

test("VBM mode falls back to calculated Fermi until a valence band is selected", () => {
  const resolved = resolveBandPlotEnergyReference(
    semiconductorBands,
    0.2,
    "scf",
    "vbm",
  );

  assert.equal(resolved.zeroEnergyReferenceMode, "calculated-fermi");
  assert.equal(resolved.valenceBandMaximum, null);
  assertApprox(resolved.zeroEnergy, 0.2);
});

test("band gap display reanchors VBM to zero in VBM mode", () => {
  const gapAtFermi = calculateDisplayedBandGap(
    semiconductorBands.energies,
    semiconductorBands.k_points,
    0.2,
    0.2,
  );
  assert.ok(gapAtFermi);
  assertApprox(gapAtFermi.value, 1.2);
  assertApprox(gapAtFermi.vbm_energy, -0.7);
  assertApprox(gapAtFermi.cbm_energy, 0.5);

  const gapAtVbm = calculateDisplayedBandGap(
    semiconductorBands.energies,
    semiconductorBands.k_points,
    0.2,
    -0.5,
  );
  assert.ok(gapAtVbm);
  assertApprox(gapAtVbm.value, 1.2);
  assertApprox(gapAtVbm.vbm_energy, 0);
  assertApprox(gapAtVbm.cbm_energy, 1.2);
  assert.equal(gapAtVbm.is_direct, true);
});

test("selected valence band recomputes the displayed band gap from that band maximum", () => {
  const gapFromSelectedBand = calculateDisplayedBandGapFromSelectedValenceBand(
    semiconductorBands.energies,
    semiconductorBands.k_points,
    0,
    -0.5,
  );

  assert.ok(gapFromSelectedBand);
  assertApprox(gapFromSelectedBand.value, 1.2);
  assertApprox(gapFromSelectedBand.vbm_energy, 0);
  assertApprox(gapFromSelectedBand.cbm_energy, 1.2);
});

test("VBM zero mode uses the highest occupied state even when the system is metallic", () => {
  const resolved = resolveBandPlotEnergyReference(
    metallicBands,
    0.2,
    "scf",
    "vbm",
    -0.2,
  );

  assert.equal(resolved.zeroEnergyReferenceMode, "vbm");
  assertApprox(resolved.zeroEnergy, -0.2);
  assertApprox(resolved.valenceBandMaximum ?? Number.NaN, -0.2);
});

test("metallic systems still report no band gap", () => {
  assert.equal(
    calculateDisplayedBandGap(metallicBands.energies, metallicBands.k_points, 0.2, 0.2),
    null,
  );
});
