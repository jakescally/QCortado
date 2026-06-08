import assert from "node:assert/strict";
import test from "node:test";
import { summarizeSelectedPseudoCutoffs } from "../../src/lib/engines/qe/pseudopotentialCutoffs";

test("uses SSSP JSON cutoffs even before UPF metadata has loaded", () => {
  const summary = summarizeSelectedPseudoCutoffs(
    ["O", "Fe"],
    {
      O: "O.pbe-n-kjpaw_psl.0.1.UPF",
      Fe: "Fe.pbe-spn-kjpaw_psl.0.2.1.UPF",
    },
    {},
    {
      O: {
        filename: "O.pbe-n-kjpaw_psl.0.1.UPF",
        cutoff_wfc: 75,
        cutoff_rho: 600,
      },
      Fe: {
        filename: "Fe.pbe-spn-kjpaw_psl.0.2.1.UPF",
        cutoff_wfc: 90,
        cutoff_rho: 1080,
      },
    },
  );

  assert.equal(summary.maxWfc, 90);
  assert.equal(summary.maxRho, 1080);
  assert.equal(summary.wfcStatus, "parsed");
  assert.equal(summary.rhoStatus, "parsed");
  assert.equal(summary.wfcProvenance, "sssp");
  assert.equal(summary.rhoProvenance, "sssp");
  assert.equal(summary.hasMissingCutoff, false);
});

test("chooses the strictest wavefunction and charge cutoffs independently across selected atoms", () => {
  const summary = summarizeSelectedPseudoCutoffs(
    ["A", "B"],
    {
      A: "A.upf",
      B: "B.upf",
    },
    {},
    {
      A: {
        filename: "A.upf",
        cutoff_wfc: 100,
        cutoff_rho: 500,
      },
      B: {
        filename: "B.upf",
        cutoff_wfc: 80,
        cutoff_rho: 800,
      },
    },
  );

  assert.equal(summary.maxWfc, 100);
  assert.equal(summary.maxRho, 800);
  assert.equal(summary.wfcProvenance, "sssp");
  assert.equal(summary.rhoProvenance, "sssp");
});
