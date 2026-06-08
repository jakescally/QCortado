import assert from "node:assert/strict";
import test from "node:test";
import {
  GENERAL_HUBBARD_U_GUESS_EV,
  getHubbardEligibilityReason,
  getHubbardRecommendations,
  getLatestHubbardLrtValue,
  getScfHubbardUDisplayValues,
  getHundJDefaultEv,
  isDudarevDftUScf,
  resolveHubbardUDefault,
} from "../../src/lib/engines/qe/hubbard";

test("recommends Hubbard manifolds for transition metals and f elements", () => {
  const recommendations = getHubbardRecommendations(["O", "Fe", "Ce", "Si"]);
  assert.deepEqual(
    recommendations.map((entry) => [entry.element, entry.manifold]),
    [
      ["Fe", "3d"],
      ["Ce", "4f"],
    ],
  );
});

test("accepts only converged pure Dudarev DFT+U SCFs for LRT", () => {
  const base = {
    calc_type: "scf",
    result: { converged: true },
    parameters: {
      lda_plus_u: true,
      hubbard_formulation: 0,
      hubbard_parameters: [{ parameter: "U", manifold: "Fe-3d", value: 6 }],
    },
  };

  assert.equal(isDudarevDftUScf(base), true);
  assert.equal(getHubbardEligibilityReason(base), null);
  assert.equal(
    isDudarevDftUScf({
      ...base,
      parameters: { ...base.parameters, hubbard_formulation: 1 },
    }),
    false,
  );
  assert.match(
    getHubbardEligibilityReason({
      ...base,
      parameters: {
        ...base.parameters,
        hubbard_formulation: 1,
        hubbard_parameters: [
          { parameter: "U", manifold: "Fe-3d", value: 6 },
          { parameter: "J", manifold: "Fe-3d", value: 0.4 },
        ],
      },
    }) || "",
    /DFT\+U\+J/,
  );
});

test("uses newest LRT U before general 6 eV guess", () => {
  const calculations = [
    {
      id: "old",
      calc_type: "hubbard_lrt",
      completed_at: "2026-01-01T00:00:00Z",
      result: {
        converged: true,
        hubbard_lrt_data: {
          u_values: [{ element: "Fe", manifold: "3d", value_ev: 4.5 }],
        },
      },
    },
    {
      id: "new",
      calc_type: "hubbard_lrt",
      completed_at: "2026-02-01T00:00:00Z",
      result: {
        converged: true,
        hubbard_lrt_data: {
          u_values: [{ target: "Fe-3d", value_ev: 5.25 }],
        },
      },
    },
  ];

  assert.deepEqual(resolveHubbardUDefault("Fe", "3d", calculations), {
    value: 5.25,
    source: "lrt",
    label: "From Hubbard LRT (5.250 eV)",
    lrtCalcId: "new",
  });
  assert.equal(resolveHubbardUDefault("Ni", "3d", calculations).value, GENERAL_HUBBARD_U_GUESS_EV);
});

test("finds the latest matching Hubbard LRT value directly", () => {
  const calculations = [
    {
      id: "old",
      calc_type: "hubbard_lrt",
      completed_at: "2026-01-01T00:00:00Z",
      result: {
        converged: true,
        hubbard_lrt_data: {
          u_values: [{ target: "Fe-3d", value_ev: 4.5 }],
        },
      },
    },
    {
      id: "new",
      calc_type: "hubbard_lrt",
      completed_at: "2026-02-01T00:00:00Z",
      result: {
        converged: true,
        hubbard_lrt_data: {
          u_values: [{ target: "Fe-3d", value_ev: 5.25 }],
        },
      },
    },
  ];

  assert.deepEqual(getLatestHubbardLrtValue("Fe", "3d", calculations), {
    element: "Fe",
    manifold: "3d",
    value: 5.25,
    calcId: "new",
    completedAt: "2026-02-01T00:00:00Z",
  });
  assert.equal(getLatestHubbardLrtValue("Ni", "3d", calculations), null);
});

test("returns conservative Hund J defaults only for covered elements", () => {
  assert.equal(getHundJDefaultEv("Fe"), 0.437);
  assert.equal(getHundJDefaultEv("Ce"), null);
});

test("formats saved SCF Hubbard U values for source selectors", () => {
  assert.deepEqual(
    getScfHubbardUDisplayValues({
      calc_type: "scf",
      result: { converged: true },
      parameters: {
        lda_plus_u: true,
        hubbard_parameters: [
          { parameter: "U", manifold: "Ho-4f", value: 6.99 },
          { parameter: "J", manifold: "Ho-4f", value: 0.2 },
        ],
      },
    }),
    [{ target: "Ho-4f", value_ev: 6.99 }],
  );

  assert.deepEqual(
    getScfHubbardUDisplayValues({
      calc_type: "scf",
      result: { converged: true },
      parameters: {
        lda_plus_u: true,
        hubbard_u: { Fe: 5.25 },
        hubbard_manifold: { Fe: "3d" },
      },
    }),
    [{ target: "Fe-3d", value_ev: 5.25 }],
  );
});
