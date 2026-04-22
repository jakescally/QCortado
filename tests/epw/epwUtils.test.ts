import assert from "node:assert/strict";
import test from "node:test";
import {
  EpwMobilityDataset,
  EpwViewerData,
  buildEpwMobilitySeries,
  buildEpwTableSeries,
  collectEpwWarnings,
  getDefaultEpwTab,
  getEpwComponentOptions,
  getEpwMobilityValue,
  mergeEpwArtifacts,
} from "../../src/lib/epw";

const mobilityDataset: EpwMobilityDataset = {
  carrier_type: "electron",
  method: "IBTE",
  iteration: 2,
  component_labels: ["xx", "xy", "yy", "zz"],
  temperature_values_k: [300, 400],
  fermi_values_ev: [5.1, 5.2],
  density_values_cm3: [1e14, 2e14],
  population_values: [1e-12, 2e-12],
  mobility_values: [
    [10, 20],
    [1, 2],
    [30, 40],
    [50, 60],
  ],
};

test("default EPW tab prioritizes transport then superconductivity", () => {
  assert.equal(getDefaultEpwTab({ transport: { mobility: [mobilityDataset] } }), "transport");
  assert.equal(getDefaultEpwTab({ superconductivity: { lambda: 0.7 } }), "superconductivity");
  assert.equal(getDefaultEpwTab({ result_summary: { completed: true } }), "overview");
});

test("component options use transport-friendly ordering and avg fallback", () => {
  const options = getEpwComponentOptions(["xy", "zz", "xx", "yy"]);
  assert.deepEqual(options.map((option) => option.value), ["xx", "yy", "zz", "xy"]);
  assert.equal(getEpwMobilityValue(mobilityDataset, "avg", 0), 30);
});

test("mobility series extracts metric-specific values", () => {
  const mobility = buildEpwMobilitySeries(mobilityDataset, "mobility", "zz");
  assert.deepEqual(mobility.map((point) => point.x), [300, 400]);
  assert.deepEqual(mobility.map((point) => point.y), [50, 60]);

  const fermi = buildEpwMobilitySeries(mobilityDataset, "fermi", "xx");
  assert.deepEqual(fermi.map((point) => point.y), [5.1, 5.2]);

  const density = buildEpwMobilitySeries(mobilityDataset, "density", "xx");
  assert.deepEqual(density.map((point) => point.y), [1e14, 2e14]);
});

test("table series and warning aggregation are stable", () => {
  const tableSeries = buildEpwTableSeries({
    file_name: "sample.a2f",
    family: "a2f",
    title: "sample.a2f",
    column_labels: ["omega", "a2f"],
    rows: [[0, 1], [2, null], [3, 4]],
  });
  assert.deepEqual(tableSeries.map((point) => point.y), [1, null, 4]);

  const data: EpwViewerData = {
    result_summary: { warnings: ["EPW warning"] },
    transport: { warnings: ["EPW warning", "transport warning"] },
    errors: [{ code: "parse-partial", message: "partial" }],
  };
  assert.deepEqual(collectEpwWarnings(data), [
    "EPW warning",
    "transport warning",
    "[parse-partial] partial",
  ]);
});

test("artifact merge deduplicates saved and generated outputs", () => {
  const merged = mergeEpwArtifacts({
    result_summary: {
      generated_outputs: [
        { file_name: "epw.out", size_bytes: 10 },
        { file_name: "a2f.dat", size_bytes: 20 },
      ],
    },
    artifacts: [
      { file_name: "epw.out", size_bytes: 12 },
      { file_name: "z.dat", size_bytes: 30 },
    ],
  });
  assert.deepEqual(merged.map((entry) => entry.file_name), ["a2f.dat", "epw.out", "z.dat"]);
  assert.equal(merged.find((entry) => entry.file_name === "epw.out")?.size_bytes, 12);
});
