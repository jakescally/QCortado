import assert from "node:assert/strict";
import test from "node:test";
import {
  TransportDataset,
  buildTransportSeries,
  buildTransportTicks,
  createDefaultTransportAxisSettings,
  formatTransportNumber,
  getDatasetComponentValue,
  getNextTransportSliceSelection,
  getTransportComponentOptions,
  getTransportMetricDefinition,
  resolveTransportXAxisLabel,
  resolveTransportYAxisLabel,
} from "../../src/lib/transport";

const gridDataset: TransportDataset = {
  file_name: "grid.dat",
  component_labels: ["xx", "xy", "yy", "zz"],
  values: [
    [[1, 2, 3], [10, 20, 30]],
    [[4, 5, 6], [40, 50, 60]],
    [[7, 8, 9], [70, 80, 90]],
    [[11, 12, 13], [110, 120, 130]],
  ],
};

const transportResult = {
  engine: "boltzwann",
  seedname: "seed",
  source_wannier_calc_id: "calc-1",
  reference_fermi_energy_ev: 5,
  mu_values_ev: [4.9, 5.0, 5.1],
  mu_offsets_ev: [-0.1, 0, 0.1],
  temperature_values_k: [300, 600],
  relaxation_time_fs: 10,
  is_2d: false,
  conductivity: gridDataset,
  sigma_s: gridDataset,
  seebeck: gridDataset,
  kappa: gridDataset,
  tdf: null,
  notes: [],
  warnings: [],
  artifact_manifest: [],
};

test("component options include avg and preferred tensor ordering", () => {
  const options = getTransportComponentOptions(gridDataset);
  assert.deepEqual(
    options.map((option) => option.value),
    ["avg", "xx", "yy", "zz", "xy"],
  );
});

test("dataset component averaging uses available diagonal entries", () => {
  assert.equal(getDatasetComponentValue(gridDataset, "avg", 0, 1), (2 + 8 + 12) / 3);
  assert.equal(getDatasetComponentValue(gridDataset, "xy", 1, 2), 60);
});

test("transport series extraction supports mu and temperature sweeps", () => {
  const muSeries = buildTransportSeries({
    result: transportResult,
    dataset: gridDataset,
    component: "yy",
    axisMode: "mu",
    selectedTemperatureIndex: 1,
    selectedMuIndex: 0,
    showAbsoluteMu: false,
  });
  assert.deepEqual(muSeries.map((point) => point.x), [-0.1, 0, 0.1]);
  assert.deepEqual(muSeries.map((point) => point.y), [70, 80, 90]);

  const tempSeries = buildTransportSeries({
    result: transportResult,
    dataset: gridDataset,
    component: "zz",
    axisMode: "temperature",
    selectedTemperatureIndex: 0,
    selectedMuIndex: 2,
    showAbsoluteMu: false,
  });
  assert.deepEqual(tempSeries.map((point) => point.x), [300, 600]);
  assert.deepEqual(tempSeries.map((point) => point.y), [13, 130]);
});

test("heatmap slice selection updates only the active sweep axis", () => {
  assert.deepEqual(
    getNextTransportSliceSelection("mu", 1, 2, {
      selectedTemperatureIndex: 0,
      selectedMuIndex: 0,
    }),
    {
      selectedTemperatureIndex: 1,
      selectedMuIndex: 0,
    },
  );

  assert.deepEqual(
    getNextTransportSliceSelection("temperature", 1, 2, {
      selectedTemperatureIndex: 0,
      selectedMuIndex: 0,
    }),
    {
      selectedTemperatureIndex: 0,
      selectedMuIndex: 2,
    },
  );
});

test("tick generation and numeric formatting stay readable", () => {
  assert.deepEqual(buildTransportTicks(-0.3, 0.3, "normal"), [-0.2, 0, 0.2]);
  assert.equal(formatTransportNumber(12.3456, "manual", 2), "12.35");
  assert.equal(formatTransportNumber(1234.567, "auto", 4), "1235");
  assert.equal(formatTransportNumber(0.00045, "auto", 4), "4.500e-4");
});

test("metric metadata and default labels resolve as expected", () => {
  const settings = createDefaultTransportAxisSettings();
  const kappa = getTransportMetricDefinition("kappa");
  assert.equal(kappa.defaultUnit, "W/(m*K)");
  assert.equal(resolveTransportXAxisLabel("mu", settings), "Δμ (eV)");
  assert.equal(resolveTransportYAxisLabel(kappa, settings), "K (W/(m*K))");
});
