import { build } from "esbuild";
import { mkdtempSync, rmSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, basename } from "node:path";
import { spawnSync } from "node:child_process";

const testEntries = [
  "tests/bandDatasetAdapters/bandDatasetAdapters.test.ts",
  "tests/bandPlot/bandEnergyReference.test.ts",
  "tests/cifSubstitution/cifSubstitution.test.ts",
  "tests/enginePlugin/enginePluginRegistry.test.ts",
  "tests/epw/epwUtils.test.ts",
  "tests/transport/transportUtils.test.ts",
  "tests/qeBravais/qeBravaisInference.test.ts",
  "tests/kPathTransforms/sourceScfUsesPrimitiveCell.test.ts",
  "tests/brillouinZoneData/tetragonalBodyCentered.test.ts",
  "tests/hpcBandLimits/hpcBandLimits.test.ts",
  "tests/hpcProfileAdapter/hpcProfileAdapter.test.ts",
  "tests/hpcUtilization/hpcUtilization.test.ts",
  "tests/hubbard/hubbard.test.ts",
  "tests/magnetism/magnetism.test.ts",
  "tests/pseudopotentialCutoffs/ssspPresetCutoffs.test.ts",
];

const tempDir = mkdtempSync(join(tmpdir(), "qcortado-unit-tests-"));
const bundledTests = testEntries.map((entry) =>
  join(tempDir, `${basename(entry, ".ts")}.mjs`),
);

try {
  for (let i = 0; i < testEntries.length; i += 1) {
    await build({
      entryPoints: [testEntries[i]],
      outfile: bundledTests[i],
      bundle: true,
      format: "esm",
      platform: "node",
      target: ["node20"],
      sourcemap: "inline",
      logLevel: "silent",
    });
  }

  const result = spawnSync(process.execPath, ["--test", ...bundledTests], {
    stdio: "inherit",
  });

  if (result.error) {
    throw result.error;
  }
  process.exit(result.status ?? 1);
} finally {
  rmSync(tempDir, { recursive: true, force: true });
}
