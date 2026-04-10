import { build } from "esbuild";
import { mkdtempSync, rmSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { spawnSync } from "node:child_process";

const tempDir = mkdtempSync(join(tmpdir(), "qcortado-transport-tests-"));
const bundledTestPath = join(tempDir, "transport-utils.test.mjs");

try {
  await build({
    entryPoints: ["tests/transport/transportUtils.test.ts"],
    outfile: bundledTestPath,
    bundle: true,
    format: "esm",
    platform: "node",
    target: ["node20"],
    sourcemap: "inline",
    logLevel: "silent",
  });

  const result = spawnSync(process.execPath, ["--test", bundledTestPath], {
    stdio: "inherit",
  });

  if (result.error) {
    throw result.error;
  }
  process.exit(result.status ?? 1);
} finally {
  rmSync(tempDir, { recursive: true, force: true });
}
