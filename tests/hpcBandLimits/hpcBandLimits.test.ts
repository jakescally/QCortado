import assert from "node:assert/strict";
import test from "node:test";
import {
  requestedHpcTaskCount,
  validateHpcTasksWithinBandCount,
} from "../../src/lib/hpcBandLimits";
import { SlurmResourceRequest } from "../../src/lib/types";

const baseResources: SlurmResourceRequest = {
  resource_type: "cpu",
  ntasks: 4,
};

test("manual band limits allow equal or fewer HPC tasks", () => {
  assert.equal(validateHpcTasksWithinBandCount(baseResources, 4), null);
  assert.equal(validateHpcTasksWithinBandCount(baseResources, 6), null);
});

test("manual band limits reject more HPC tasks than bands", () => {
  assert.match(
    validateHpcTasksWithinBandCount(baseResources, 3, "bands") ?? "",
    /HPC tasks \(4\) cannot exceed manually requested bands \(3\)/,
  );
});

test("manual band limits ignore auto or missing band counts", () => {
  assert.equal(validateHpcTasksWithinBandCount(baseResources, null), null);
  assert.equal(validateHpcTasksWithinBandCount(baseResources, undefined), null);
});

test("requested HPC task count falls back to one", () => {
  assert.equal(requestedHpcTaskCount(null), 1);
  assert.equal(requestedHpcTaskCount({ ...baseResources, ntasks: 0 }), 1);
});
