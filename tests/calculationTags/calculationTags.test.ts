import test from "node:test";
import assert from "node:assert/strict";

import { buildCalculationTagList, normalizeCalculationTagBadges } from "../../src/lib/calculationTags";

test("WIEN2k SCF calculations expose the shared DFT+U feature tag", () => {
  const tags = buildCalculationTagList({
    engine_id: "wien2k",
    calc_type: "scf",
    parameters: {
      initialization: { kMesh: [6, 6, 6] },
      run: { dftU: { enabled: true } },
    },
  });

  assert.ok(tags.includes("DFT+U"));
  assert.deepEqual(
    normalizeCalculationTagBadges(tags).find((tag) => tag.label === "DFT+U"),
    { label: "DFT+U", type: "feature" },
  );
});

test("WIEN2k DFT+U tag inference accepts persisted snake_case settings", () => {
  const tags = buildCalculationTagList({
    engine_id: "wien2k",
    calc_type: "scf",
    parameters: { run: { dft_u: { enabled: true } } },
  });

  assert.ok(tags.includes("DFT+U"));
});

test("WIEN2k SCF calculations without DFT+U do not gain the tag", () => {
  const tags = buildCalculationTagList({
    engine_id: "wien2k",
    calc_type: "scf",
    parameters: { run: { dftU: { enabled: false } } },
  });

  assert.ok(!tags.includes("DFT+U"));
});
