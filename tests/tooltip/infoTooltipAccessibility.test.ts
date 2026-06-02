import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { join } from "node:path";
import test from "node:test";

const infoTooltipSource = readFileSync(join(process.cwd(), "src/components/InfoTooltip.tsx"), "utf8");

test("icon-only info tooltips are not tabbable", () => {
  assert.equal(infoTooltipSource.includes("tabIndex={children ? undefined : 0}"), false);
  assert.equal(infoTooltipSource.includes('role={children ? undefined : "button"}'), false);
});
