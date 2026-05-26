import assert from "node:assert/strict";
import test from "node:test";
import { getTooltipOpenDelay, TOOLTIP_OPEN_DELAY_MS } from "../../src/components/tooltipTiming";

test("question mark tooltip triggers open immediately", () => {
  assert.equal(getTooltipOpenDelay(false), 0);
});

test("custom tooltip triggers wait briefly before opening", () => {
  assert.equal(getTooltipOpenDelay(true), TOOLTIP_OPEN_DELAY_MS);
});
