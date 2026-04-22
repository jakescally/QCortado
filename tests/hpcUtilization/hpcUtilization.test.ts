import assert from "node:assert/strict";
import test from "node:test";
import {
  calculateCpuUtilizationPercent,
  formatTelemetryBytes,
  formatTelemetryPercent,
  resolveTelemetryResourceType,
} from "../../src/lib/hpcUtilization";

test("CPU utilization uses CPU-time delta over allocated CPU wall time", () => {
  const percent = calculateCpuUtilizationPercent(
    {
      jobId: "123",
      capturedAt: "2026-04-22T12:00:00Z",
      totalCpuSeconds: 120,
      allocatedCpus: 8,
    },
    {
      jobId: "123",
      capturedAt: "2026-04-22T12:00:05Z",
      totalCpuSeconds: 140,
      allocatedCpus: 8,
    },
  );
  assert.equal(percent, 50);
});

test("CPU utilization rejects incompatible or invalid samples", () => {
  assert.equal(calculateCpuUtilizationPercent(null, null), null);
  assert.equal(
    calculateCpuUtilizationPercent(
      { jobId: "123", capturedAt: "2026-04-22T12:00:00Z", totalCpuSeconds: 120, allocatedCpus: 4 },
      { jobId: "456", capturedAt: "2026-04-22T12:00:05Z", totalCpuSeconds: 140, allocatedCpus: 4 },
    ),
    null,
  );
  assert.equal(
    calculateCpuUtilizationPercent(
      { jobId: "123", capturedAt: "2026-04-22T12:00:05Z", totalCpuSeconds: 140, allocatedCpus: 4 },
      { jobId: "123", capturedAt: "2026-04-22T12:00:00Z", totalCpuSeconds: 120, allocatedCpus: 4 },
    ),
    null,
  );
});

test("resource type gating prefers task metadata with fallback", () => {
  assert.equal(resolveTelemetryResourceType("gpu", "cpu"), "gpu");
  assert.equal(resolveTelemetryResourceType("cpu", "gpu"), "cpu");
  assert.equal(resolveTelemetryResourceType(null, "gpu"), "gpu");
  assert.equal(resolveTelemetryResourceType(null, null), "cpu");
});

test("telemetry formatters handle missing and finite values", () => {
  assert.equal(formatTelemetryBytes(null), "Unavailable");
  assert.equal(formatTelemetryBytes(1024 * 1024 * 3), "3.0 MB");
  assert.equal(formatTelemetryPercent(null), "Unavailable");
  assert.equal(formatTelemetryPercent(42.2), "42%");
});
