import { detectBravaisLattice } from "./brillouinZone";
import { BravaisLatticeType, getBrillouinZoneData } from "./brillouinZoneData";

interface Fixture {
  id: string;
  latticeType: BravaisLatticeType;
  params: {
    a: number;
    b: number;
    c: number;
    alpha: number;
    beta: number;
    gamma: number;
  };
  expectedLabels: string[];
  expectedMissingLabels?: string[];
  expectedSegments?: [string, string][];
}

const FIXTURES: Fixture[] = [
  {
    id: "hR1-alpha-75",
    latticeType: "hR",
    params: { a: 1, b: 1, c: 1, alpha: 75, beta: 75, gamma: 75 },
    expectedLabels: ["B", "B₁", "P₂"],
    expectedMissingLabels: ["Q₁"],
    expectedSegments: [["L", "B₁"], ["B₁", "B"]],
  },
  {
    id: "hR2-alpha-102",
    latticeType: "hR",
    params: { a: 1, b: 1, c: 1, alpha: 102, beta: 102, gamma: 102 },
    expectedLabels: ["Q", "Q₁", "P₁"],
    expectedMissingLabels: ["B", "B₁"],
    expectedSegments: [["Γ", "P"], ["P", "Z"]],
  },
  {
    id: "oC",
    latticeType: "oC",
    params: { a: 1.0, b: 1.2, c: 1.4, alpha: 90, beta: 90, gamma: 90 },
    expectedLabels: ["A", "A₁", "X₁", "T"],
    expectedSegments: [["Γ", "X"], ["X₁", "A₁"]],
  },
  {
    id: "oI",
    latticeType: "oI",
    params: { a: 1.0, b: 1.15, c: 1.4, alpha: 90, beta: 90, gamma: 90 },
    expectedLabels: ["L₂", "X₁", "Y₁", "W"],
    expectedSegments: [["Γ", "X"], ["Y₁", "Z"]],
  },
  {
    id: "oF1",
    latticeType: "oF",
    params: { a: 0.8, b: 1.1, c: 1.3, alpha: 90, beta: 90, gamma: 90 },
    expectedLabels: ["A", "A₁", "X₁", "T"],
    expectedSegments: [["Γ", "Y"], ["T", "X₁"]],
  },
  {
    id: "oF2",
    latticeType: "oF",
    params: { a: 1.2, b: 1.3, c: 1.6, alpha: 90, beta: 90, gamma: 90 },
    expectedLabels: ["C", "C₁", "D₁", "H₁"],
    expectedSegments: [["Γ", "Y"], ["X", "H₁"]],
  },
  {
    id: "oF3",
    latticeType: "oF",
    params: { a: 0.70710678118, b: 1.0, c: 1.0, alpha: 90, beta: 90, gamma: 90 },
    expectedLabels: ["A", "A₁", "X"],
    expectedMissingLabels: ["X₁"],
    expectedSegments: [["Γ", "Y"], ["X", "A"]],
  },
  {
    id: "mC",
    latticeType: "mC",
    params: { a: 1.0, b: 1.1, c: 1.3, alpha: 90, beta: 105, gamma: 90 },
    expectedLabels: ["N", "N₁", "F", "Y₁"],
  },
];

let fixtureValidationRan = false;

export function validateBrillouinZoneFixtures(): void {
  if (fixtureValidationRan) {
    return;
  }
  fixtureValidationRan = true;

  const errors: string[] = [];

  if (detectBravaisLattice(166) !== "trigonal-R") {
    errors.push("Space group 166 must map to trigonal-R.");
  }

  for (const fixture of FIXTURES) {
    const data = getBrillouinZoneData(fixture.latticeType, fixture.params);
    const labels = new Set(data.points.map((point) => point.label));

    // 1) No duplicate labels.
    if (labels.size !== data.points.length) {
      errors.push(`${fixture.id}: duplicate point labels detected.`);
    }

    // 2) Required labels must exist.
    for (const label of fixture.expectedLabels) {
      if (!labels.has(label)) {
        errors.push(`${fixture.id}: missing required label '${label}'.`);
      }
    }

    // 2b) Coordinates should be finite.
    for (const point of data.points) {
      if (!Number.isFinite(point.coords[0]) || !Number.isFinite(point.coords[1]) || !Number.isFinite(point.coords[2])) {
        errors.push(`${fixture.id}: non-finite coordinate found at '${point.label}'.`);
      }
    }

    // 3) Forbidden labels should be absent for branch-discriminating fixtures.
    for (const label of fixture.expectedMissingLabels ?? []) {
      if (labels.has(label)) {
        errors.push(`${fixture.id}: unexpected label '${label}' found.`);
      }
    }

    // 4) Path segments must refer to valid labels.
    for (const [from, to] of data.recommendedPath) {
      if (!labels.has(from)) {
        errors.push(`${fixture.id}: recommendedPath references unknown start '${from}'.`);
      }
      if (!labels.has(to)) {
        errors.push(`${fixture.id}: recommendedPath references unknown end '${to}'.`);
      }
    }

    // 5) Fixture-specific segment presence checks.
    const segmentSet = new Set(data.recommendedPath.map(([from, to]) => `${from}->${to}`));
    for (const [from, to] of fixture.expectedSegments ?? []) {
      if (!segmentSet.has(`${from}->${to}`)) {
        errors.push(`${fixture.id}: expected segment '${from}->${to}' not found.`);
      }
    }
  }

  if (errors.length > 0) {
    throw new Error(`Brillouin-zone fixture validation failed:\n- ${errors.join("\n- ")}`);
  }
}
