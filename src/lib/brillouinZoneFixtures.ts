import { detectBravaisLattice } from "./brillouinZone";
import {
  BravaisLatticeType,
  BrillouinZoneData,
  BrillouinZoneDataOptions,
  findHighSymmetryPoint,
  getBrillouinZoneData,
  getHighSymmetryPointId,
} from "./brillouinZoneData";

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
  options?: BrillouinZoneDataOptions;
  expectedLabels: string[];
  expectedMissingLabels?: string[];
  expectedSegments?: [string, string][];
}

const FIXTURES: Fixture[] = [
  {
    id: "hR1-alpha-75",
    latticeType: "hR",
    params: { a: 1, b: 1, c: 1, alpha: 75, beta: 75, gamma: 75 },
    options: { rhombohedralConvention: "sc_primitive" },
    expectedLabels: ["B", "B₁", "P₂"],
    expectedMissingLabels: ["Q₁"],
    expectedSegments: [["L", "B₁"], ["B₁", "B"]],
  },
  {
    id: "hR1-alpha-75-bilbao",
    latticeType: "hR",
    params: { a: 1, b: 1, c: 1, alpha: 75, beta: 75, gamma: 75 },
    options: { rhombohedralConvention: "bilbao_hex" },
    expectedLabels: ["F", "F₁", "L", "T"],
    expectedMissingLabels: ["Z", "L₁"],
    expectedSegments: [["Γ", "F"], ["B", "T"], ["L", "P₁"]],
  },
  {
    id: "hR2-alpha-102",
    latticeType: "hR",
    params: { a: 1, b: 1, c: 1, alpha: 102, beta: 102, gamma: 102 },
    options: { rhombohedralConvention: "sc_primitive" },
    expectedLabels: ["Q", "Q₁", "P₁"],
    expectedMissingLabels: ["B", "B₁"],
    expectedSegments: [["Γ", "P"], ["P", "Z"]],
  },
  {
    id: "hR2-alpha-102-bilbao",
    latticeType: "hR",
    params: { a: 1, b: 1, c: 1, alpha: 102, beta: 102, gamma: 102 },
    options: { rhombohedralConvention: "bilbao_hex" },
    expectedLabels: ["F", "L", "Q₁", "T"],
    expectedMissingLabels: ["Z"],
    expectedSegments: [["Γ", "P"], ["P", "T"], ["Q₁", "F"]],
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
    const data = getBrillouinZoneData(fixture.latticeType, fixture.params, fixture.options);
    const labels = new Set(data.points.map((point) => point.label));
    const pointIds = new Set(data.points.map((point) => getHighSymmetryPointId(point)));

    // 1) No duplicate labels.
    if (labels.size !== data.points.length) {
      errors.push(`${fixture.id}: duplicate point labels detected.`);
    }

    // 1b) No duplicate point ids.
    if (pointIds.size !== data.points.length) {
      errors.push(`${fixture.id}: duplicate point ids detected.`);
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

  // Cross-convention checks for rhombohedral branch families.
  const rhombohedralParams = [
    { id: "hR1-cross", params: { a: 1, b: 1, c: 1, alpha: 75, beta: 75, gamma: 75 } },
    { id: "hR2-cross", params: { a: 1, b: 1, c: 1, alpha: 102, beta: 102, gamma: 102 } },
  ];
  for (const entry of rhombohedralParams) {
    const primitiveData = getBrillouinZoneData("hR", entry.params, {
      rhombohedralConvention: "sc_primitive",
    });
    const bilbaoData = getBrillouinZoneData("hR", entry.params, {
      rhombohedralConvention: "bilbao_hex",
    });

    const pointById = (data: BrillouinZoneData, id: string) =>
      data.points.find((point) => getHighSymmetryPointId(point) === id) ?? null;

    for (const point of primitiveData.points) {
      const pointId = getHighSymmetryPointId(point);
      const inBilbao = pointById(bilbaoData, pointId);
      if (!inBilbao) {
        errors.push(`${entry.id}: point id '${pointId}' missing in bilbao convention.`);
        continue;
      }
      if (point.label !== inBilbao.label && !(inBilbao.aliases ?? []).includes(point.label)) {
        errors.push(`${entry.id}: bilbao point '${inBilbao.label}' missing alias '${point.label}'.`);
      }
    }

    for (const point of bilbaoData.points) {
      const pointId = getHighSymmetryPointId(point);
      const inPrimitive = pointById(primitiveData, pointId);
      if (!inPrimitive) {
        errors.push(`${entry.id}: point id '${pointId}' missing in primitive convention.`);
        continue;
      }
      if (point.label !== inPrimitive.label && !(inPrimitive.aliases ?? []).includes(point.label)) {
        errors.push(`${entry.id}: primitive point '${inPrimitive.label}' missing alias '${point.label}'.`);
      }
    }

    const toIdSegmentSet = (data: BrillouinZoneData): Set<string> => {
      const ids = new Set<string>();
      for (const [from, to] of data.recommendedPath) {
        const fromPoint = findHighSymmetryPoint(data, from);
        const toPoint = findHighSymmetryPoint(data, to);
        if (!fromPoint || !toPoint) {
          errors.push(`${entry.id}: failed to resolve segment '${from}->${to}' by label.`);
          continue;
        }
        ids.add(`${getHighSymmetryPointId(fromPoint)}->${getHighSymmetryPointId(toPoint)}`);
      }
      return ids;
    };

    const primitiveSegments = toIdSegmentSet(primitiveData);
    const bilbaoSegments = toIdSegmentSet(bilbaoData);
    if (primitiveSegments.size !== bilbaoSegments.size) {
      errors.push(`${entry.id}: convention segment count mismatch after id normalization.`);
      continue;
    }
    for (const segment of primitiveSegments) {
      if (!bilbaoSegments.has(segment)) {
        errors.push(`${entry.id}: missing cross-convention segment '${segment}'.`);
      }
    }

    interface PathNode {
      pointId: string;
      label: string;
      npoints: number;
    }

    const recommendedPathToNodes = (data: BrillouinZoneData): PathNode[] => {
      const nodes: PathNode[] = [];
      const segmentPoints = 20;
      for (const [from, to] of data.recommendedPath) {
        const fromPoint = findHighSymmetryPoint(data, from);
        const toPoint = findHighSymmetryPoint(data, to);
        if (!fromPoint || !toPoint) {
          errors.push(`${entry.id}: cannot build point-id path for segment '${from}->${to}'.`);
          continue;
        }
        const fromId = getHighSymmetryPointId(fromPoint);
        const toId = getHighSymmetryPointId(toPoint);
        if (nodes.length === 0) {
          nodes.push({
            pointId: fromId,
            label: fromPoint.label,
            npoints: segmentPoints,
          });
        } else if (nodes[nodes.length - 1].pointId !== fromId) {
          nodes[nodes.length - 1].npoints = 0;
          nodes.push({
            pointId: fromId,
            label: fromPoint.label,
            npoints: segmentPoints,
          });
        }
        if (nodes[nodes.length - 1].pointId !== toId) {
          nodes.push({
            pointId: toId,
            label: toPoint.label,
            npoints: segmentPoints,
          });
        } else {
          nodes[nodes.length - 1].npoints = segmentPoints;
        }
      }
      if (nodes.length > 0) {
        nodes[nodes.length - 1].npoints = 0;
      }
      return nodes;
    };

    const remapById = (nodes: PathNode[], target: BrillouinZoneData): PathNode[] => {
      const remapped: PathNode[] = [];
      for (const node of nodes) {
        const match = pointById(target, node.pointId);
        if (!match) {
          errors.push(`${entry.id}: point-id remap failed for '${node.pointId}'.`);
          continue;
        }
        remapped.push({
          pointId: node.pointId,
          label: match.label,
          npoints: node.npoints,
        });
      }
      if (remapped.length > 0) {
        remapped[remapped.length - 1].npoints = 0;
      }
      return remapped;
    };

    const primitiveNodes = recommendedPathToNodes(primitiveData);
    const bilbaoNodes = remapById(primitiveNodes, bilbaoData);
    const roundTrippedNodes = remapById(bilbaoNodes, primitiveData);
    if (primitiveNodes.length !== roundTrippedNodes.length) {
      errors.push(`${entry.id}: point-id round-trip changed path length.`);
      continue;
    }
    for (let i = 0; i < primitiveNodes.length; i++) {
      if (primitiveNodes[i].pointId !== roundTrippedNodes[i].pointId) {
        errors.push(`${entry.id}: point-id round-trip mismatch at index ${i}.`);
      }
      if (primitiveNodes[i].npoints !== roundTrippedNodes[i].npoints) {
        errors.push(`${entry.id}: npoints changed after point-id round-trip at index ${i}.`);
      }
    }
  }

  if (errors.length > 0) {
    throw new Error(`Brillouin-zone fixture validation failed:\n- ${errors.join("\n- ")}`);
  }
}
