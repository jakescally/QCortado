import type { SavedStructureData } from "./types";

export type MagneticMomentVector = [number, number, number];

export interface AtomicMagneticMoment {
  atomIndex: number;
  vector: MagneticMomentVector;
  magnitude: number;
}

export interface MagnetismViewerData {
  structure: SavedStructureData;
  moments: AtomicMagneticMoment[];
}

export interface MagneticSpeciesFields {
  starting_magnetization?: number;
  theta?: number;
  phi?: number;
}

const NUMBER_PATTERN = "[-+]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[Ee][-+]?\\d+)?";

function asRecord(value: unknown): Record<string, unknown> | null {
  return value && typeof value === "object" && !Array.isArray(value)
    ? value as Record<string, unknown>
    : null;
}

function getFiniteElementValue(params: Record<string, unknown>, keys: string[], element: string): number | undefined {
  for (const key of keys) {
    const map = asRecord(params[key]);
    if (!map) continue;
    const value = Number(map[element]);
    if (Number.isFinite(value)) {
      return value;
    }
  }
  return undefined;
}

export function getMagneticSpeciesFields(
  params: Record<string, unknown>,
  element: string,
): MagneticSpeciesFields {
  const fields: MagneticSpeciesFields = {};
  const startingMagnetization = getFiniteElementValue(params, ["starting_magnetization"], element);
  if (startingMagnetization !== undefined) {
    fields.starting_magnetization = startingMagnetization;
  }

  const theta = getFiniteElementValue(
    params,
    ["starting_magnetization_theta", "starting_magnetization_angle1", "theta", "angle1"],
    element,
  );
  if (theta !== undefined) {
    fields.theta = theta;
  }

  const phi = getFiniteElementValue(
    params,
    ["starting_magnetization_phi", "starting_magnetization_angle2", "phi", "angle2"],
    element,
  );
  if (phi !== undefined) {
    fields.phi = phi;
  }

  return fields;
}

function isMomentVector(value: unknown): value is MagneticMomentVector {
  return (
    Array.isArray(value) &&
    value.length === 3 &&
    value.every((entry) => Number.isFinite(Number(entry)))
  );
}

function vectorMagnitude(vector: MagneticMomentVector): number {
  return Math.hypot(vector[0], vector[1], vector[2]);
}

function isSavedStructureData(value: unknown): value is SavedStructureData {
  const record = asRecord(value);
  return Boolean(
    record
    && Array.isArray(record.atoms)
    && Array.isArray(record.cell_parameters)
    && record.cell_parameters.length === 3,
  );
}

function toCartesianPosition(
  structure: SavedStructureData,
  position: [number, number, number],
): [number, number, number] {
  const cell = structure.cell_parameters;
  if (!cell || cell.length !== 3) {
    return position;
  }
  const [a, b, c] = cell;
  if (structure.position_units === "crystal") {
    return [
      a[0] * position[0] + b[0] * position[1] + c[0] * position[2],
      a[1] * position[0] + b[1] * position[1] + c[1] * position[2],
      a[2] * position[0] + b[2] * position[1] + c[2] * position[2],
    ];
  }
  return position;
}

function invert3x3(m: [[number, number, number], [number, number, number], [number, number, number]]) {
  const [
    [a, b, c],
    [d, e, f],
    [g, h, i],
  ] = m;
  const A = (e * i) - (f * h);
  const B = -((d * i) - (f * g));
  const C = (d * h) - (e * g);
  const D = -((b * i) - (c * h));
  const E = (a * i) - (c * g);
  const F = -((a * h) - (b * g));
  const G = (b * f) - (c * e);
  const H = -((a * f) - (c * d));
  const I = (a * e) - (b * d);
  const det = (a * A) + (b * B) + (c * C);
  if (Math.abs(det) < 1e-12) return null;
  const invDet = 1 / det;
  return [
    [A * invDet, D * invDet, G * invDet],
    [B * invDet, E * invDet, H * invDet],
    [C * invDet, F * invDet, I * invDet],
  ] as [[number, number, number], [number, number, number], [number, number, number]];
}

function cartesianToFractional(
  structure: SavedStructureData,
  cartesian: [number, number, number],
): [number, number, number] | null {
  const cell = structure.cell_parameters;
  if (!cell || cell.length !== 3) {
    return null;
  }
  const [a, b, c] = cell;
  const basis = [
    [a[0], b[0], c[0]],
    [a[1], b[1], c[1]],
    [a[2], b[2], c[2]],
  ] as [[number, number, number], [number, number, number], [number, number, number]];
  const inverse = invert3x3(basis);
  if (!inverse) return null;
  return [
    (inverse[0][0] * cartesian[0]) + (inverse[0][1] * cartesian[1]) + (inverse[0][2] * cartesian[2]),
    (inverse[1][0] * cartesian[0]) + (inverse[1][1] * cartesian[1]) + (inverse[1][2] * cartesian[2]),
    (inverse[2][0] * cartesian[0]) + (inverse[2][1] * cartesian[1]) + (inverse[2][2] * cartesian[2]),
  ];
}

function wrapFractional(value: number): number {
  const wrapped = value - Math.floor(value);
  return wrapped >= 1 - 1e-10 ? 0 : wrapped;
}

function wrappedFractionalDistance(a: [number, number, number], b: [number, number, number]): number {
  const periodicDelta = (x: number) => {
    const abs = Math.abs(x);
    return Math.min(abs, 1 - abs);
  };
  const dx = periodicDelta(a[0] - b[0]);
  const dy = periodicDelta(a[1] - b[1]);
  const dz = periodicDelta(a[2] - b[2]);
  return Math.hypot(dx, dy, dz);
}

function remapMomentsToDisplayStructure(
  displayStructure: SavedStructureData,
  referenceStructure: SavedStructureData,
  referenceMoments: AtomicMagneticMoment[],
): AtomicMagneticMoment[] | null {
  const referenceMomentByAtom = new Map(referenceMoments.map((moment) => [moment.atomIndex, moment]));
  if (referenceMomentByAtom.size === 0) return null;

  const refFractional = referenceStructure.atoms.map((atom) => {
    if (referenceStructure.position_units === "crystal") {
      return [
        wrapFractional(atom.position[0]),
        wrapFractional(atom.position[1]),
        wrapFractional(atom.position[2]),
      ] as [number, number, number];
    }
    const cartesian = toCartesianPosition(referenceStructure, atom.position);
    const fractional = cartesianToFractional(referenceStructure, cartesian);
    if (!fractional) return null;
    return [
      wrapFractional(fractional[0]),
      wrapFractional(fractional[1]),
      wrapFractional(fractional[2]),
    ] as [number, number, number];
  });

  const mappedMoments: AtomicMagneticMoment[] = [];
  const MATCH_TOLERANCE = 5e-3;

  for (let atomIndex = 0; atomIndex < displayStructure.atoms.length; atomIndex += 1) {
    const displayAtom = displayStructure.atoms[atomIndex];
    const displayCartesian = toCartesianPosition(displayStructure, displayAtom.position);
    const displayFractionalInRef = cartesianToFractional(referenceStructure, displayCartesian);
    if (!displayFractionalInRef) return null;
    const wrappedDisplay = [
      wrapFractional(displayFractionalInRef[0]),
      wrapFractional(displayFractionalInRef[1]),
      wrapFractional(displayFractionalInRef[2]),
    ] as [number, number, number];

    let bestMoment: AtomicMagneticMoment | null = null;
    let bestDistance = Number.POSITIVE_INFINITY;

    for (let refIndex = 0; refIndex < referenceStructure.atoms.length; refIndex += 1) {
      const refAtom = referenceStructure.atoms[refIndex];
      if (refAtom.symbol !== displayAtom.symbol) continue;
      const refMoment = referenceMomentByAtom.get(refIndex);
      if (!refMoment) continue;
      const refPos = refFractional[refIndex];
      if (!refPos) continue;
      const distance = wrappedFractionalDistance(wrappedDisplay, refPos);
      if (distance < bestDistance) {
        bestDistance = distance;
        bestMoment = refMoment;
      }
    }

    if (!bestMoment || bestDistance > MATCH_TOLERANCE) {
      return null;
    }

    mappedMoments.push({
      atomIndex,
      vector: bestMoment.vector,
      magnitude: bestMoment.magnitude,
    });
  }

  return mappedMoments;
}

function normalizeMoments(value: unknown): AtomicMagneticMoment[] {
  if (!Array.isArray(value)) return [];
  return value
    .map((entry, index): AtomicMagneticMoment | null => {
      const vector = isMomentVector(entry)
        ? entry.map(Number) as MagneticMomentVector
        : isMomentVector((entry as Record<string, unknown> | null)?.vector)
          ? ((entry as { vector: unknown[] }).vector.map(Number) as MagneticMomentVector)
          : null;
      if (!vector) return null;
      return {
        atomIndex: index,
        vector,
        magnitude: vectorMagnitude(vector),
      };
    })
    .filter((entry): entry is AtomicMagneticMoment => entry !== null);
}

export function parseAtomicMagneticMomentsFromOutput(rawOutput: string): AtomicMagneticMoment[] {
  const vectorRegex = new RegExp(
    `\\batom\\s*:?\\s*(\\d+).*?\\bm(?:agn(?:etization)?)?\\s*=\\s*\\(?\\s*(${NUMBER_PATTERN})\\s+(${NUMBER_PATTERN})\\s+(${NUMBER_PATTERN})\\s*\\)?`,
    "gi",
  );
  const atomNumberRegex = /^\s*atom\s+number\s+(\d+)/i;
  const magnetizationRegex = new RegExp(
    `^\\s*magnetization\\s*:\\s*(${NUMBER_PATTERN})\\s+(${NUMBER_PATTERN})\\s+(${NUMBER_PATTERN})(?:\\s|$)`,
    "i",
  );
  const collinearRegex = new RegExp(
    `\\batom\\s*:?\\s*(\\d+).*?\\bmagn\\s*[:=]\\s*(${NUMBER_PATTERN})`,
    "gi",
  );

  const momentsByAtom = new Map<number, MagneticMomentVector>();
  let pendingAtomIndex: number | null = null;
  for (const line of rawOutput.split(/\r?\n/)) {
    vectorRegex.lastIndex = 0;
    collinearRegex.lastIndex = 0;
    const vectorMatch = vectorRegex.exec(line);
    if (vectorMatch) {
      const atomIndex = Number(vectorMatch[1]) - 1;
      if (atomIndex >= 0) {
        momentsByAtom.set(atomIndex, [
          Number(vectorMatch[2]),
          Number(vectorMatch[3]),
          Number(vectorMatch[4]),
        ]);
      }
      continue;
    }

    const atomNumberMatch = atomNumberRegex.exec(line);
    if (atomNumberMatch) {
      const atomIndex = Number(atomNumberMatch[1]) - 1;
      pendingAtomIndex = atomIndex >= 0 ? atomIndex : null;
      continue;
    }

    const magnetizationMatch = magnetizationRegex.exec(line);
    if (magnetizationMatch && pendingAtomIndex != null) {
      momentsByAtom.set(pendingAtomIndex, [
        Number(magnetizationMatch[1]),
        Number(magnetizationMatch[2]),
        Number(magnetizationMatch[3]),
      ]);
      pendingAtomIndex = null;
      continue;
    }

    const collinearMatch = collinearRegex.exec(line);
    if (collinearMatch) {
      const atomIndex = Number(collinearMatch[1]) - 1;
      if (atomIndex >= 0) {
        momentsByAtom.set(atomIndex, [0, 0, Number(collinearMatch[2])]);
      }
    }
  }

  return Array.from(momentsByAtom.entries())
    .sort(([a], [b]) => a - b)
    .map(([atomIndex, vector]) => ({
      atomIndex,
      vector,
      magnitude: vectorMagnitude(vector),
    }));
}

export function isMagneticScfCalculation(calc: {
  calc_type: string;
  parameters?: Record<string, unknown> | null;
  result?: { converged?: boolean } | null;
}): boolean {
  if (calc.calc_type !== "scf") return false;
  if (calc.result?.converged === false) return false;
  const params = calc.parameters ?? {};
  const nspin = Number(params.nspin);
  return nspin === 2 || nspin === 4 || params.noncolin === true || params.lspinorb === true;
}

export function getMagnetismViewerData(calc: {
  parameters?: Record<string, unknown> | null;
  result?: {
    atomic_magnetic_moments?: unknown;
    raw_output?: string;
  } | null;
}): MagnetismViewerData | null {
  const params = calc.parameters ?? {};
  const sourceStructure = isSavedStructureData(params.source_structure) ? params.source_structure : null;
  const viewerStructure = isSavedStructureData(params.magnetism_viewer_structure)
    ? params.magnetism_viewer_structure
    : null;
  const displayStructure = sourceStructure ?? viewerStructure;
  if (!displayStructure) return null;

  const storedMoments = normalizeMoments(calc.result?.atomic_magnetic_moments);
  const parsedMoments = storedMoments.length > 0
    ? storedMoments
    : parseAtomicMagneticMomentsFromOutput(calc.result?.raw_output ?? "");
  if (parsedMoments.length === 0) return null;

  const referenceStructure = (() => {
    if (viewerStructure && parsedMoments.some((moment) => moment.atomIndex < viewerStructure.atoms.length)) {
      return viewerStructure;
    }
    if (sourceStructure && parsedMoments.some((moment) => moment.atomIndex < sourceStructure.atoms.length)) {
      return sourceStructure;
    }
    return viewerStructure ?? sourceStructure;
  })() ?? displayStructure;

  const referenceMoments = parsedMoments.filter((moment) => moment.atomIndex < referenceStructure.atoms.length);
  if (referenceMoments.length === 0) return null;

  const directDisplayMoments = referenceStructure === displayStructure
    ? referenceMoments
    : null;
  const remappedMoments = directDisplayMoments ?? remapMomentsToDisplayStructure(
    displayStructure,
    referenceStructure,
    referenceMoments,
  );
  if (!remappedMoments || remappedMoments.length === 0) return null;

  return {
    structure: displayStructure,
    moments: remappedMoments,
  };
}
