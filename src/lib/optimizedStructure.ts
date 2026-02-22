import { QePositionUnit, SavedCellSummary, SavedStructureData } from "./types";

const BOHR_TO_ANGSTROM = 0.529177210903;

type CellMatrix = [[number, number, number], [number, number, number], [number, number, number]];

function getBaseElement(symbol: string): string {
  return symbol.replace(/[\d+-]+$/, "");
}

function isQeUnit(value: string): value is QePositionUnit {
  return value === "alat" || value === "bohr" || value === "angstrom" || value === "crystal";
}

function parseUnitFromHeader(line: string): QePositionUnit | null {
  const match = line.match(/[\(\{]\s*([a-zA-Z]+)/);
  if (!match) return null;
  const unit = match[1].toLowerCase();
  return isQeUnit(unit) ? unit : null;
}

function parseAlatFromHeader(line: string): number | null {
  const match = line.match(/[\(\{]\s*alat\s*=\s*([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)/i);
  if (!match) return null;
  const alatValue = Number.parseFloat(match[1]);
  if (!Number.isFinite(alatValue) || alatValue <= 0) return null;
  return alatValue;
}

function parseTriplet(line: string): [number, number, number] | null {
  const numbers = line.match(/[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?/g);
  if (!numbers || numbers.length < 3) return null;

  const x = Number.parseFloat(numbers[0]);
  const y = Number.parseFloat(numbers[1]);
  const z = Number.parseFloat(numbers[2]);
  if ([x, y, z].some((entry) => Number.isNaN(entry))) return null;
  return [x, y, z];
}

function parseCellBlock(lines: string[], index: number): {
  cellParameters: CellMatrix | null;
  cellUnits: QePositionUnit | null;
  cellAlat: number | null;
} {
  const vectors: [number, number, number][] = [];

  for (let i = index + 1; i < Math.min(lines.length, index + 8); i += 1) {
    const vector = parseTriplet(lines[i]);
    if (!vector) break;
    vectors.push(vector);
    if (vectors.length === 3) break;
  }

  return {
    cellParameters: vectors.length === 3 ? [vectors[0], vectors[1], vectors[2]] : null,
    cellUnits: parseUnitFromHeader(lines[index]),
    cellAlat: parseAlatFromHeader(lines[index]),
  };
}

function parseAtomsBlock(lines: string[], index: number): {
  atoms: { symbol: string; position: [number, number, number] }[];
  positionUnits: QePositionUnit;
} {
  const atoms: { symbol: string; position: [number, number, number] }[] = [];
  const parsedUnit = parseUnitFromHeader(lines[index]);
  const positionUnits: QePositionUnit = parsedUnit || "crystal";

  for (let i = index + 1; i < lines.length; i += 1) {
    const line = lines[i].trim();
    if (!line) break;
    if (
      line.startsWith("CELL_PARAMETERS") ||
      line.startsWith("ATOMIC_POSITIONS") ||
      line.startsWith("End final coordinates")
    ) {
      break;
    }

    const match = line.match(/^([A-Za-z]{1,3}[A-Za-z0-9]*)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)/);
    if (!match) continue;

    const x = Number.parseFloat(match[2]);
    const y = Number.parseFloat(match[3]);
    const z = Number.parseFloat(match[4]);
    if ([x, y, z].some((entry) => Number.isNaN(entry))) continue;

    atoms.push({
      symbol: getBaseElement(match[1]),
      position: [x, y, z],
    });
  }

  return { atoms, positionUnits };
}

function convertCellToAngstrom(
  cellParameters: CellMatrix | null,
  cellUnits: QePositionUnit | null,
  alatInBohr: number | null,
): {
  cellParameters: CellMatrix | null;
  cellUnits: QePositionUnit | null;
} {
  if (!cellParameters || !cellUnits) {
    return { cellParameters, cellUnits };
  }

  if (cellUnits === "angstrom") {
    return { cellParameters, cellUnits };
  }

  let scale: number | null = null;
  if (cellUnits === "bohr") {
    scale = BOHR_TO_ANGSTROM;
  } else if (cellUnits === "alat" && alatInBohr && Number.isFinite(alatInBohr) && alatInBohr > 0) {
    scale = alatInBohr * BOHR_TO_ANGSTROM;
  }

  if (scale === null) {
    return { cellParameters, cellUnits };
  }

  const scaled = cellParameters.map((vector) => [
    vector[0] * scale,
    vector[1] * scale,
    vector[2] * scale,
  ]) as CellMatrix;

  return {
    cellParameters: scaled,
    cellUnits: "angstrom",
  };
}

function asCellMatrix(value: unknown): CellMatrix | null {
  if (!Array.isArray(value) || value.length !== 3) return null;

  const rows: [number, number, number][] = [];
  for (const row of value) {
    if (!Array.isArray(row) || row.length !== 3) return null;
    const x = Number(row[0]);
    const y = Number(row[1]);
    const z = Number(row[2]);
    if (![x, y, z].every((entry) => Number.isFinite(entry))) return null;
    rows.push([x, y, z]);
  }

  return [rows[0], rows[1], rows[2]];
}

export function isSavedStructureData(value: unknown): value is SavedStructureData {
  if (!value || typeof value !== "object") return false;

  const candidate = value as Record<string, unknown>;
  if (typeof candidate.position_units !== "string" || !isQeUnit(candidate.position_units)) {
    return false;
  }

  if (!Array.isArray(candidate.atoms)) {
    return false;
  }
  for (const atom of candidate.atoms) {
    if (!atom || typeof atom !== "object") return false;
    const entry = atom as Record<string, unknown>;
    if (typeof entry.symbol !== "string") return false;
    if (!Array.isArray(entry.position) || entry.position.length !== 3) return false;
    if (!entry.position.every((component) => Number.isFinite(Number(component)))) return false;
  }

  if (candidate.cell_parameters !== null && asCellMatrix(candidate.cell_parameters) === null) {
    return false;
  }

  if (
    candidate.cell_units !== null &&
    (typeof candidate.cell_units !== "string" || !isQeUnit(candidate.cell_units))
  ) {
    return false;
  }

  return true;
}

export function extractOptimizedStructure(
  rawOutput: string,
  fallback: SavedStructureData | null,
): SavedStructureData | null {
  const safeFallback = fallback && isSavedStructureData(fallback) ? fallback : null;
  if (!rawOutput) return safeFallback;

  const lines = rawOutput.split(/\r?\n/);
  const begin = lines.findIndex((line) => line.includes("Begin final coordinates"));
  const end = lines.findIndex((line) => line.includes("End final coordinates"));
  const scope = begin >= 0 && end > begin ? lines.slice(begin, end + 1) : lines;

  let lastCell = safeFallback?.cell_parameters ?? null;
  let lastCellUnits = safeFallback?.cell_units ?? null;
  let lastCellAlat: number | null = null;
  let lastAtoms = safeFallback?.atoms ?? [];
  let lastPositionUnits = safeFallback?.position_units ?? "crystal";
  let sawStructureBlock = false;

  for (let i = 0; i < scope.length; i += 1) {
    const line = scope[i];

    if (line.includes("CELL_PARAMETERS")) {
      const parsed = parseCellBlock(scope, i);
      if (parsed.cellParameters) {
        lastCell = parsed.cellParameters;
        sawStructureBlock = true;
      }
      if (parsed.cellUnits) {
        lastCellUnits = parsed.cellUnits;
        lastCellAlat = parsed.cellUnits === "alat" ? parsed.cellAlat : null;
      }
      if (parsed.cellAlat !== null) {
        lastCellAlat = parsed.cellAlat;
      }
    }

    if (line.includes("ATOMIC_POSITIONS")) {
      const parsed = parseAtomsBlock(scope, i);
      if (parsed.atoms.length > 0) {
        lastAtoms = parsed.atoms;
        lastPositionUnits = parsed.positionUnits;
        sawStructureBlock = true;
      }
    }
  }

  if (!sawStructureBlock && !safeFallback) {
    return null;
  }

  const normalizedCell = convertCellToAngstrom(lastCell, lastCellUnits, lastCellAlat);

  return {
    position_units: lastPositionUnits,
    cell_units: normalizedCell.cellUnits,
    cell_parameters: normalizedCell.cellParameters,
    atoms: lastAtoms,
  };
}

export function summarizeCell(structure: SavedStructureData): SavedCellSummary | null {
  if (!structure.cell_parameters) return null;

  const [aVec, bVec, cVec] = structure.cell_parameters;
  const norm = (v: [number, number, number]) => Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  const dot = (u: [number, number, number], v: [number, number, number]) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  const cross = (u: [number, number, number], v: [number, number, number]): [number, number, number] => [
    u[1] * v[2] - u[2] * v[1],
    u[2] * v[0] - u[0] * v[2],
    u[0] * v[1] - u[1] * v[0],
  ];
  const angle = (u: [number, number, number], v: [number, number, number]) => {
    const cosine = dot(u, v) / (norm(u) * norm(v));
    const clamped = Math.max(-1, Math.min(1, cosine));
    return (Math.acos(clamped) * 180) / Math.PI;
  };

  const a = norm(aVec);
  const b = norm(bVec);
  const c = norm(cVec);
  const alpha = angle(bVec, cVec);
  const beta = angle(aVec, cVec);
  const gamma = angle(aVec, bVec);
  const volume = Math.abs(dot(aVec, cross(bVec, cVec)));

  return {
    a,
    b,
    c,
    alpha,
    beta,
    gamma,
    volume,
    units: structure.cell_units || "angstrom",
  };
}
