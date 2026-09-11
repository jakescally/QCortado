export interface HubbardOccupationAtomBlock {
  atomIndex: number;
  label: string;
  text: string;
}

export interface HubbardOccupationSection {
  text: string;
  atoms: HubbardOccupationAtomBlock[];
}

interface ComplexValue {
  re: number;
  im: number;
}

interface Wien2kDmatBlock {
  atomIndex: number;
  orbitalL: number;
  spin: "up" | "down";
  matrix: ComplexValue[][];
}

export interface Wien2kDmatFile {
  path: string;
  contents: string;
}

const HUBBARD_OCCUPATIONS_MARKER = "HUBBARD OCCUPATIONS";
const ATOM_HEADER_RE = /^\s*-+\s*ATOM\s+(\d+)\s+-+\s*$/i;
const HUBBARD_ENERGY_RE = /^\s*HUBBARD ENERGY\b/i;

function getLastMarkerIndex(text: string, marker: string): number {
  return text.lastIndexOf(marker);
}

export function parseLatestHubbardOccupations(output: string): HubbardOccupationSection | null {
  const markerIndex = getLastMarkerIndex(output, HUBBARD_OCCUPATIONS_MARKER);
  if (markerIndex < 0) return null;

  const rawSectionText = output.slice(markerIndex);
  const rawLines = rawSectionText.split(/\r?\n/);
  const endIndex = rawLines.findIndex((line) => HUBBARD_ENERGY_RE.test(line));
  const lines = (endIndex >= 0 ? rawLines.slice(0, endIndex + 1) : rawLines).filter((line, index) => {
    if (index === 0) return true;
    return line.length > 0 || rawLines[index - 1].length > 0;
  });
  const sectionText = lines.join("\n").trimEnd();

  const atoms: HubbardOccupationAtomBlock[] = [];
  let currentAtomIndex: number | null = null;
  let currentLines: string[] = [];

  function flushCurrentBlock() {
    if (currentAtomIndex == null || currentLines.length === 0) return;
    const text = currentLines.join("\n").trimEnd();
    atoms.push({
      atomIndex: currentAtomIndex,
      label: `Atom ${currentAtomIndex}`,
      text,
    });
  }

  for (const line of lines) {
    const atomMatch = line.match(ATOM_HEADER_RE);
    if (atomMatch) {
      flushCurrentBlock();
      currentAtomIndex = Number(atomMatch[1]);
      currentLines = [line];
      continue;
    }

    if (currentAtomIndex != null) {
      currentLines.push(line);
    }
  }

  flushCurrentBlock();

  if (atoms.length === 0) return null;

  return {
    text: sectionText,
    atoms,
  };
}

const WIEN2K_DMAT_HEADER_RE = /^\s*(\d+)\s+atom\s+density\s+matrix\s*$/i;
const WIEN2K_NUMBER_RE = /[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?/g;

function parseWien2kNumber(value: string): number {
  return Number(value.replace(/[Dd]/g, "E"));
}

function parseWien2kDmat(contents: string, spin: "up" | "down"): Wien2kDmatBlock[] {
  const lines = contents.split(/\r?\n/);
  const headers = lines.flatMap((line, index) => (
    WIEN2K_DMAT_HEADER_RE.test(line) ? [index] : []
  ));
  const blocks: Wien2kDmatBlock[] = [];

  for (let headerPosition = 0; headerPosition < headers.length; headerPosition += 1) {
    const headerIndex = headers[headerPosition];
    const headerMatch = lines[headerIndex].match(WIEN2K_DMAT_HEADER_RE);
    const metadataLine = lines[headerIndex + 1] ?? "";
    const metadataValues = metadataLine.match(WIEN2K_NUMBER_RE) ?? [];
    const atomIndex = Number(headerMatch?.[1]);
    const orbitalL = metadataValues.length > 0
      ? Math.trunc(parseWien2kNumber(metadataValues[0] ?? ""))
      : -1;
    const dimension = 2 * orbitalL + 1;
    if (!Number.isInteger(atomIndex) || atomIndex <= 0 || orbitalL < 0 || dimension <= 0) continue;

    const nextHeaderIndex = headers[headerPosition + 1] ?? lines.length;
    const payload = lines.slice(headerIndex + 2, nextHeaderIndex).join("\n");
    const values = (payload.match(WIEN2K_NUMBER_RE) ?? []).map(parseWien2kNumber);
    const expectedValues = 2 * dimension * dimension;
    if (values.length < expectedValues || values.slice(0, expectedValues).some((value) => !Number.isFinite(value))) {
      continue;
    }

    const matrix: ComplexValue[][] = [];
    let valueIndex = 0;
    for (let row = 0; row < dimension; row += 1) {
      const matrixRow: ComplexValue[] = [];
      for (let column = 0; column < dimension; column += 1) {
        matrixRow.push({ re: values[valueIndex], im: values[valueIndex + 1] });
        valueIndex += 2;
      }
      matrix.push(matrixRow);
    }
    blocks.push({ atomIndex, orbitalL, spin, matrix });
  }

  return blocks;
}

function jacobiEigenvalues(matrix: number[][]): number[] {
  const size = matrix.length;
  const values = matrix.map((row) => [...row]);
  const tolerance = 1e-12;
  const maxIterations = Math.max(50, 100 * size * size);

  for (let iteration = 0; iteration < maxIterations; iteration += 1) {
    let p = 0;
    let q = 0;
    let largest = 0;
    for (let row = 0; row < size; row += 1) {
      for (let column = row + 1; column < size; column += 1) {
        const magnitude = Math.abs(values[row][column]);
        if (magnitude > largest) {
          largest = magnitude;
          p = row;
          q = column;
        }
      }
    }
    if (largest < tolerance) break;

    const angle = 0.5 * Math.atan2(2 * values[p][q], values[q][q] - values[p][p]);
    const cosine = Math.cos(angle);
    const sine = Math.sin(angle);
    const app = values[p][p];
    const aqq = values[q][q];
    const apq = values[p][q];

    for (let index = 0; index < size; index += 1) {
      if (index === p || index === q) continue;
      const aip = values[index][p];
      const aiq = values[index][q];
      values[index][p] = cosine * aip - sine * aiq;
      values[p][index] = values[index][p];
      values[index][q] = sine * aip + cosine * aiq;
      values[q][index] = values[index][q];
    }
    values[p][p] = cosine * cosine * app - 2 * sine * cosine * apq + sine * sine * aqq;
    values[q][q] = sine * sine * app + 2 * sine * cosine * apq + cosine * cosine * aqq;
    values[p][q] = 0;
    values[q][p] = 0;
  }

  return values.map((row, index) => row[index]).sort((left, right) => left - right);
}

function hermitianEigenvalues(matrix: ComplexValue[][]): number[] {
  const size = matrix.length;
  const realSize = 2 * size;
  const represented = Array.from({ length: realSize }, () => Array(realSize).fill(0));

  for (let row = 0; row < size; row += 1) {
    for (let column = 0; column < size; column += 1) {
      const direct = matrix[row][column];
      const conjugate = matrix[column][row];
      const real = 0.5 * (direct.re + conjugate.re);
      const imaginary = 0.5 * (direct.im - conjugate.im);
      represented[row][column] = real;
      represented[row][column + size] = -imaginary;
      represented[row + size][column] = imaginary;
      represented[row + size][column + size] = real;
    }
  }

  const duplicated = jacobiEigenvalues(represented);
  return Array.from({ length: size }, (_, index) => duplicated[index * 2]);
}

function formatWien2kMatrixBlock(block: Wien2kDmatBlock): string {
  const eigenvalues = hermitianEigenvalues(block.matrix);
  const trace = block.matrix.reduce((sum, row, index) => sum + row[index].re, 0);
  const matrixRows = block.matrix.map((row) => row.map((value) => (
    `${value.re.toFixed(6)} ${value.im < 0 ? "-" : "+"} ${Math.abs(value.im).toFixed(6)}i`
  )).join("   "));

  return [
    `WIEN2k density matrix: l = ${block.orbitalL}, spin ${block.spin}`,
    `Tr[n] = ${trace.toFixed(6)}`,
    "eigenvalues:",
    `  ${eigenvalues.map((value) => value.toFixed(6)).join("  ")}`,
    "occupation matrix (real + imaginary i):",
    ...matrixRows.map((row) => `  ${row}`),
  ].join("\n");
}

export function parseWien2kHubbardOccupations(files: Wien2kDmatFile[]): HubbardOccupationSection | null {
  const blocks = files.flatMap((file) => {
    const lowerPath = file.path.toLowerCase();
    if (lowerPath.endsWith(".dmatup")) return parseWien2kDmat(file.contents, "up");
    if (lowerPath.endsWith(".dmatdn")) return parseWien2kDmat(file.contents, "down");
    return [];
  });
  if (blocks.length === 0) return null;

  blocks.sort((left, right) => (
    left.atomIndex - right.atomIndex
    || left.orbitalL - right.orbitalL
    || (left.spin === right.spin ? 0 : left.spin === "up" ? -1 : 1)
  ));
  const grouped = new Map<number, Wien2kDmatBlock[]>();
  for (const block of blocks) {
    const atomBlocks = grouped.get(block.atomIndex) ?? [];
    atomBlocks.push(block);
    grouped.set(block.atomIndex, atomBlocks);
  }

  const atoms = [...grouped.entries()].map(([atomIndex, atomBlocks]) => ({
    atomIndex,
    label: `Atom ${atomIndex}`,
    text: [
      `------------------------ ATOM ${atomIndex.toString().padStart(4, " ")} ------------------------`,
      ...atomBlocks.map(formatWien2kMatrixBlock),
    ].join("\n\n"),
  }));

  return {
    text: atoms.map((atom) => atom.text).join("\n\n"),
    atoms,
  };
}
