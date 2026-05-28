export interface HubbardOccupationAtomBlock {
  atomIndex: number;
  label: string;
  text: string;
}

export interface HubbardOccupationSection {
  text: string;
  atoms: HubbardOccupationAtomBlock[];
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
