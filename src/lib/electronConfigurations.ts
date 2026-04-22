import { ELEMENT_SYMBOLS, normalizeElementSymbol } from "./elements";

export type SubshellOrbital = "s" | "p" | "d" | "f";

export interface ElectronConfiguration {
  compact: string;
  full: string;
}

export interface ElectronSubshell {
  shell: number;
  orbital: SubshellOrbital;
  electrons: number;
  label: string;
}

interface SubshellDefinition {
  shell: number;
  orbital: SubshellOrbital;
  capacity: number;
}

const SUBSHELL_ORDER: SubshellDefinition[] = [
  { shell: 1, orbital: "s", capacity: 2 },
  { shell: 2, orbital: "s", capacity: 2 },
  { shell: 2, orbital: "p", capacity: 6 },
  { shell: 3, orbital: "s", capacity: 2 },
  { shell: 3, orbital: "p", capacity: 6 },
  { shell: 4, orbital: "s", capacity: 2 },
  { shell: 3, orbital: "d", capacity: 10 },
  { shell: 4, orbital: "p", capacity: 6 },
  { shell: 5, orbital: "s", capacity: 2 },
  { shell: 4, orbital: "d", capacity: 10 },
  { shell: 5, orbital: "p", capacity: 6 },
  { shell: 6, orbital: "s", capacity: 2 },
  { shell: 4, orbital: "f", capacity: 14 },
  { shell: 5, orbital: "d", capacity: 10 },
  { shell: 6, orbital: "p", capacity: 6 },
  { shell: 7, orbital: "s", capacity: 2 },
  { shell: 5, orbital: "f", capacity: 14 },
  { shell: 6, orbital: "d", capacity: 10 },
  { shell: 7, orbital: "p", capacity: 6 },
];

const NOBLE_GAS_CORES = [
  { symbol: "He", atomicNumber: 2 },
  { symbol: "Ne", atomicNumber: 10 },
  { symbol: "Ar", atomicNumber: 18 },
  { symbol: "Kr", atomicNumber: 36 },
  { symbol: "Xe", atomicNumber: 54 },
  { symbol: "Rn", atomicNumber: 86 },
] as const;

const EXCEPTIONAL_COMPACT_CONFIGS: Partial<Record<number, string>> = {
  24: "[Ar] 3d5 4s1",
  29: "[Ar] 3d10 4s1",
  41: "[Kr] 4d4 5s1",
  42: "[Kr] 4d5 5s1",
  44: "[Kr] 4d7 5s1",
  45: "[Kr] 4d8 5s1",
  46: "[Kr] 4d10",
  47: "[Kr] 4d10 5s1",
  57: "[Xe] 5d1 6s2",
  58: "[Xe] 4f1 5d1 6s2",
  64: "[Xe] 4f7 5d1 6s2",
  78: "[Xe] 4f14 5d9 6s1",
  79: "[Xe] 4f14 5d10 6s1",
  89: "[Rn] 6d1 7s2",
  90: "[Rn] 6d2 7s2",
  91: "[Rn] 5f2 6d1 7s2",
  92: "[Rn] 5f3 6d1 7s2",
  93: "[Rn] 5f4 6d1 7s2",
  96: "[Rn] 5f7 6d1 7s2",
};

const symbolToAtomicNumber = new Map<string, number>(
  ELEMENT_SYMBOLS.map((symbol, index) => [symbol, index + 1]),
);

function getAtomicNumber(symbol: string): number | null {
  const normalized = normalizeElementSymbol(symbol);
  return symbolToAtomicNumber.get(normalized) ?? null;
}

function formatSubshellToken(shell: number, orbital: SubshellOrbital, electrons: number): string {
  return `${shell}${orbital}${electrons}`;
}

function buildAufbauConfigurationTokens(atomicNumber: number): string[] {
  let electronsRemaining = atomicNumber;
  const tokens: string[] = [];

  for (const subshell of SUBSHELL_ORDER) {
    if (electronsRemaining <= 0) {
      break;
    }
    const electrons = Math.min(subshell.capacity, electronsRemaining);
    tokens.push(formatSubshellToken(subshell.shell, subshell.orbital, electrons));
    electronsRemaining -= electrons;
  }

  return tokens;
}

function expandCompactConfiguration(compact: string): string[] {
  const tokens: string[] = [];

  for (const piece of compact.split(/\s+/).filter(Boolean)) {
    const coreMatch = piece.match(/^\[([A-Z][a-z]?)\]$/);
    if (coreMatch) {
      const core = getNeutralElectronConfiguration(coreMatch[1]);
      if (core) {
        tokens.push(...core.full.split(/\s+/).filter(Boolean));
      }
      continue;
    }
    tokens.push(piece);
  }

  return tokens;
}

function buildCompactConfiguration(fullTokens: string[], atomicNumber: number): string {
  const nobleGas = [...NOBLE_GAS_CORES]
    .reverse()
    .find((core) => core.atomicNumber < atomicNumber);

  if (!nobleGas) {
    return fullTokens.join(" ");
  }

  const coreTokens = buildAufbauConfigurationTokens(nobleGas.atomicNumber);
  const matchesCore = coreTokens.every((token, index) => fullTokens[index] === token);
  if (!matchesCore || coreTokens.length >= fullTokens.length) {
    return fullTokens.join(" ");
  }

  return [`[${nobleGas.symbol}]`, ...fullTokens.slice(coreTokens.length)].join(" ");
}

function parseSubshellToken(token: string): ElectronSubshell | null {
  const match = token.match(/^(\d+)([spdf])(\d+)$/);
  if (!match) return null;

  const shell = Number(match[1]);
  const orbital = match[2] as SubshellOrbital;
  const electrons = Number(match[3]);
  if (!Number.isFinite(shell) || !Number.isFinite(electrons) || electrons <= 0) {
    return null;
  }

  return {
    shell,
    orbital,
    electrons,
    label: `${shell}${orbital}`,
  };
}

export function getNeutralElectronConfiguration(symbol: string): ElectronConfiguration | null {
  const normalized = normalizeElementSymbol(symbol);
  const atomicNumber = symbolToAtomicNumber.get(normalized);
  if (!atomicNumber) {
    return null;
  }

  const compactException = EXCEPTIONAL_COMPACT_CONFIGS[atomicNumber];
  if (compactException) {
    return {
      compact: compactException,
      full: expandCompactConfiguration(compactException).join(" "),
    };
  }

  const fullTokens = buildAufbauConfigurationTokens(atomicNumber);
  return {
    compact: buildCompactConfiguration(fullTokens, atomicNumber),
    full: fullTokens.join(" "),
  };
}

export function getNeutralElectronSubshells(symbol: string): ElectronSubshell[] {
  const config = getNeutralElectronConfiguration(symbol);
  if (!config) return [];

  return config.full
    .split(/\s+/)
    .map(parseSubshellToken)
    .filter((subshell): subshell is ElectronSubshell => Boolean(subshell));
}

export function getOutermostOccupiedOrbitalManifold(symbol: string): string | null {
  const subshells = getNeutralElectronSubshells(symbol);
  if (subshells.length === 0) return null;

  const highestShell = Math.max(...subshells.map((subshell) => subshell.shell));
  const outerShellSubshells = subshells.filter((subshell) => subshell.shell === highestShell);
  const orbitalRank: Record<SubshellOrbital, number> = { s: 0, p: 1, d: 2, f: 3 };
  outerShellSubshells.sort((a, b) => orbitalRank[b.orbital] - orbitalRank[a.orbital]);

  return outerShellSubshells[0]?.label ?? null;
}

export function getOutermostOccupiedOrbital(symbol: string): SubshellOrbital | null {
  const manifold = getOutermostOccupiedOrbitalManifold(symbol);
  return manifold ? (manifold.slice(-1) as SubshellOrbital) : null;
}

export function getDefaultHubbardManifold(symbol: string): string | null {
  const atomicNumber = getAtomicNumber(symbol);
  if (!atomicNumber) return null;

  if (atomicNumber >= 57 && atomicNumber <= 71) return "4f";
  if (atomicNumber >= 89 && atomicNumber <= 103) return "5f";
  if (atomicNumber >= 21 && atomicNumber <= 30) return "3d";
  if (atomicNumber >= 39 && atomicNumber <= 48) return "4d";
  if (atomicNumber >= 72 && atomicNumber <= 80) return "5d";
  if (atomicNumber >= 104 && atomicNumber <= 112) return "6d";

  return null;
}

export function getDefaultHubbardOrbital(symbol: string): SubshellOrbital | null {
  const manifold = getDefaultHubbardManifold(symbol);
  return manifold ? (manifold.slice(-1) as SubshellOrbital) : null;
}
