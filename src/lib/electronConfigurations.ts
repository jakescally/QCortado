type SubshellOrbital = "s" | "p" | "d" | "f";

export interface ElectronConfiguration {
  compact: string;
  full: string;
}

interface SubshellDefinition {
  shell: number;
  orbital: SubshellOrbital;
  capacity: number;
}

const ELEMENT_SYMBOLS = [
  "H", "He",
  "Li", "Be", "B", "C", "N", "O", "F", "Ne",
  "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
  "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
  "Ga", "Ge", "As", "Se", "Br", "Kr",
  "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
  "In", "Sn", "Sb", "Te", "I", "Xe",
  "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
  "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
  "Tl", "Pb", "Bi", "Po", "At", "Rn",
  "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
  "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
  "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
] as const;

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

function normalizeElementSymbol(symbol: string): string {
  const trimmed = String(symbol || "").trim();
  if (!trimmed) {
    return "";
  }
  if (trimmed.length === 1) {
    return trimmed.toUpperCase();
  }
  return `${trimmed[0].toUpperCase()}${trimmed.slice(1).toLowerCase()}`;
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
