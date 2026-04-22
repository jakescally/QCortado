export const ELEMENT_SYMBOLS = [
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

const ELEMENT_SYMBOL_SET = new Set<string>(ELEMENT_SYMBOLS);

const ELEMENT_SYMBOLS_LONGEST_FIRST = [...ELEMENT_SYMBOLS].sort((a, b) => b.length - a.length);

export function normalizeElementSymbol(symbol: string): string {
  const trimmed = String(symbol || "").trim();
  if (!trimmed) {
    return "";
  }
  if (trimmed.length === 1) {
    return trimmed.toUpperCase();
  }
  return `${trimmed[0].toUpperCase()}${trimmed.slice(1).toLowerCase()}`;
}

export function isElementSymbol(symbol: string): boolean {
  return ELEMENT_SYMBOL_SET.has(normalizeElementSymbol(symbol));
}

export function getLeadingElementSymbol(value: string): string | null {
  const raw = String(value || "").trim();
  for (const symbol of ELEMENT_SYMBOLS_LONGEST_FIRST) {
    if (!raw.startsWith(symbol)) continue;
    const next = raw[symbol.length] ?? "";
    if (!next || !/[a-z]/.test(next)) {
      return symbol;
    }
  }
  return null;
}

export function replaceLeadingElementSymbol(value: string, from: string, to: string): string {
  const normalizedFrom = normalizeElementSymbol(from);
  const normalizedTo = normalizeElementSymbol(to);
  const leading = getLeadingElementSymbol(value);
  if (leading !== normalizedFrom) {
    return value;
  }
  return `${normalizedTo}${String(value).trim().slice(leading.length)}`;
}
