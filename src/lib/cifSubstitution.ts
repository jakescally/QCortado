import { parseCIF } from "./cifParser";
import {
  getLeadingElementSymbol,
  isElementSymbol,
  normalizeElementSymbol,
  replaceLeadingElementSymbol,
} from "./elements";

export interface CifSubstitutionMapping {
  from: string;
  to: string;
}

export interface CifAtomSitePreview {
  label: string;
  typeSymbol: string;
  element: string;
  fractX: string;
  fractY: string;
  fractZ: string;
  occupancy: string;
}

export interface CifAtomElementGroup {
  element: string;
  count: number;
  sites: CifAtomSitePreview[];
}

export interface CifSubstitutionChangedSite {
  label: string;
  newLabel: string;
  typeSymbol: string;
  newTypeSymbol: string;
  fractX: string;
  fractY: string;
  fractZ: string;
  occupancy: string;
}

export interface CifSubstitutionResult {
  content: string;
  originalFormula: string | null;
  newFormula: string | null;
  changedSites: CifSubstitutionChangedSite[];
  warnings: string[];
  suggestedFilename: string;
}

export interface CifSubstitutionOptions {
  sourceFilename?: string;
}

interface CifToken {
  value: string;
  rawStart: number;
  rawEnd: number;
  valueStart: number;
  valueEnd: number;
  line: number;
}

interface CifLoop {
  loopTokenIndex: number;
  columns: string[];
  columnTokens: CifToken[];
  rows: CifToken[][];
}

interface Replacement {
  start: number;
  end: number;
  value: string;
}

const FORMULA_TAGS = new Set([
  "_chemical_formula_structural",
  "_chemical_formula_sum",
]);

const TYPE_SYMBOL_COLUMNS = new Set([
  "_atom_site_type_symbol",
  "_atom_type_symbol",
  "_atom_site_aniso_type_symbol",
]);

function tokenizeCif(content: string): CifToken[] {
  const tokens: CifToken[] = [];
  let i = 0;
  let line = 1;

  while (i < content.length) {
    const char = content[i];

    if (char === "\n") {
      line += 1;
      i += 1;
      continue;
    }

    if (/\s/.test(char)) {
      i += 1;
      continue;
    }

    const isLineStart = i === 0 || content[i - 1] === "\n";
    if (isLineStart && char === ";") {
      i += 1;
      while (i < content.length) {
        if (content[i] === "\n") {
          line += 1;
          if (content[i + 1] === ";") {
            i += 2;
            break;
          }
        }
        i += 1;
      }
      continue;
    }

    if (char === "#") {
      while (i < content.length && content[i] !== "\n") {
        i += 1;
      }
      continue;
    }

    if (char === "'" || char === "\"") {
      const quote = char;
      const rawStart = i;
      const valueStart = i + 1;
      i += 1;
      while (i < content.length && content[i] !== quote) {
        if (content[i] === "\n") line += 1;
        i += 1;
      }
      const valueEnd = i;
      const rawEnd = i < content.length ? i + 1 : i;
      tokens.push({
        value: content.slice(valueStart, valueEnd),
        rawStart,
        rawEnd,
        valueStart,
        valueEnd,
        line,
      });
      i = rawEnd;
      continue;
    }

    const rawStart = i;
    while (i < content.length && !/\s/.test(content[i]) && content[i] !== "#") {
      i += 1;
    }
    tokens.push({
      value: content.slice(rawStart, i),
      rawStart,
      rawEnd: i,
      valueStart: rawStart,
      valueEnd: i,
      line,
    });
  }

  return tokens;
}

function parseLoops(tokens: CifToken[]): CifLoop[] {
  const loops: CifLoop[] = [];

  for (let i = 0; i < tokens.length; i += 1) {
    if (tokens[i].value.toLowerCase() !== "loop_") continue;

    const columnTokens: CifToken[] = [];
    let cursor = i + 1;
    while (cursor < tokens.length && tokens[cursor].value.startsWith("_")) {
      columnTokens.push(tokens[cursor]);
      cursor += 1;
    }

    if (columnTokens.length === 0) continue;

    const rows: CifToken[][] = [];
    const columnCount = columnTokens.length;
    while (cursor < tokens.length) {
      const value = tokens[cursor].value;
      const lower = value.toLowerCase();
      if (lower === "loop_" || lower.startsWith("data_") || value.startsWith("_")) {
        break;
      }
      if (cursor + columnCount > tokens.length) {
        break;
      }
      rows.push(tokens.slice(cursor, cursor + columnCount));
      cursor += columnCount;
    }

    loops.push({
      loopTokenIndex: i,
      columns: columnTokens.map((token) => token.value),
      columnTokens,
      rows,
    });
    i = Math.max(i, cursor - 1);
  }

  return loops;
}

function columnIndex(loop: CifLoop, name: string): number {
  const target = name.toLowerCase();
  return loop.columns.findIndex((column) => column.toLowerCase() === target);
}

function getAtomSiteLoop(loops: CifLoop[]): CifLoop | null {
  return loops.find((loop) => columnIndex(loop, "_atom_site_label") >= 0) ?? null;
}

function getSiteElement(row: CifToken[], labelIdx: number, typeIdx: number): string | null {
  if (typeIdx >= 0) {
    const typeElement = getLeadingElementSymbol(row[typeIdx]?.value ?? "");
    if (typeElement) return typeElement;
  }
  if (labelIdx >= 0) {
    return getLeadingElementSymbol(row[labelIdx]?.value ?? "");
  }
  return null;
}

function normalizeMappings(mappings: CifSubstitutionMapping[]): Map<string, string> {
  const normalized = new Map<string, string>();

  for (const mapping of mappings) {
    const from = normalizeElementSymbol(mapping.from);
    const to = normalizeElementSymbol(mapping.to);
    if (!isElementSymbol(from)) {
      throw new Error(`Invalid source element symbol: ${mapping.from}`);
    }
    if (!isElementSymbol(to)) {
      throw new Error(`Invalid replacement element symbol: ${mapping.to}`);
    }
    if (from === to) {
      throw new Error(`Choose a different replacement for ${from}.`);
    }
    if (normalized.has(from)) {
      throw new Error(`Duplicate substitution for ${from}.`);
    }
    normalized.set(from, to);
  }

  if (normalized.size === 0) {
    throw new Error("Add at least one element substitution.");
  }

  return normalized;
}

function replaceFormulaElements(value: string, mapping: Map<string, string>): string {
  return value.replace(/(^|[^A-Za-z])([A-Z][a-z]?)(?=([0-9.+-]|\s|$))/g, (match, prefix: string, symbol: string) => {
    const replacement = mapping.get(symbol);
    if (!replacement) return match;
    return `${prefix}${replacement}`;
  });
}

function sanitizedFormulaFilename(formula: string | null, fallback: string): string {
  const stem = (formula || fallback || "modified-structure")
    .replace(/['"]/g, "")
    .replace(/\s+/g, "")
    .replace(/[^A-Za-z0-9._-]/g, "_")
    .replace(/^_+|_+$/g, "");
  return `${stem || "modified-structure"}.cif`;
}

function addReplacement(replacements: Replacement[], token: CifToken, value: string): void {
  if (token.value === value) return;
  replacements.push({
    start: token.valueStart,
    end: token.valueEnd,
    value,
  });
}

function applyReplacements(content: string, replacements: Replacement[]): string {
  const sorted = [...replacements].sort((a, b) => b.start - a.start);
  let next = content;
  let previousStart = Number.POSITIVE_INFINITY;

  for (const replacement of sorted) {
    if (replacement.end > previousStart) {
      throw new Error("Internal CIF substitution error: overlapping edits.");
    }
    next = `${next.slice(0, replacement.start)}${replacement.value}${next.slice(replacement.end)}`;
    previousStart = replacement.start;
  }

  return next;
}

function insertGeneratedComment(content: string, mappings: Map<string, string>, sourceFilename?: string): string {
  const mappingText = [...mappings.entries()].map(([from, to]) => `${from}->${to}`).join(", ");
  const cleanSource = (sourceFilename || "source CIF").replace(/[\r\n]/g, " ").trim();
  const comment = `# QCortado generated CIF substitution: ${mappingText} from ${cleanSource}\n`;
  const dataMatch = content.match(/^data_[^\r\n]*(?:\r?\n)?/m);
  if (!dataMatch || dataMatch.index == null) {
    return `${comment}${content}`;
  }
  const insertAt = dataMatch.index + dataMatch[0].length;
  return `${content.slice(0, insertAt)}${comment}${content.slice(insertAt)}`;
}

function isAtomSiteLabelReferenceColumn(column: string): boolean {
  const lower = column.toLowerCase();
  return lower === "_atom_site_label"
    || lower === "_atom_site_aniso_label"
    || lower.includes("_atom_site_label");
}

function validateAtomTypeMerges(loops: CifLoop[], mapping: Map<string, string>): void {
  const atomTypeLoop = loops.find((loop) => columnIndex(loop, "_atom_type_symbol") >= 0);
  if (!atomTypeLoop) return;

  const symbolIdx = columnIndex(atomTypeLoop, "_atom_type_symbol");
  const finalToOriginal = new Map<string, string>();

  for (const row of atomTypeLoop.rows) {
    const token = row[symbolIdx];
    const from = getLeadingElementSymbol(token.value);
    const to = from ? mapping.get(from) : undefined;
    const finalValue = to && from ? replaceLeadingElementSymbol(token.value, from, to) : token.value;
    const existing = finalToOriginal.get(finalValue);
    if (existing != null && existing !== token.value) {
      throw new Error(
        `Substitution would create duplicate atom type "${finalValue}". Remove or merge duplicate atom-type rows first.`,
      );
    }
    finalToOriginal.set(finalValue, token.value);
  }
}

function findFormulaTokens(tokens: CifToken[]): Array<{ tag: CifToken; value: CifToken }> {
  const formulas: Array<{ tag: CifToken; value: CifToken }> = [];
  for (let i = 0; i < tokens.length - 1; i += 1) {
    if (FORMULA_TAGS.has(tokens[i].value.toLowerCase())) {
      formulas.push({ tag: tokens[i], value: tokens[i + 1] });
    }
  }
  return formulas;
}

export function inspectCifAtomGroups(content: string): CifAtomElementGroup[] {
  const tokens = tokenizeCif(content);
  const loops = parseLoops(tokens);
  const atomSiteLoop = getAtomSiteLoop(loops);
  if (!atomSiteLoop) {
    return [];
  }

  const labelIdx = columnIndex(atomSiteLoop, "_atom_site_label");
  const typeIdx = columnIndex(atomSiteLoop, "_atom_site_type_symbol");
  const xIdx = columnIndex(atomSiteLoop, "_atom_site_fract_x");
  const yIdx = columnIndex(atomSiteLoop, "_atom_site_fract_y");
  const zIdx = columnIndex(atomSiteLoop, "_atom_site_fract_z");
  const occupancyIdx = columnIndex(atomSiteLoop, "_atom_site_occupancy");

  const groups = new Map<string, CifAtomSitePreview[]>();
  for (const row of atomSiteLoop.rows) {
    const element = getSiteElement(row, labelIdx, typeIdx);
    if (!element) continue;
    const sites = groups.get(element) ?? [];
    sites.push({
      label: labelIdx >= 0 ? row[labelIdx].value : "",
      typeSymbol: typeIdx >= 0 ? row[typeIdx].value : "",
      element,
      fractX: xIdx >= 0 ? row[xIdx].value : "",
      fractY: yIdx >= 0 ? row[yIdx].value : "",
      fractZ: zIdx >= 0 ? row[zIdx].value : "",
      occupancy: occupancyIdx >= 0 ? row[occupancyIdx].value : "",
    });
    groups.set(element, sites);
  }

  return [...groups.entries()]
    .map(([element, sites]) => ({ element, count: sites.length, sites }))
    .sort((a, b) => a.element.localeCompare(b.element));
}

export function substituteCifElements(
  content: string,
  mappings: CifSubstitutionMapping[],
  options: CifSubstitutionOptions = {},
): CifSubstitutionResult {
  const normalizedMappings = normalizeMappings(mappings);
  const tokens = tokenizeCif(content);
  const loops = parseLoops(tokens);
  const atomSiteLoop = getAtomSiteLoop(loops);
  if (!atomSiteLoop) {
    throw new Error("No _atom_site loop was found in this CIF.");
  }

  validateAtomTypeMerges(loops, normalizedMappings);

  const labelIdx = columnIndex(atomSiteLoop, "_atom_site_label");
  const typeIdx = columnIndex(atomSiteLoop, "_atom_site_type_symbol");
  const xIdx = columnIndex(atomSiteLoop, "_atom_site_fract_x");
  const yIdx = columnIndex(atomSiteLoop, "_atom_site_fract_y");
  const zIdx = columnIndex(atomSiteLoop, "_atom_site_fract_z");
  const occupancyIdx = columnIndex(atomSiteLoop, "_atom_site_occupancy");

  const labelRenameMap = new Map<string, string>();
  const changedSites: CifSubstitutionChangedSite[] = [];

  for (const row of atomSiteLoop.rows) {
    const element = getSiteElement(row, labelIdx, typeIdx);
    if (!element) continue;
    const replacement = normalizedMappings.get(element);
    if (!replacement) continue;

    const oldLabel = labelIdx >= 0 ? row[labelIdx].value : "";
    const oldType = typeIdx >= 0 ? row[typeIdx].value : "";
    const newLabel = oldLabel ? replaceLeadingElementSymbol(oldLabel, element, replacement) : oldLabel;
    const newType = oldType ? replaceLeadingElementSymbol(oldType, element, replacement) : oldType;
    if (oldLabel && newLabel !== oldLabel) {
      labelRenameMap.set(oldLabel, newLabel);
    }
    changedSites.push({
      label: oldLabel,
      newLabel,
      typeSymbol: oldType,
      newTypeSymbol: newType,
      fractX: xIdx >= 0 ? row[xIdx].value : "",
      fractY: yIdx >= 0 ? row[yIdx].value : "",
      fractZ: zIdx >= 0 ? row[zIdx].value : "",
      occupancy: occupancyIdx >= 0 ? row[occupancyIdx].value : "",
    });
  }

  if (changedSites.length === 0) {
    throw new Error("None of the requested source elements were found in the CIF atom sites.");
  }

  const replacements: Replacement[] = [];
  for (const loop of loops) {
    loop.rows.forEach((row) => {
      loop.columns.forEach((column, colIdx) => {
        const token = row[colIdx];
        if (!token) return;
        const lower = column.toLowerCase();

        if (isAtomSiteLabelReferenceColumn(lower)) {
          const renamed = labelRenameMap.get(token.value);
          if (renamed) {
            addReplacement(replacements, token, renamed);
          }
        }

        if (TYPE_SYMBOL_COLUMNS.has(lower)) {
          const element = getLeadingElementSymbol(token.value);
          const replacement = element ? normalizedMappings.get(element) : undefined;
          if (replacement && element) {
            addReplacement(replacements, token, replaceLeadingElementSymbol(token.value, element, replacement));
          }
        }
      });
    });
  }

  const formulaTokens = findFormulaTokens(tokens);
  const formulas = new Map<string, { oldValue: string; newValue: string }>();
  for (const { tag, value } of formulaTokens) {
    const newValue = replaceFormulaElements(value.value, normalizedMappings);
    formulas.set(tag.value.toLowerCase(), { oldValue: value.value, newValue });
    addReplacement(replacements, value, newValue);
  }

  const originalFormula =
    formulas.get("_chemical_formula_sum")?.oldValue
    ?? formulas.get("_chemical_formula_structural")?.oldValue
    ?? null;

  let modifiedContent = applyReplacements(content, replacements);
  modifiedContent = insertGeneratedComment(modifiedContent, normalizedMappings, options.sourceFilename);

  const parsed = parseCIF(modifiedContent);
  if (parsed.atom_sites.length === 0) {
    throw new Error("Modified CIF parsed successfully, but no atom sites were found.");
  }

  const newFormula =
    parsed.formula_sum
    ?? parsed.formula_structural
    ?? formulas.get("_chemical_formula_sum")?.newValue
    ?? formulas.get("_chemical_formula_structural")?.newValue
    ?? null;

  const warnings: string[] = [];
  if (formulaTokens.length === 0) {
    warnings.push("No chemical formula tag was found, so only atom-site records were updated.");
  }
  if (labelRenameMap.size === 0) {
    warnings.push("Atom-site labels did not use element prefixes, so labels were left unchanged.");
  }

  return {
    content: modifiedContent,
    originalFormula,
    newFormula,
    changedSites,
    warnings,
    suggestedFilename: sanitizedFormulaFilename(newFormula, options.sourceFilename || "modified-structure"),
  };
}
