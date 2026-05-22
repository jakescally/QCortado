import { getDefaultHubbardManifold } from "../../electronConfigurations";

export interface HubbardParameterLike {
  parameter?: string | null;
  manifold?: string | null;
  value?: number | string | null;
}

export interface HubbardCalculationLike {
  id?: string;
  calc_type?: string;
  parameters?: Record<string, any> | null;
  result?: Record<string, any> | null;
  started_at?: string | null;
  completed_at?: string | null;
}

export interface HubbardRecommendation {
  element: string;
  manifold: string;
  reason: string;
}

export interface HubbardUDefault {
  value: number;
  source: "lrt" | "general_guess";
  label: string;
  lrtCalcId?: string;
}

export interface HubbardLrtValue {
  element: string;
  manifold: string;
  value: number;
  calcId: string;
  completedAt: string;
}

export interface HubbardUDisplayValue {
  target: string;
  value_ev: number;
}

export const GENERAL_HUBBARD_U_GUESS_EV = 6.0;

export const HUBBARD_J_SOURCE =
  "Moore et al., Phys. Rev. Materials 8, 014409 (2024), Table IX";

// Conservative subset of transition-metal oxide Hund J values from the cited
// source. Uncovered elements intentionally do not get an automatic J.
const HUND_J_DEFAULTS_EV: Record<string, number> = {
  V: 0.849,
  Cr: 0.544,
  Mn: 0.767,
  Fe: 0.437,
  Co: 0.677,
  Ni: 0.770,
  Cu: 0.735,
};

function normalizeElementSymbol(symbol: string): string {
  const trimmed = String(symbol || "").trim().replace(/[\d+-]+$/, "");
  if (!trimmed) return "";
  return trimmed.charAt(0).toUpperCase() + trimmed.slice(1).toLowerCase();
}

export function normalizeHubbardTarget(raw: string): string {
  return String(raw || "")
    .trim()
    .replace(/\s+/g, "")
    .replace(/^U\(/i, "")
    .replace(/\)$/, "");
}

export function splitHubbardTarget(target: string): { element: string; manifold: string } | null {
  const normalized = normalizeHubbardTarget(target);
  const match = normalized.match(/^([A-Z][a-z]?)-?(\d+[spdf])$/);
  if (!match) return null;
  return {
    element: normalizeElementSymbol(match[1]),
    manifold: match[2],
  };
}

export function normalizeHubbardLrtUValues(rawValues: any[] | null | undefined): Array<{
  element: string;
  manifold: string;
  target: string;
  value_ev: number;
}> {
  if (!Array.isArray(rawValues)) return [];
  const values = new Map<string, { element: string; manifold: string; target: string; value_ev: number }>();
  for (const raw of rawValues) {
    const target =
      splitHubbardTarget(String(raw?.target || ""))
      ?? splitHubbardTarget(`${raw?.element || ""}-${raw?.manifold || ""}`)
      ?? splitHubbardTarget(String(raw?.manifold || ""));
    const value = Number(raw?.value_ev ?? raw?.value ?? raw?.u_ev);
    if (!target || !Number.isFinite(value)) continue;
    const normalized = {
      element: target.element,
      manifold: target.manifold,
      target: `${target.element}-${target.manifold}`,
      value_ev: value,
    };
    values.set(normalized.target, normalized);
  }
  return Array.from(values.values());
}

export function getHubbardRecommendations(elements: string[]): HubbardRecommendation[] {
  const seen = new Set<string>();
  const recommendations: HubbardRecommendation[] = [];
  for (const rawElement of elements) {
    const element = normalizeElementSymbol(rawElement);
    if (!element || seen.has(element)) continue;
    seen.add(element);
    const manifold = getDefaultHubbardManifold(element);
    if (!manifold) continue;
    const orbital = manifold.slice(-1);
    const reason =
      orbital === "f"
        ? `${element} has localized f states where DFT+U is commonly used.`
        : `${element} is a transition-metal element with localized d states.`;
    recommendations.push({ element, manifold, reason });
  }
  return recommendations;
}

export function getRecommendedHubbardElements(elements: string[]): string[] {
  return getHubbardRecommendations(elements).map((entry) => entry.element);
}

export function getHundJDefaultEv(element: string): number | null {
  const normalized = normalizeElementSymbol(element);
  return HUND_J_DEFAULTS_EV[normalized] ?? null;
}

function getHubbardParameters(calc: HubbardCalculationLike): HubbardParameterLike[] {
  const raw = calc.parameters?.hubbard_parameters;
  return Array.isArray(raw) ? raw : [];
}

export function getScfHubbardUDisplayValues(calc: HubbardCalculationLike): HubbardUDisplayValue[] {
  const params = calc.parameters || {};
  if (!params.lda_plus_u) return [];

  const values = new Map<string, HubbardUDisplayValue>();
  for (const param of getHubbardParameters(calc)) {
    if (String(param.parameter || "").trim().toUpperCase() !== "U") continue;
    const target = normalizeHubbardTarget(String(param.manifold || ""));
    const value = Number(param.value);
    if (!target || !Number.isFinite(value)) continue;
    values.set(target, { target, value_ev: value });
  }

  if (values.size > 0) {
    return Array.from(values.values());
  }

  const hubbardU = params.hubbard_u;
  if (!hubbardU || typeof hubbardU !== "object" || Array.isArray(hubbardU)) {
    return [];
  }

  const manifolds = params.hubbard_manifold && typeof params.hubbard_manifold === "object" && !Array.isArray(params.hubbard_manifold)
    ? params.hubbard_manifold
    : {};
  for (const [rawElement, rawValue] of Object.entries(hubbardU)) {
    const element = normalizeElementSymbol(rawElement);
    const value = Number(rawValue);
    if (!element || !Number.isFinite(value)) continue;
    const manifold = normalizeHubbardTarget(String((manifolds as Record<string, any>)[rawElement] || ""));
    const target = manifold ? `${element}-${manifold}` : element;
    values.set(target, { target, value_ev: value });
  }

  return Array.from(values.values());
}

export function isDudarevDftUScf(calc: HubbardCalculationLike): boolean {
  if (calc.calc_type !== "scf") return false;
  if (!calc.result?.converged) return false;
  const params = calc.parameters || {};
  if (params.lda_plus_u !== true) return false;
  if (Number(params.hubbard_formulation ?? 0) !== 0) return false;
  const parameters = getHubbardParameters(calc);
  if (parameters.length === 0) return false;
  return parameters.every((param) => String(param.parameter || "").trim().toUpperCase() === "U");
}

export function getHubbardLrtValues(calculations: HubbardCalculationLike[]): HubbardLrtValue[] {
  const values: HubbardLrtValue[] = [];
  for (const calc of calculations) {
    if (calc.calc_type !== "hubbard_lrt" || !calc.result?.converged) continue;
    const completedAt = calc.completed_at || calc.started_at || "";
    const data = calc.result.hubbard_lrt_data ?? calc.result;
    const rawValues = normalizeHubbardLrtUValues(data?.u_values);
    for (const raw of rawValues) {
      values.push({
        element: raw.element,
        manifold: raw.manifold,
        value: raw.value_ev,
        calcId: calc.id || "",
        completedAt,
      });
    }
  }
  return values;
}

export function resolveHubbardUDefault(
  element: string,
  manifold: string,
  calculations: HubbardCalculationLike[],
): HubbardUDefault {
  const normalizedElement = normalizeElementSymbol(element);
  const normalizedManifold = normalizeHubbardTarget(manifold);
  const matching = getHubbardLrtValues(calculations)
    .filter((entry) => entry.element === normalizedElement && entry.manifold === normalizedManifold)
    .sort((a, b) => Date.parse(b.completedAt || "0") - Date.parse(a.completedAt || "0"));

  if (matching[0]) {
    return {
      value: matching[0].value,
      source: "lrt",
      label: `From Hubbard LRT (${matching[0].value.toFixed(3)} eV)`,
      lrtCalcId: matching[0].calcId,
    };
  }

  return {
    value: GENERAL_HUBBARD_U_GUESS_EV,
    source: "general_guess",
    label: "General guess: 6.0 eV",
  };
}

export function getHubbardEligibilityReason(calc: HubbardCalculationLike): string | null {
  if (calc.calc_type !== "scf") return "Only SCF calculations can seed Hubbard LRT.";
  if (!calc.result?.converged) return "SCF did not converge.";
  const params = calc.parameters || {};
  if (params.lda_plus_u !== true) return "Requires an SCF with Dudarev DFT+U enabled.";
  const formulation = Number(params.hubbard_formulation ?? 0);
  if (formulation === 1) return "DFT+U+J SCFs are not supported by hp.x LRT.";
  if (formulation === 2) return "DFT+U+J0 SCFs are not supported by hp.x LRT.";
  if (formulation !== 0) return "Unsupported Hubbard formulation.";
  const parameters = getHubbardParameters(calc);
  if (parameters.length === 0) return "No saved Hubbard U parameters were found.";
  if (!parameters.every((param) => String(param.parameter || "").trim().toUpperCase() === "U")) {
    return "Only pure U Hubbard parameters are supported.";
  }
  return null;
}
