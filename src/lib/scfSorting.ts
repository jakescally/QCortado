export type ScfSortMode = "recent" | "best";

export interface ScfLikeCalculation {
  parameters?: any;
  result?: {
    converged?: boolean | null;
  } | null;
  started_at?: string;
  completed_at?: string | null;
}

const SOC_PRIORITY_BOOST = 250;

function getRecencyTimestamp(calc: ScfLikeCalculation): number {
  const completed = calc.completed_at ? Date.parse(calc.completed_at) : Number.NaN;
  if (Number.isFinite(completed)) return completed;
  const started = calc.started_at ? Date.parse(calc.started_at) : Number.NaN;
  if (Number.isFinite(started)) return started;
  return 0;
}

function getMeshProduct(mesh: unknown): number {
  if (!Array.isArray(mesh) || mesh.length !== 3) return 0;
  const values = mesh.map((entry) => Number(entry));
  if (!values.every((value) => Number.isFinite(value) && value > 0)) return 0;
  return values[0] * values[1] * values[2];
}

function getThresholdTightness(value: unknown, maxScore = 20): number {
  const numeric = Number(value);
  if (!Number.isFinite(numeric) || numeric <= 0) return 0;
  return Math.max(0, Math.min(maxScore, -Math.log10(numeric)));
}

function getScfBestScore(calc: ScfLikeCalculation): number {
  const params = calc.parameters || {};
  const convergedBonus = calc.result?.converged ? 100 : 0;
  const socBonus = params.lspinorb ? SOC_PRIORITY_BOOST : 0;
  const kScore = Math.log2(Math.max(1, getMeshProduct(params.kgrid)));
  const convScore = getThresholdTightness(params.conv_thr);
  const ecutScore = Math.log2(Math.max(1, Number(params.ecutwfc) || 1));

  return convergedBonus + (4 * convScore) + (3 * kScore) + socBonus + ecutScore;
}

export function sortScfByMode<T extends ScfLikeCalculation>(
  calculations: T[],
  sortMode: ScfSortMode,
): T[] {
  const sorted = [...calculations];
  sorted.sort((a, b) => {
    if (sortMode === "best") {
      const bestDiff = getScfBestScore(b) - getScfBestScore(a);
      if (Math.abs(bestDiff) > 1e-9) return bestDiff;
    }
    return getRecencyTimestamp(b) - getRecencyTimestamp(a);
  });
  return sorted;
}
