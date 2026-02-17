type CalculationParams = Record<string, unknown> | null | undefined;

function usesOptimizedStructure(
  params: CalculationParams,
  tags?: string[] | null,
): boolean {
  const structureSource = params?.["structure_source"];

  if (structureSource === "optimization") {
    return true;
  }

  if (typeof structureSource === "object" && structureSource !== null) {
    const sourceType = (structureSource as Record<string, unknown>)["type"];
    if (sourceType === "optimization") {
      return true;
    }
  }

  // Backward-compatible fallback for legacy records.
  if (Array.isArray(tags)) {
    return tags.includes("geometry") || tags.includes("structure-optimized");
  }

  return false;
}

export function isPhononReadyScf(
  params: CalculationParams,
  tags?: string[] | null,
): boolean {
  const convThr = Number(params?.["conv_thr"]);
  if (!Number.isFinite(convThr) || convThr <= 0) {
    return false;
  }

  return convThr <= 1e-12 && usesOptimizedStructure(params, tags);
}
