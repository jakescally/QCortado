export interface NamedCalculation {
  id: string;
  name?: string | null;
}

export function getCalculationName(calc: NamedCalculation | null | undefined): string | null {
  const name = calc?.name?.trim();
  return name ? name : null;
}

export function formatCalculationSourceLabel(
  calc: NamedCalculation | null | undefined,
  fallbackPrefix = "SCF",
): string {
  if (!calc) return "N/A";
  const shortId = calc.id.slice(0, 8);
  const name = getCalculationName(calc);
  return name ? `${name} (${shortId})` : `${fallbackPrefix} ${shortId}`;
}
