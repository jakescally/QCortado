export type WannierIssueSeverity = "info" | "warning" | "error";

export interface WannierQualityIssue {
  code: string;
  severity: WannierIssueSeverity;
  message: string;
}

function severityRank(severity: WannierIssueSeverity): number {
  switch (severity) {
    case "error":
      return 3;
    case "warning":
      return 2;
    default:
      return 1;
  }
}

function pushIssue(
  issues: WannierQualityIssue[],
  code: string,
  severity: WannierIssueSeverity,
  message: string,
): void {
  const normalizedMessage = message.trim();
  if (!normalizedMessage) return;
  if (issues.some((issue) => issue.code === code && issue.message === normalizedMessage)) {
    return;
  }
  issues.push({ code, severity, message: normalizedMessage });
}

function getFiniteNumber(value: unknown): number | null {
  const numeric = Number(value);
  return Number.isFinite(numeric) ? numeric : null;
}

function getWannierBandMinDistance(wannierData: any, fermiEnergy: number): number | null {
  const bandData = wannierData?.band_data;
  if (!bandData || !Array.isArray(bandData.energies)) return null;

  let minDistance = Number.POSITIVE_INFINITY;
  for (const band of bandData.energies) {
    if (!Array.isArray(band)) continue;
    for (const energy of band) {
      const numeric = getFiniteNumber(energy);
      if (numeric == null) continue;
      minDistance = Math.min(minDistance, Math.abs(numeric - fermiEnergy));
    }
  }

  return Number.isFinite(minDistance) ? minDistance : null;
}

function parseRawOutputIssues(rawOutput: string, issues: WannierQualityIssue[]): void {
  if (!rawOutput.trim()) return;

  if (/Maximum number of disentanglement iterations reached/i.test(rawOutput)) {
    pushIssue(
      issues,
      "nonconverged",
      "error",
      "Maximum number of disentanglement iterations reached.",
    );
  }
  if (/Disentanglement convergence criteria not satisfied/i.test(rawOutput)) {
    pushIssue(
      issues,
      "nonconverged",
      "error",
      "Disentanglement convergence criteria not satisfied.",
    );
  }
  if (
    /Maximum number of wannieri[sz]ation iterations reached/i.test(rawOutput)
    || /Maximum number of iterations reached/i.test(rawOutput)
  ) {
    pushIssue(
      issues,
      "nonconverged",
      "error",
      "Maximum number of minimization iterations reached.",
    );
  }
}

export function getWannierQualityIssues(
  wannierData: any,
  rawOutput?: string | null,
  scfFermiEnergy?: number | null,
): WannierQualityIssue[] {
  const issues: WannierQualityIssue[] = [];

  for (const issue of Array.isArray(wannierData?.quality_issues) ? wannierData.quality_issues : []) {
    const severity = String(issue?.severity || "warning").toLowerCase() as WannierIssueSeverity;
    const normalizedSeverity: WannierIssueSeverity =
      severity === "error" || severity === "warning" || severity === "info" ? severity : "warning";
    pushIssue(
      issues,
      String(issue?.code || "warning"),
      normalizedSeverity,
      String(issue?.message || ""),
    );
  }

  const convergence = wannierData?.convergence;
  for (const reason of Array.isArray(convergence?.failure_reasons) ? convergence.failure_reasons : []) {
    pushIssue(issues, "nonconverged", "error", String(reason || ""));
  }
  for (const warning of Array.isArray(convergence?.warnings) ? convergence.warnings : []) {
    pushIssue(issues, "warning", "warning", String(warning || ""));
  }

  const combinedRawOutput = [rawOutput, wannierData?.raw_output]
    .filter((value): value is string => typeof value === "string" && value.trim().length > 0)
    .join("\n");
  parseRawOutputIssues(combinedRawOutput, issues);

  const fermiAlignment = wannierData?.fermi_alignment;
  const alignmentDistance = getFiniteNumber(fermiAlignment?.wannier_min_distance_ev);
  const alignmentBrackets = fermiAlignment?.wannier_brackets_fermi === true;
  if (fermiAlignment?.source_brackets_fermi === true && alignmentBrackets === false) {
    pushIssue(
      issues,
      "misses_source_fermi",
      "error",
      alignmentDistance != null
        ? `Interpolated Wannier manifold misses source states near E_F; nearest Wannier band stays ${alignmentDistance.toFixed(3)} eV away.`
        : "Interpolated Wannier manifold misses source states near E_F.",
    );
  } else if (alignmentBrackets === false && alignmentDistance != null && alignmentDistance > 0.25) {
    pushIssue(
      issues,
      "far_from_fermi",
      "warning",
      `Interpolated Wannier manifold does not reach E_F; nearest Wannier band is ${alignmentDistance.toFixed(3)} eV away.`,
    );
  } else if (!fermiAlignment && scfFermiEnergy != null && Number.isFinite(scfFermiEnergy)) {
    const minDistance = getWannierBandMinDistance(wannierData, scfFermiEnergy);
    if (minDistance != null && minDistance > 0.25) {
      pushIssue(
        issues,
        "far_from_fermi",
        "warning",
        `Interpolated Wannier manifold does not reach E_F; nearest Wannier band is ${minDistance.toFixed(3)} eV away.`,
      );
    }
  }

  if (convergence?.converged === false && !issues.some((issue) => issue.code === "nonconverged")) {
    pushIssue(
      issues,
      "nonconverged",
      "error",
      "Wannier90 did not report a converged solution.",
    );
  }

  return [...issues].sort((left, right) => {
    const severityDelta = severityRank(right.severity) - severityRank(left.severity);
    if (severityDelta !== 0) return severityDelta;
    return left.message.localeCompare(right.message);
  });
}

export function getWannierIssueCounts(
  wannierData: any,
  rawOutput?: string | null,
  scfFermiEnergy?: number | null,
): { errors: number; warnings: number } {
  const issues = getWannierQualityIssues(wannierData, rawOutput, scfFermiEnergy);
  return {
    errors: issues.filter((issue) => issue.severity === "error").length,
    warnings: issues.filter((issue) => issue.severity === "warning").length,
  };
}

export function formatWannierConvergenceFlag(value: boolean | null | undefined): string {
  if (value === true) return "Yes";
  if (value === false) return "No";
  return "N/A";
}
