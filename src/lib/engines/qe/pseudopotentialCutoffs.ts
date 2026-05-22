import { PseudopotentialMetadata } from "../../types";

export interface SSSPElementData {
  filename: string;
  md5?: string;
  pseudopotential?: string;
  cutoff_wfc: number;
  cutoff_rho: number;
}

export type CutoffStatusKind = "parsed" | "inferred" | "missing";
export type CutoffStatus = "idle" | CutoffStatusKind;
export type CutoffProvenance = "sssp" | "djrepo" | "upf" | "upf_info" | "upf_fallback" | "unknown" | "missing";
export type CutoffDerivation = "direct" | "from_wfc" | "from_rho" | "missing";

export interface ResolvedPseudoCutoffs {
  wfc: number | null;
  rho: number | null;
  wfcStatus: CutoffStatusKind;
  rhoStatus: CutoffStatusKind;
  wfcProvenance: CutoffProvenance;
  rhoProvenance: CutoffProvenance;
  wfcDerivation: CutoffDerivation;
  rhoDerivation: CutoffDerivation;
  wfcRatio: number | null;
  rhoRatio: number | null;
}

export interface SelectedPseudoCutoffSummary {
  maxWfc: number;
  maxRho: number;
  wfcStatus: CutoffStatus;
  rhoStatus: CutoffStatus;
  wfcProvenance: CutoffProvenance;
  rhoProvenance: CutoffProvenance;
  wfcDerivation: CutoffDerivation;
  rhoDerivation: CutoffDerivation;
  wfcRatio: number | null;
  rhoRatio: number | null;
  hasInferredCutoff: boolean;
  hasMissingCutoff: boolean;
}

export function getPseudoCutoffRatio(pseudo: PseudopotentialMetadata | undefined): number {
  const pseudoType = (pseudo?.pseudo_type || "").toUpperCase();
  if (pseudoType.includes("US")) {
    return 8;
  }
  return 4;
}

export function inferChargeDensityCutoff(
  pseudo: PseudopotentialMetadata | undefined,
  wavefunctionCutoff: number | null,
): number | null {
  if (!pseudo || wavefunctionCutoff == null || !Number.isFinite(wavefunctionCutoff) || wavefunctionCutoff <= 0) {
    return null;
  }

  return wavefunctionCutoff * getPseudoCutoffRatio(pseudo);
}

export function inferWavefunctionCutoff(
  pseudo: PseudopotentialMetadata | undefined,
  chargeDensityCutoff: number | null,
): number | null {
  if (!pseudo || chargeDensityCutoff == null || !Number.isFinite(chargeDensityCutoff) || chargeDensityCutoff <= 0) {
    return null;
  }

  return chargeDensityCutoff / getPseudoCutoffRatio(pseudo);
}

export function normalizeCutoff(value: number | null | undefined): number | null {
  return Number.isFinite(value) && (value ?? 0) > 0 ? value as number : null;
}

export function normalizeCutoffProvenance(value: string | null | undefined): CutoffProvenance {
  switch (value) {
    case "djrepo":
    case "upf":
    case "upf_info":
    case "upf_fallback":
      return value;
    default:
      return value ? "unknown" : "missing";
  }
}

export function classifyPreferredCutoffStatus(provenance: CutoffProvenance): CutoffStatusKind {
  return provenance === "upf_fallback" ? "inferred" : "parsed";
}

export function resolvePseudoCutoffs(
  selectedFilename: string | null | undefined,
  pseudo: PseudopotentialMetadata | undefined,
  ssspEntry?: SSSPElementData | null,
): ResolvedPseudoCutoffs {
  const ssspMatchesPseudo = Boolean(
    ssspEntry
    && (
      (pseudo && pseudo.filename.toLowerCase() === ssspEntry.filename.toLowerCase())
      || (selectedFilename && selectedFilename.toLowerCase() === ssspEntry.filename.toLowerCase())
    ),
  );
  const preferredWfc = normalizeCutoff((ssspMatchesPseudo ? ssspEntry?.cutoff_wfc : null) ?? pseudo?.cutoff_wfc ?? null);
  const preferredRho = normalizeCutoff((ssspMatchesPseudo ? ssspEntry?.cutoff_rho : null) ?? pseudo?.cutoff_rho ?? null);
  const preferredWfcProvenance = ssspMatchesPseudo && normalizeCutoff(ssspEntry?.cutoff_wfc) != null
    ? "sssp"
    : normalizeCutoffProvenance(pseudo?.cutoff_wfc_source);
  const preferredRhoProvenance = ssspMatchesPseudo && normalizeCutoff(ssspEntry?.cutoff_rho) != null
    ? "sssp"
    : normalizeCutoffProvenance(pseudo?.cutoff_rho_source);
  const wfc = preferredWfc ?? inferWavefunctionCutoff(pseudo, preferredRho);
  const rho = preferredRho ?? inferChargeDensityCutoff(pseudo, preferredWfc);
  return {
    wfc,
    rho,
    wfcStatus: preferredWfc != null
      ? classifyPreferredCutoffStatus(preferredWfcProvenance)
      : wfc != null ? "inferred" : "missing",
    rhoStatus: preferredRho != null
      ? classifyPreferredCutoffStatus(preferredRhoProvenance)
      : rho != null ? "inferred" : "missing",
    wfcProvenance: preferredWfc != null
      ? preferredWfcProvenance
      : wfc != null ? preferredRhoProvenance : "missing",
    rhoProvenance: preferredRho != null
      ? preferredRhoProvenance
      : rho != null ? preferredWfcProvenance : "missing",
    wfcDerivation: preferredWfc != null
      ? "direct"
      : wfc != null ? "from_rho" : "missing",
    rhoDerivation: preferredRho != null
      ? "direct"
      : rho != null ? "from_wfc" : "missing",
    wfcRatio: preferredWfc != null || wfc == null ? null : getPseudoCutoffRatio(pseudo),
    rhoRatio: preferredRho != null || rho == null ? null : getPseudoCutoffRatio(pseudo),
  };
}

export function cutoffStatusRank(status: CutoffStatusKind): number {
  switch (status) {
    case "parsed":
      return 2;
    case "inferred":
      return 1;
    default:
      return 0;
  }
}

export function cutoffProvenanceRank(provenance: CutoffProvenance): number {
  switch (provenance) {
    case "sssp":
      return 6;
    case "djrepo":
      return 5;
    case "upf_info":
      return 4;
    case "upf":
      return 3;
    case "upf_fallback":
      return 2;
    case "unknown":
      return 1;
    default:
      return 0;
  }
}

export function summarizeSelectedPseudoCutoffs(
  elements: string[],
  selectedPseudos: Record<string, string>,
  pseudoMetadataByFilename: Record<string, PseudopotentialMetadata>,
  ssspData: Record<string, SSSPElementData> | null,
): SelectedPseudoCutoffSummary {
  if (elements.length === 0) {
    return {
      maxWfc: 0,
      maxRho: 0,
      wfcStatus: "idle",
      rhoStatus: "idle",
      wfcProvenance: "missing",
      rhoProvenance: "missing",
      wfcDerivation: "missing",
      rhoDerivation: "missing",
      wfcRatio: null,
      rhoRatio: null,
      hasInferredCutoff: false,
      hasMissingCutoff: false,
    };
  }

  let selectedCount = 0;
  let maxWfc = 0;
  let maxRho = 0;
  let maxWfcStatus: CutoffStatusKind = "missing";
  let maxRhoStatus: CutoffStatusKind = "missing";
  let maxWfcProvenance: CutoffProvenance = "missing";
  let maxRhoProvenance: CutoffProvenance = "missing";
  let maxWfcDerivation: CutoffDerivation = "missing";
  let maxRhoDerivation: CutoffDerivation = "missing";
  let maxWfcRatio: number | null = null;
  let maxRhoRatio: number | null = null;
  let hasInferredCutoff = false;
  let hasMissingCutoff = false;

  for (const element of elements) {
    const selectedFilename = selectedPseudos[element];
    if (!selectedFilename) {
      continue;
    }

    selectedCount += 1;
    const selectedPseudo = pseudoMetadataByFilename[selectedFilename];
    const {
      wfc,
      rho,
      wfcStatus,
      rhoStatus,
      wfcProvenance,
      rhoProvenance,
      wfcDerivation,
      rhoDerivation,
      wfcRatio,
      rhoRatio,
    } = resolvePseudoCutoffs(selectedFilename, selectedPseudo, ssspData?.[element] ?? null);
    if (wfc != null) {
      if (
        wfc > maxWfc
        || (
          wfc === maxWfc
          && (
            cutoffStatusRank(wfcStatus) > cutoffStatusRank(maxWfcStatus)
            || (
              cutoffStatusRank(wfcStatus) === cutoffStatusRank(maxWfcStatus)
              && cutoffProvenanceRank(wfcProvenance) > cutoffProvenanceRank(maxWfcProvenance)
            )
          )
        )
      ) {
        maxWfc = wfc;
        maxWfcStatus = wfcStatus;
        maxWfcProvenance = wfcProvenance;
        maxWfcDerivation = wfcDerivation;
        maxWfcRatio = wfcRatio;
      }
    }
    if (rho != null) {
      if (
        rho > maxRho
        || (
          rho === maxRho
          && (
            cutoffStatusRank(rhoStatus) > cutoffStatusRank(maxRhoStatus)
            || (
              cutoffStatusRank(rhoStatus) === cutoffStatusRank(maxRhoStatus)
              && cutoffProvenanceRank(rhoProvenance) > cutoffProvenanceRank(maxRhoProvenance)
            )
          )
        )
      ) {
        maxRho = rho;
        maxRhoStatus = rhoStatus;
        maxRhoProvenance = rhoProvenance;
        maxRhoDerivation = rhoDerivation;
        maxRhoRatio = rhoRatio;
      }
    }

    if (wfcStatus === "inferred" || rhoStatus === "inferred") {
      hasInferredCutoff = true;
    }
    if (wfcStatus === "missing" || rhoStatus === "missing") {
      hasMissingCutoff = true;
    }
  }

  return {
    maxWfc,
    maxRho,
    wfcStatus: selectedCount === 0 ? "idle" : maxWfc > 0 ? maxWfcStatus : "missing",
    rhoStatus: selectedCount === 0 ? "idle" : maxRho > 0 ? maxRhoStatus : "missing",
    wfcProvenance: maxWfc > 0 ? maxWfcProvenance : "missing",
    rhoProvenance: maxRho > 0 ? maxRhoProvenance : "missing",
    wfcDerivation: maxWfc > 0 ? maxWfcDerivation : "missing",
    rhoDerivation: maxRho > 0 ? maxRhoDerivation : "missing",
    wfcRatio: maxWfc > 0 ? maxWfcRatio : null,
    rhoRatio: maxRho > 0 ? maxRhoRatio : null,
    hasInferredCutoff,
    hasMissingCutoff,
  };
}
