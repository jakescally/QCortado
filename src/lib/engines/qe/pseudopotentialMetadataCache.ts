import { PseudopotentialMetadata } from "../../types";
import { SSSPElementData } from "./pseudopotentialCutoffs";

const PSEUDO_METADATA_CACHE_PREFIX = "qcortado.pseudopotential_metadata_cache.v1";

export type PseudopotentialCacheKind = "local" | "remote";

export interface PseudopotentialCacheScope {
  kind: PseudopotentialCacheKind;
  pseudoDir: string;
  profileId?: string | null;
}

export interface PseudopotentialMetadataCacheRecord {
  filenames: string[];
  signature: string;
  metadata: PseudopotentialMetadata[];
  ssspData: Record<string, SSSPElementData> | null;
  ssspMissing: boolean;
  updatedAt: string;
  checkedAt: string;
}

export interface PseudopotentialInventoryEntry {
  filename: string;
  size_bytes: number;
  modified_at_epoch: number;
}

function normalizePseudoDir(pseudoDir: string): string {
  return pseudoDir.trim().replace(/\/+$/, "");
}

function cacheKey(scope: PseudopotentialCacheScope): string | null {
  const pseudoDir = normalizePseudoDir(scope.pseudoDir);
  if (!pseudoDir) return null;
  const profileSegment = scope.kind === "remote" ? (scope.profileId || "active") : "local";
  return [
    PSEUDO_METADATA_CACHE_PREFIX,
    scope.kind,
    encodeURIComponent(profileSegment),
    encodeURIComponent(pseudoDir),
  ].join(".");
}

export function buildPseudopotentialFilenameSignature(filenames: string[]): string {
  return [...new Set(filenames.map((filename) => filename.trim()).filter(Boolean))]
    .sort((left, right) => left.localeCompare(right))
    .join("\n");
}

export function buildPseudopotentialInventorySignature(entries: PseudopotentialInventoryEntry[]): string {
  return [...entries]
    .sort((left, right) => left.filename.localeCompare(right.filename))
    .map((entry) => `${entry.filename}\t${entry.size_bytes}\t${entry.modified_at_epoch}`)
    .join("\n");
}

export function readPseudopotentialMetadataCache(
  scope: PseudopotentialCacheScope,
): PseudopotentialMetadataCacheRecord | null {
  if (typeof window === "undefined") return null;
  const key = cacheKey(scope);
  if (!key) return null;

  try {
    const raw = window.localStorage.getItem(key);
    if (!raw) return null;
    const parsed = JSON.parse(raw) as Partial<PseudopotentialMetadataCacheRecord> | null;
    if (!parsed || !Array.isArray(parsed.filenames) || typeof parsed.signature !== "string" || !Array.isArray(parsed.metadata)) {
      return null;
    }
    return {
      filenames: parsed.filenames,
      signature: parsed.signature,
      metadata: parsed.metadata,
      ssspData: parsed.ssspData ?? null,
      ssspMissing: parsed.ssspMissing ?? false,
      updatedAt: parsed.updatedAt ?? new Date(0).toISOString(),
      checkedAt: parsed.checkedAt ?? parsed.updatedAt ?? new Date(0).toISOString(),
    };
  } catch {
    return null;
  }
}

export function writePseudopotentialMetadataCache(
  scope: PseudopotentialCacheScope,
  record: PseudopotentialMetadataCacheRecord,
): void {
  if (typeof window === "undefined") return;
  const key = cacheKey(scope);
  if (!key) return;

  try {
    window.localStorage.setItem(key, JSON.stringify(record));
  } catch {
    // Ignore cache failures; direct parsing remains available.
  }
}

export function updateCachedPseudopotentialMetadata(
  scope: PseudopotentialCacheScope,
  filename: string,
  metadata: PseudopotentialMetadata,
): PseudopotentialMetadataCacheRecord | null {
  const current = readPseudopotentialMetadataCache(scope);
  if (!current) return null;

  const nextMetadata = current.metadata.filter((entry) => entry.filename !== filename);
  nextMetadata.push(metadata);
  nextMetadata.sort((left, right) => left.filename.localeCompare(right.filename));
  const nextRecord = {
    ...current,
    metadata: nextMetadata,
    updatedAt: new Date().toISOString(),
    checkedAt: new Date().toISOString(),
  };
  writePseudopotentialMetadataCache(scope, nextRecord);
  return nextRecord;
}

export function isPseudopotentialMetadataCacheFresh(
  record: PseudopotentialMetadataCacheRecord | null,
  ttlMs: number,
): boolean {
  if (!record) return false;
  const checkedAtMs = Date.parse(record.checkedAt);
  if (!Number.isFinite(checkedAtMs)) return false;
  return Date.now() - checkedAtMs < ttlMs;
}
