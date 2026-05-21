#!/usr/bin/env node
import fs from "node:fs";
import os from "node:os";
import path from "node:path";

const NUMBER_PATTERN = "[-+]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[Ee][-+]?\\d+)?";

function parseArgs(argv) {
  const args = {
    projectsDirs: [],
    write: false,
    backup: true,
  };

  for (let i = 0; i < argv.length; i += 1) {
    const arg = argv[i];
    if (arg === "--projects-dir") {
      const next = argv[++i];
      if (!next) throw new Error("--projects-dir requires a path");
      args.projectsDirs.push(path.resolve(next));
    } else if (arg === "--write") {
      args.write = true;
    } else if (arg === "--dry-run") {
      args.write = false;
    } else if (arg === "--no-backup") {
      args.backup = false;
    } else if (arg === "--help" || arg === "-h") {
      printHelp();
      process.exit(0);
    } else {
      throw new Error(`Unknown argument: ${arg}`);
    }
  }

  if (args.projectsDirs.length === 0) {
    args.projectsDirs.push(
      path.join(os.homedir(), "Library/Application Support/com.qcortado.app/projects"),
      path.join(os.homedir(), "Library/Application Support/com.qcortado.viewer/projects"),
    );
  }

  return args;
}

function printHelp() {
  console.log(`Patch saved magnetic SCF calculations with per-atom magnetic moments.

Usage:
  node scripts/patch-magnetic-scf-moments.mjs [--write] [--projects-dir PATH]

Options:
  --dry-run            Report what would change. This is the default.
  --write              Update calc.json and project.json files in place.
  --projects-dir PATH  Patch one QCortado projects directory. Can be repeated.
  --no-backup          Do not write .bak-* backups before modifying files.
`);
}

function readJson(filePath) {
  return JSON.parse(fs.readFileSync(filePath, "utf8"));
}

function writeJson(filePath, value, options) {
  if (!options.write) return;
  if (options.backup) {
    const backupPath = `${filePath}.bak-${options.backupStamp}`;
    if (!fs.existsSync(backupPath)) {
      fs.copyFileSync(filePath, backupPath);
    }
  }
  fs.writeFileSync(filePath, `${JSON.stringify(value, null, 2)}\n`);
}

function isObject(value) {
  return value !== null && typeof value === "object" && !Array.isArray(value);
}

function isValidStructure(value) {
  return (
    isObject(value) &&
    Array.isArray(value.atoms) &&
    value.atoms.length > 0 &&
    Array.isArray(value.cell_parameters) &&
    value.cell_parameters.length === 3
  );
}

function hasMomentVectors(value) {
  return (
    Array.isArray(value) &&
    value.length > 0 &&
    value.every((entry) => (
      Array.isArray(entry) &&
      entry.length === 3 &&
      entry.every((component) => Number.isFinite(Number(component)))
    ))
  );
}

function isMagneticScf(calc) {
  if (calc?.calc_type !== "scf") return false;
  const params = calc.parameters ?? {};
  const nspin = Number(params.nspin);
  return nspin === 2 || nspin === 4 || params.noncolin === true || params.lspinorb === true;
}

function parseAtomicMagneticMoments(rawOutput) {
  const vectorRegex = new RegExp(
    `\\batom\\s*:?\\s*(\\d+).*?\\bm(?:agn(?:etization)?)?\\s*=\\s*\\(?\\s*(${NUMBER_PATTERN})\\s+(${NUMBER_PATTERN})\\s+(${NUMBER_PATTERN})\\s*\\)?`,
    "gi",
  );
  const atomNumberRegex = /^\s*atom\s+number\s+(\d+)/i;
  const magnetizationRegex = new RegExp(
    `^\\s*magnetization\\s*:\\s*(${NUMBER_PATTERN})\\s+(${NUMBER_PATTERN})\\s+(${NUMBER_PATTERN})(?:\\s|$)`,
    "i",
  );
  const collinearRegex = new RegExp(
    `\\batom\\s*:?\\s*(\\d+).*?\\bmagn\\s*[:=]\\s*(${NUMBER_PATTERN})`,
    "gi",
  );
  const momentsByAtom = new Map();
  let pendingAtomIndex = null;

  for (const line of rawOutput.split(/\r?\n/)) {
    vectorRegex.lastIndex = 0;
    collinearRegex.lastIndex = 0;
    const vectorMatch = vectorRegex.exec(line);
    if (vectorMatch) {
      const atomIndex = Number(vectorMatch[1]) - 1;
      if (atomIndex >= 0) {
        momentsByAtom.set(atomIndex, [
          Number(vectorMatch[2]),
          Number(vectorMatch[3]),
          Number(vectorMatch[4]),
        ]);
      }
      continue;
    }

    const atomNumberMatch = atomNumberRegex.exec(line);
    if (atomNumberMatch) {
      const atomIndex = Number(atomNumberMatch[1]) - 1;
      pendingAtomIndex = atomIndex >= 0 ? atomIndex : null;
      continue;
    }

    const magnetizationMatch = magnetizationRegex.exec(line);
    if (magnetizationMatch && pendingAtomIndex != null) {
      momentsByAtom.set(pendingAtomIndex, [
        Number(magnetizationMatch[1]),
        Number(magnetizationMatch[2]),
        Number(magnetizationMatch[3]),
      ]);
      pendingAtomIndex = null;
      continue;
    }

    const collinearMatch = collinearRegex.exec(line);
    if (collinearMatch) {
      const atomIndex = Number(collinearMatch[1]) - 1;
      if (atomIndex >= 0) {
        momentsByAtom.set(atomIndex, [0, 0, Number(collinearMatch[2])]);
      }
    }
  }

  return Array.from(momentsByAtom.entries())
    .sort(([a], [b]) => a - b)
    .map(([, vector]) => vector);
}

function parseTotalMagnetization(rawOutput) {
  const regex = new RegExp(`total magnetization\\s*=\\s*(${NUMBER_PATTERN})`, "gi");
  let value = null;
  let match;
  while ((match = regex.exec(rawOutput)) !== null) {
    value = Number(match[1]);
  }
  return Number.isFinite(value) ? value : null;
}

function readCalculationOutput(calcDir, calc) {
  const embedded = calc?.result?.raw_output;
  if (typeof embedded === "string" && embedded.trim().length > 0) {
    return embedded;
  }

  for (const filename of ["pw.out", "scf.out"]) {
    const outputPath = path.join(calcDir, filename);
    if (fs.existsSync(outputPath)) {
      return fs.readFileSync(outputPath, "utf8");
    }
  }

  return "";
}

function findSummaryCalculation(project, calcId) {
  for (const variant of project.cif_variants ?? []) {
    const calc = (variant.calculations ?? []).find((entry) => entry.id === calcId);
    if (calc) return calc;
  }
  return null;
}

function applyMagnetismPatch(calc, moments, totalMagnetization) {
  let changed = false;
  calc.result ??= {};

  if (!hasMomentVectors(calc.result.atomic_magnetic_moments)) {
    calc.result.atomic_magnetic_moments = moments;
    changed = true;
  }

  if (
    totalMagnetization !== null &&
    !Number.isFinite(Number(calc.result.total_magnetization))
  ) {
    calc.result.total_magnetization = totalMagnetization;
    changed = true;
  }

  calc.parameters ??= {};
  if (
    !isValidStructure(calc.parameters.magnetism_viewer_structure) &&
    isValidStructure(calc.parameters.source_structure)
  ) {
    calc.parameters.magnetism_viewer_structure = calc.parameters.source_structure;
    changed = true;
  }

  return changed;
}

function patchProjectsDir(projectsDir, options) {
  const stats = {
    projects: 0,
    magneticScfs: 0,
    patched: 0,
    alreadyPatched: 0,
    noOutput: 0,
    noMoments: 0,
    errors: 0,
  };

  if (!fs.existsSync(projectsDir)) {
    return stats;
  }

  for (const projectEntry of fs.readdirSync(projectsDir, { withFileTypes: true })) {
    if (!projectEntry.isDirectory()) continue;
    const projectDir = path.join(projectsDir, projectEntry.name);
    const projectJsonPath = path.join(projectDir, "project.json");
    if (!fs.existsSync(projectJsonPath)) continue;

    stats.projects += 1;
    let project;
    try {
      project = readJson(projectJsonPath);
    } catch (error) {
      stats.errors += 1;
      console.error(`Failed to read ${projectJsonPath}: ${error.message}`);
      continue;
    }

    let projectChanged = false;
    const calculationsDir = path.join(projectDir, "calculations");
    if (!fs.existsSync(calculationsDir)) continue;

    for (const calcEntry of fs.readdirSync(calculationsDir, { withFileTypes: true })) {
      if (!calcEntry.isDirectory()) continue;
      const calcDir = path.join(calculationsDir, calcEntry.name);
      const calcJsonPath = path.join(calcDir, "calc.json");
      if (!fs.existsSync(calcJsonPath)) continue;

      let calc;
      try {
        calc = readJson(calcJsonPath);
      } catch (error) {
        stats.errors += 1;
        console.error(`Failed to read ${calcJsonPath}: ${error.message}`);
        continue;
      }

      if (!isMagneticScf(calc)) continue;
      stats.magneticScfs += 1;

      const output = readCalculationOutput(calcDir, calc);
      if (!output) {
        stats.noOutput += 1;
        continue;
      }

      const moments = parseAtomicMagneticMoments(output);
      if (moments.length === 0) {
        stats.noMoments += 1;
        continue;
      }

      const totalMagnetization = parseTotalMagnetization(output);
      const calcChanged = applyMagnetismPatch(calc, moments, totalMagnetization);
      const summaryCalc = findSummaryCalculation(project, calc.id);
      const summaryChanged = summaryCalc
        ? applyMagnetismPatch(summaryCalc, moments, totalMagnetization)
        : false;

      if (calcChanged || summaryChanged) {
        stats.patched += 1;
        projectChanged ||= summaryChanged;
        writeJson(calcJsonPath, calc, options);
        console.log(`${options.write ? "Patched" : "Would patch"} ${calcJsonPath} (${moments.length} moments)`);
      } else {
        stats.alreadyPatched += 1;
      }
    }

    if (projectChanged) {
      writeJson(projectJsonPath, project, options);
    }
  }

  return stats;
}

function mergeStats(statsList) {
  return statsList.reduce((acc, stats) => {
    for (const [key, value] of Object.entries(stats)) {
      acc[key] = (acc[key] ?? 0) + value;
    }
    return acc;
  }, {});
}

const options = parseArgs(process.argv.slice(2));
options.backupStamp = new Date().toISOString().replace(/[:.]/g, "-");

if (!options.write) {
  console.log("Dry run only. Re-run with --write to update saved calculations.");
}

const allStats = options.projectsDirs.map((projectsDir) => {
  console.log(`Scanning ${projectsDir}`);
  return patchProjectsDir(projectsDir, options);
});

const total = mergeStats(allStats);
console.log(JSON.stringify(total, null, 2));
