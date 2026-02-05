// Symmetry operations for crystal structures - ported from Reciprocal project

import { AtomSite } from "./types";

interface SymmetryComponent {
  xCoeff: number;
  yCoeff: number;
  zCoeff: number;
  offset: number;
}

interface SymmetryOperation {
  x: SymmetryComponent;
  y: SymmetryComponent;
  z: SymmetryComponent;
}

function parseSymmetryComponent(component: string): SymmetryComponent {
  const result: SymmetryComponent = {
    xCoeff: 0,
    yCoeff: 0,
    zCoeff: 0,
    offset: 0,
  };

  let str = component.trim().toLowerCase();
  if (!str) return result;

  let i = 0;
  let currentSign = 1;

  while (i < str.length) {
    const char = str[i];

    if (char === "+") {
      currentSign = 1;
      i++;
    } else if (char === "-") {
      currentSign = -1;
      i++;
    } else if (char === "x" || char === "y" || char === "z") {
      if (char === "x") result.xCoeff = currentSign;
      else if (char === "y") result.yCoeff = currentSign;
      else if (char === "z") result.zCoeff = currentSign;
      currentSign = 1;
      i++;
    } else if (char >= "0" && char <= "9") {
      let numStr = "";
      while (
        i < str.length &&
        ((str[i] >= "0" && str[i] <= "9") || str[i] === "/" || str[i] === ".")
      ) {
        numStr += str[i];
        i++;
      }

      let value: number;
      if (numStr.includes("/")) {
        const [num, denom] = numStr.split("/");
        value = parseFloat(num) / parseFloat(denom);
      } else {
        value = parseFloat(numStr);
      }

      result.offset += currentSign * value;
      currentSign = 1;
    } else {
      i++;
    }
  }

  return result;
}

function parseSymmetryOperation(opString: string): SymmetryOperation {
  const parts = opString.split(",").map((s) => s.trim());

  return {
    x: parseSymmetryComponent(parts[0] || ""),
    y: parseSymmetryComponent(parts[1] || ""),
    z: parseSymmetryComponent(parts[2] || ""),
  };
}

function applySymmetryOperation(
  op: SymmetryOperation,
  x: number,
  y: number,
  z: number
): { x: number; y: number; z: number } {
  return {
    x: op.x.xCoeff * x + op.x.yCoeff * y + op.x.zCoeff * z + op.x.offset,
    y: op.y.xCoeff * x + op.y.yCoeff * y + op.y.zCoeff * z + op.y.offset,
    z: op.z.xCoeff * x + op.z.yCoeff * y + op.z.zCoeff * z + op.z.offset,
  };
}

function wrapCoordinate(val: number): number {
  let result = val % 1;
  if (result < 0) result += 1;
  if (result > 0.9999) result = 0;
  return result;
}

function positionsEqual(
  p1: { x: number; y: number; z: number },
  p2: { x: number; y: number; z: number },
  tolerance: number = 0.001
): boolean {
  const dx = Math.abs(p1.x - p2.x);
  const dy = Math.abs(p1.y - p2.y);
  const dz = Math.abs(p1.z - p2.z);

  const dxPeriodic = Math.min(dx, 1 - dx);
  const dyPeriodic = Math.min(dy, 1 - dy);
  const dzPeriodic = Math.min(dz, 1 - dz);

  return (
    dxPeriodic < tolerance && dyPeriodic < tolerance && dzPeriodic < tolerance
  );
}

export function generateEquivalentPositions(
  atom: AtomSite,
  symmetryOperations: string[]
): AtomSite[] {
  const positions: AtomSite[] = [];
  const uniquePositions: { x: number; y: number; z: number }[] = [];

  if (symmetryOperations.length === 0) {
    return [atom];
  }

  for (const opString of symmetryOperations) {
    const op = parseSymmetryOperation(opString);
    const newPos = applySymmetryOperation(
      op,
      atom.fract_x,
      atom.fract_y,
      atom.fract_z
    );

    const wrappedPos = {
      x: wrapCoordinate(newPos.x),
      y: wrapCoordinate(newPos.y),
      z: wrapCoordinate(newPos.z),
    };

    const isDuplicate = uniquePositions.some((p) =>
      positionsEqual(p, wrappedPos)
    );

    if (!isDuplicate) {
      uniquePositions.push(wrappedPos);
      positions.push({
        ...atom,
        fract_x: wrappedPos.x,
        fract_y: wrappedPos.y,
        fract_z: wrappedPos.z,
        label: `${atom.label}`,
      });
    }
  }

  return positions;
}

function isNearBoundary(
  val: number,
  tolerance: number = 0.001
): { near0: boolean; near1: boolean } {
  return {
    near0: val < tolerance,
    near1: val > 1 - tolerance,
  };
}

function duplicateBoundaryAtoms(atoms: AtomSite[]): AtomSite[] {
  const result: AtomSite[] = [];
  const tolerance = 0.001;

  for (const atom of atoms) {
    const xBoundary = isNearBoundary(atom.fract_x, tolerance);
    const yBoundary = isNearBoundary(atom.fract_y, tolerance);
    const zBoundary = isNearBoundary(atom.fract_z, tolerance);

    const xPositions = [atom.fract_x];
    const yPositions = [atom.fract_y];
    const zPositions = [atom.fract_z];

    if (xBoundary.near0) xPositions.push(1);
    if (xBoundary.near1) xPositions.push(0);
    if (yBoundary.near0) yPositions.push(1);
    if (yBoundary.near1) yPositions.push(0);
    if (zBoundary.near0) zPositions.push(1);
    if (zBoundary.near1) zPositions.push(0);

    for (const x of xPositions) {
      for (const y of yPositions) {
        for (const z of zPositions) {
          result.push({
            ...atom,
            fract_x: x,
            fract_y: y,
            fract_z: z,
          });
        }
      }
    }
  }

  return result;
}

/**
 * Expand all atom sites using symmetry operations and duplicate boundary atoms
 */
export function expandAtomSites(
  atomSites: AtomSite[],
  symmetryOperations: string[]
): AtomSite[] {
  const expandedSites: AtomSite[] = [];

  for (const atom of atomSites) {
    const equivalentPositions = generateEquivalentPositions(
      atom,
      symmetryOperations
    );
    expandedSites.push(...equivalentPositions);
  }

  return duplicateBoundaryAtoms(expandedSites);
}
