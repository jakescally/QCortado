import {
  CenteringType, Matrix3x3, conventionalToPrimitive, realSpaceLatticeVectors,
} from "./reciprocalLattice";
import type { CrystalData } from "./types";
import type { SymmetryTransformResult } from "./symmetryTransform";

/** Resolve the actual bases in which the viewer's table coordinates are expressed. */
export function getBrillouinZoneLattices(
  crystal: CrystalData,
  centering: CenteringType,
  isMonoclinicC: boolean,
  symmetry: SymmetryTransformResult | null,
): { conventionalLattice: Matrix3x3; primitiveLattice: Matrix3x3 } {
  // mC tables normalize the input setting internally and return coordinates
  // in its canonical C primitive basis. Retain that input frame even when
  // spglib chooses different conventional axes or a different primitive basis.
  if (symmetry && !isMonoclinicC) {
    return {
      conventionalLattice: symmetry.standardizedConventionalLattice,
      primitiveLattice: symmetry.standardizedPrimitiveLattice,
    };
  }
  const conventionalLattice = realSpaceLatticeVectors(
    crystal.cell_length_a.value, crystal.cell_length_b.value, crystal.cell_length_c.value,
    crystal.cell_angle_alpha.value, crystal.cell_angle_beta.value, crystal.cell_angle_gamma.value,
  );
  return {
    conventionalLattice,
    primitiveLattice: conventionalToPrimitive(conventionalLattice, centering),
  };
}
