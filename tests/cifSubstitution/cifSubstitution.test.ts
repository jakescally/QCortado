import assert from "node:assert/strict";
import test from "node:test";
import {
  inspectCifAtomGroups,
  substituteCifElements,
} from "../../src/lib/cifSubstitution";

const eu5In2Sb6 = `
data_175763-ICSD
_chemical_name_common 'Pentaeuropium(II) diindium hexaantimonide'
_chemical_formula_structural 'Eu5 In2 Sb6'
_chemical_formula_sum 'Eu5 In2 Sb6'
_cell_length_a 12.5069(9)
_cell_length_b 14.5985(9)
_cell_length_c 4.6354(3)
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_H-M_alt 'P b a m'
_space_group_IT_number 55
loop_
_space_group_symop_id
_space_group_symop_operation_xyz
1 'x+1/2, -y+1/2, z'
2 '-x+1/2, y+1/2, z'
loop_
_atom_type_symbol
_atom_type_oxidation_number
Eu0+ 0
Sb0+ 0
In0+ 0
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_symmetry_multiplicity
_atom_site_Wyckoff_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_U_iso_or_equiv
_atom_site_occupancy
Eu1 Eu0+ 2 a 0 0 0 0.01219(5) 1
Sb1 Sb0+ 4 g 0.66503(2) 0.67342(2) 0 0.01168(4) 1
In1 In0+ 4 h 0.66978(2) 0.78537(2) 0.5 0.01381(5) 1
loop_
_atom_site_aniso_label
_atom_site_aniso_type_symbol
_atom_site_aniso_U_11
_atom_site_aniso_U_22
_atom_site_aniso_U_33
_atom_site_aniso_U_12
_atom_site_aniso_U_13
_atom_site_aniso_U_23
Eu1 Eu0+ 0.01018(13) 0.01031(8) 0.01608(11) -0.00086(7) 0 0
Sb1 Sb0+ 0.00954(12) 0.01134(8) 0.01416(9) -0.00112(6) 0 0
In1 In0+ 0.01426(14) 0.01443(9) 0.01276(10) -0.00189(8) 0 0
`;

test("inspects atom-site groups from raw CIF rows", () => {
  const groups = inspectCifAtomGroups(eu5In2Sb6);
  assert.deepEqual(groups.map((group) => [group.element, group.count]), [
    ["Eu", 1],
    ["In", 1],
    ["Sb", 1],
  ]);
});

test("substitutes Eu5In2Sb6 In sites to Tl while preserving structure parameters", () => {
  const result = substituteCifElements(eu5In2Sb6, [{ from: "In", to: "Tl" }], {
    sourceFilename: "Eu5In2Sb6.cif",
  });

  assert.match(result.content, /# QCortado generated CIF substitution: In->Tl from Eu5In2Sb6\.cif/);
  assert.match(result.content, /_chemical_formula_sum 'Eu5 Tl2 Sb6'/);
  assert.match(result.content, /^Tl0\+ 0$/m);
  assert.match(result.content, /^Tl1 Tl0\+ 4 h 0\.66978\(2\) 0\.78537\(2\) 0\.5 0\.01381\(5\) 1$/m);
  assert.match(result.content, /^Tl1 Tl0\+ 0\.01426\(14\) 0\.01443\(9\) 0\.01276\(10\) -0\.00189\(8\) 0 0$/m);
  assert.match(result.content, /^_cell_length_a 12\.5069\(9\)$/m);
  assert.match(result.content, /^1 'x\+1\/2, -y\+1\/2, z'$/m);
  assert.equal(result.originalFormula, "Eu5 In2 Sb6");
  assert.equal(result.newFormula, "Eu5 Tl2 Sb6");
  assert.equal(result.suggestedFilename, "Eu5Tl2Sb6.cif");
  assert.deepEqual(result.changedSites.map((site) => [site.label, site.newLabel, site.typeSymbol, site.newTypeSymbol]), [
    ["In1", "Tl1", "In0+", "Tl0+"],
  ]);
});

test("substitutes VESTA-style atom-site loops with type symbol after coordinates", () => {
  const source = `
data_VESTA_phase_1
_chemical_name_common 'Palladium nickel oxide'
_cell_length_a 2.8300(2)
_cell_angle_gamma 120.000000
loop_
_atom_site_label
_atom_site_occupancy
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_adp_type
_atom_site_B_iso_or_equiv
_atom_site_type_symbol
Pd1 1.0 0.000000 0.000000 0.000000 Biso 0.160000 Pd
Ni1 1.0 0.000000 0.000000 0.500000 Biso 0.150000 Ni
O1 1.0 0.000000 0.000000 0.1112(8) Biso 0.340000 O
`;
  const result = substituteCifElements(source, [{ from: "Ni", to: "Co" }]);
  assert.match(result.content, /^Co1 1\.0 0\.000000 0\.000000 0\.500000 Biso 0\.150000 Co$/m);
  assert.match(result.content, /^_cell_angle_gamma 120\.000000$/m);
});

test("uses atom-site labels when type symbols are absent", () => {
  const source = `
data_labels_only
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Si1 0 0 0 1
O1 0.25 0.25 0.25 1
`;
  const result = substituteCifElements(source, [{ from: "Si", to: "Ge" }]);
  assert.match(result.content, /^Ge1 0 0 0 1$/m);
  assert.match(result.content, /^O1 0\.25 0\.25 0\.25 1$/m);
});

test("validates symbols and duplicate atom-type merges", () => {
  assert.throws(
    () => substituteCifElements(eu5In2Sb6, [{ from: "In", to: "Xx" }]),
    /Invalid replacement element symbol/,
  );
  assert.throws(
    () => substituteCifElements(eu5In2Sb6, [{ from: "In", to: "In" }]),
    /different replacement/,
  );

  const duplicateTarget = `
data_duplicate_target
loop_
_atom_type_symbol
_atom_type_oxidation_number
In0+ 0
Tl0+ 0
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
In1 In0+ 0 0 0
Tl1 Tl0+ 0.5 0.5 0.5
`;
  assert.throws(
    () => substituteCifElements(duplicateTarget, [{ from: "In", to: "Tl" }]),
    /duplicate atom type/,
  );
});

test("supports multiple substitutions in formulas and atom sites", () => {
  const result = substituteCifElements(eu5In2Sb6, [
    { from: "Eu", to: "Ba" },
    { from: "Sb", to: "As" },
  ]);
  assert.match(result.content, /_chemical_formula_sum 'Ba5 In2 As6'/);
  assert.match(result.content, /^Ba1 Ba0\+ 2 a 0 0 0 0\.01219\(5\) 1$/m);
  assert.match(result.content, /^As1 As0\+ 4 g 0\.66503\(2\) 0\.67342\(2\) 0 0\.01168\(4\) 1$/m);
});
