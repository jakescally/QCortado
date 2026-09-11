import assert from "node:assert/strict";
import test from "node:test";
import { parseLatestHubbardOccupations, parseWien2kHubbardOccupations } from "../../src/lib/hubbardOccupations";

test("parses the last Hubbard occupations section and atom blocks", () => {
  const parsed = parseLatestHubbardOccupations(`
some earlier output
================= HUBBARD OCCUPATIONS ================
     ------------------------ ATOM    1 ------------------------
     Tr[ns(  1)] (up, down, total) =   6.0  4.0 10.0
     eigenvalues:
       1.000  0.000
     occupations, | n_(i1, i2)^(sigma1, sigma2) |:
       0.900  0.100

later output
================= HUBBARD OCCUPATIONS ================
     ------------------------ ATOM    2 ------------------------
     Tr[ns(  2)] (up, down, total) =   1.0  2.0  3.0
     eigenvalues:
       0.250  0.750
     occupations, | n_(i1, i2)^(sigma1, sigma2) |:
       0.500  0.500
     HUBBARD ENERGY =    0.0423  (Ry)
     extra trailing output
`);

  assert.ok(parsed);
  assert.equal(parsed?.atoms.length, 1);
  assert.equal(parsed?.atoms[0].atomIndex, 2);
  assert.match(parsed?.atoms[0].text ?? "", /eigenvalues:/);
  assert.match(parsed?.atoms[0].text ?? "", /occupations,\s+\| n_\(/);
  assert.match(parsed?.atoms[0].text ?? "", /HUBBARD ENERGY =/);
  assert.doesNotMatch(parsed?.atoms[0].text ?? "", /extra trailing output/);
  assert.doesNotMatch(parsed?.atoms[0].text ?? "", /ATOM\s+1/);
});

test("returns null when no Hubbard occupations are present", () => {
  assert.equal(parseLatestHubbardOccupations("plain pw.out text"), null);
});

test("parses WIEN2k spin density matrices and reports their eigenvalues", () => {
  const up = `
  1 atom density matrix
  1  0.000000  0.000000  0.000000 L, Lx,Ly,Lz in global orthogonal system
  8.00000000E-01  0.00000000E+00   1.00000000E-01  0.00000000E+00   0.00000000E+00  0.00000000E+00
  1.00000000E-01  0.00000000E+00   2.00000000E-01  0.00000000E+00   0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00   0.00000000E+00  0.00000000E+00   4.00000000E-01  0.00000000E+00
`;
  const down = `
  1 atom density matrix
  1  0.000000  0.000000  0.000000 L, Lx,Ly,Lz in global orthogonal system
  3.00000000E-01  0.00000000E+00   0.00000000E+00  0.00000000E+00   0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00   1.00000000E-01  0.00000000E+00   0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00   0.00000000E+00  0.00000000E+00   2.00000000E-01  0.00000000E+00
`;

  const parsed = parseWien2kHubbardOccupations([
    { path: "NiO.dmatup", contents: up },
    { path: "NiO.dmatdn", contents: down },
  ]);

  assert.ok(parsed);
  assert.equal(parsed?.atoms.length, 1);
  assert.equal(parsed?.atoms[0].atomIndex, 1);
  assert.match(parsed?.atoms[0].text ?? "", /spin up/);
  assert.match(parsed?.atoms[0].text ?? "", /spin down/);
  assert.match(parsed?.atoms[0].text ?? "", /Tr\[n\] = 1\.400000/);
  assert.match(parsed?.atoms[0].text ?? "", /0\.183772\s+0\.400000\s+0\.816228/);
  assert.match(parsed?.atoms[0].text ?? "", /occupation matrix/);
});

test("returns null when WIEN2k dmat artifacts are absent or malformed", () => {
  assert.equal(parseWien2kHubbardOccupations([{ path: "NiO.scf", contents: "SCF output" }]), null);
  assert.equal(parseWien2kHubbardOccupations([{ path: "NiO.dmatup", contents: "incomplete" }]), null);
});
