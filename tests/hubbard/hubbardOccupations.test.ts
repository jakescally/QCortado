import assert from "node:assert/strict";
import test from "node:test";
import { parseLatestHubbardOccupations } from "../../src/lib/hubbardOccupations";

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
