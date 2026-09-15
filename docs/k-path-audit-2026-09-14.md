# K-path audit — 14 September 2026

## Confirmed cause in the saved TaP₂ and NbP₂ runs

The saved QE calculations from 12 September were inspected directly, including their `tmp/bands.in` files and stored band distances. Both use `ibrav=0`, have 503 samples, and contain two internal path breaks. They do not use the faulty inferred-ibrav route described below.

QE retains both endpoints of a disconnected jump. An internal `npoints=0` therefore advances the sample index by **one**, although the plotted distance can stay unchanged. QCortado previously advanced the marker index by zero. This leaves subsequent labels one sample early per preceding break; it does not move the calculated energies.

For TaP₂, the saved I₁ marker was at 4.8888 instead of 4.9131, and I was at 5.5327 instead of 5.5822. For NbP₂, I₁ was at 4.8590 instead of 4.8834, and I was at 5.5032 instead of 5.5528. These are the original QE plot-distance units. The correct final sample is index 502; the old final marker referred to index 500. The first marker immediately after a break can conceal the error because adjacent samples share the same plotted distance.

Wien2k's QCortado path expander also retains both endpoints, and its marker code had the same indexing defect.

**Fixed:** both marker implementations now advance by `max(npoints, 1)`. New calculations save the path counts explicitly. Existing band calculations receive a display-time marker correction when loaded from disk, using saved metadata or an explicit QE line-mode input. Correction requires the complete sample count and marker count to match. It preserves labels, energy arrays and distance arrays, and does not itself rewrite the calculation files.

Reopen the saved runs in an app built from this revision to obtain the correction. No electronic-structure recalculation is needed for this marker error. This identifies a concrete contributor to the observed visual shift; it does not establish that every difference from a literature band plot has the same cause.

## Separate QE primitive-basis handoff

The band wizard could infer a positive QE `ibrav` cell, transform the atoms into that cell, but leave the k-points in spglib's different primitive reciprocal basis. “Primitive” alone does not specify a basis: the vector ordering and signs must agree too.

For the supplied Sr₂RuO₄ cell (`a=3.871 Å`, `c=12.702 Å`), the selected M point is `(0.5, 0.5, -0.5)` in the table/spglib basis. Its Cartesian vector is approximately `(0, 0, 0.4946611) Å⁻¹`. Reusing those coefficients with QE `ibrav=7` instead sends `(1.6231427, 0, 0) Å⁻¹`. This changes the path physically, including interior points.

**Fixed:** Bravais inference now exposes its reciprocal-coordinate transformation. The band and Wannier wizards apply it whenever they use the inferred QE cell. With row-wise direct vectors `A_QE = T A_spglib`, reciprocal fractional columns obey `k_QE = T k_spglib`. Sr₂RuO₄ M is consequently `(0.5, 0.5, 0.5)` in the QE basis, preserving its physical vector. Explicit spglib cells retain their original coordinates.

Previously calculated bands affected by this physical path error need recalculation; relabeling cannot repair their energies. The inspected TaP₂/NbP₂ runs are not in that category.

## Validation and scope

- Parsed all 26 CIFs in the Desktop alias's target folder using QCortado's CIF parser and analyzed them using the actual Rust/spglib implementation.
- Added 24 synthetic cells, covering all 14 Bravais families and a total of 20 detected space groups across the combined 50 cases. Included both tetragonal-I branches, all three orthorhombic-F branches, acute/obtuse rhombohedral cells, axis permutations, and monoclinic/triclinic cells.
- Compiled QE 7.5's actual `generate_k_along_lines.f90` in a standalone driver. Compared its expanded coordinates with QCortado's actual Wien2k Rust path expander. Maximum coordinate discrepancy was `2.22 × 10⁻¹⁶` in the supplied fractional coordinates; internal breaks introduce no interpolation across the jump.
- Repeated the 50-case marker audit after the change: zero marker-index discrepancies.
- Checked Cartesian-coordinate preservation for centered QE cells, including segment interiors. Added regressions for fcc, bcc, tI1, Sr₂RuO₄/tI2, oC, oI and oF, plus the explicit-cell fallback. The corrected Sr₂RuO₄ viewer/export discrepancy was approximately `1.33 × 10⁻¹³ Å⁻¹`.
- Preserved the already-correct TaP₂/NbP₂ monoclinic-C conversions; their checked viewer/export discrepancies were below `3 × 10⁻¹² Å⁻¹`.
- Added fixtures from the exact saved TaP₂/NbP₂ paths and plot-distance arrays, plus their recommended paths with four breaks. Tested old-result repair, idempotence, refusal of incomplete data, and reopening saved QE/Wien2k calculations without rewriting their files.
- Passed **139 frontend tests**, **34 application Rust tests selected by `bands`**, TypeScript checking through the production build, and `git diff --check`. Existing monoclinic-C tests cover all five metric branches and acute/obtuse axis settings.

The Rust tests run with `cargo test --offline --manifest-path src-tauri/Cargo.toml --lib bands`; frontend validation uses `npm run test:unit` and `npm run build`. The standalone investigative harness and detailed before/after JSON remain in `/tmp/qcortado-k-audit` for this session. No new DFT calculations or native Wien2k jobs were run.

## Remaining review items, outside these two fixes

The audit is not a certification of every cell-setting convention. Nonstandard/permuted input cells exposed discrepancies between the viewer's standardized basis and the exporter's input-basis assumptions. Optimized source cells also warrant an explicit comparison against the exact saved calculation basis. Those routes were not changed here.

The primitive-monoclinic table has a separate formula/setting problem: its eta expression uses `b/a` where Setyawan–Curtarolo Table 16 requires `b/c`, and that table assumes acute unique-a axes. The current dispatcher passes beta-based parameters without the normalization used by the monoclinic-C implementation. This does not affect TaP₂/NbP₂, which use monoclinic C, and was left unchanged.

Some triclinic table points lie outside the first Wigner–Seitz cell for tested metrics. That alone does **not** prove a physical k-point error: reciprocal translations can give equivalent points, while paths and labels depend on the chosen convention. Likewise, the fcc Cartesian differences seen before the QE basis fix were consistent with a global cubic symmetry operation in the tested high-symmetry cubic cases; they were not counted as evidence of incorrect energies.

The viewer's current selection code retains coordinates and zero-count breaks through selection/recommended-path construction. It renders fractional points in its primitive reciprocal basis. Its pre-existing local edit was preserved. No interactive mouse/3D rendering test was performed. Native Wien2k interpretation after `sgroup` changes to the case structure remains an integration check, beyond the verified exporter sampling and marker behavior.

## References

- [QE input conventions](https://www.quantum-espresso.org/Doc/INPUT_PW.html): `crystal_b` coordinates belong to the actual calculation cell.
- Local QE 7.5 sources: `Modules/generate_k_along_lines.f90` (retains points at zero-count breaks), `Modules/latgen.f90` (positive-ibrav vectors), and `PP/src/bands.f90` (plot distances). The latter also uses a jump-length heuristic to collapse disconnected segments; no change to that behavior is made here.
- [Setyawan and Curtarolo, high-throughput band structures](https://arxiv.org/abs/1004.2974): primitive-basis definitions and monoclinic Table 16.
