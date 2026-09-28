# Tropical-cycle computation summaries

`results.json` preserves the numerical evidence reported in the examples
section. `python3 census/check_census.py` checks record consistency and
recomputes the 29 associated manuscript totals. This is a summary check,
not an independent reconstruction of the bases, chain complexes or cycles.

The records have the following provenance and limits:

- `hexagon_free`: 375 per-polytope results, of which 364 satisfy the
  projectivity hypothesis. Each exported result was reconciled with the
  retained detailed result. The reported strict-cycle lattices equal K in
  340 cases and have smaller rank in 24.
- `hexagonal`: 20 A1 results and 71 A2 or mixed-profile results; strict-cycle
  lattices equal K in 19 and 66 cases, respectively.
- `relative_bases`: 367 per-base results, checked against individual retained
  result files. All report that balanced coefficients span K and that the
  relevant networks bound integral relative chains; 191 require more than
  facet networks. These are records of the original computations.
- `relative_hexagonal_bases`: 71 structured per-base records and one additional
  A1 base, the product of two hexagons (list 589). The latter survives as a
  written computation report: 36 generators and 21 lifts bound relative
  chains. Its separate structured result file is unavailable. It is marked
  accordingly in the JSON and must not be counted as a newly replayed result.
- `enumerated_boundary_summary`: the number 327 is retained from the aggregate
  enumeration report. Individual enumeration records are unavailable; the
  complementary number 37 is derived from 364 minus 327.
- `belt_summary`: retained conformance counts for 3,389 belts, including 28
  meeting the discriminant in their interiors. Belt lattices span K in 179
  polytopes, and do so using conforming belts in 171. The excluded eight are
  not silently included in the strict-belt total.

The source hashes identify the preserved input records; they do not replace a
geometric verification. The five full worked examples are provided separately
in `examples/`, where the checker states exactly which assertions it verifies.
No mathematical claim of the paper is strengthened by this summary export.
