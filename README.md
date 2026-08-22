# Smoothability of Calabi–Yau threefolds

Code, data, and paper sources for a series of four papers on when a
singular Calabi–Yau threefold admits a smoothing. The setting is
Batyrev's: anticanonical hypersurfaces in Gorenstein toric Fano
fourfolds, whose singular points are cones over the two-dimensional
faces of a reflexive 4-polytope. Paper 1 classifies the local germs and
exhibits non-smoothable threefolds; paper 2 constructs smoothings;
paper 3 takes the census of mirror pairs in which both members have
isolated singularities; paper 4 proves a global obstruction, a
vanishing-cycle necessity criterion valid for any projective Calabi–Yau
threefold with nodes and del Pezzo cone points, and uses it to settle
the distinguished example the census singles out.

**Author:** Bernd Johannes Wuebben (wuebben@gmail.com)

| paper | directory | compiled PDF |
|---|---|---|
| 1. Non-smoothable Calabi–Yau threefolds from reflexive polytopes | `paper1/` | `cy-non-smoothable.pdf` |
| 2. Smoothing Calabi–Yau threefolds in Gorenstein toric Fano fourfolds | `paper2/` | `cy-smoothing.pdf` |
| 3. Doubly isolated Batyrev mirror pairs and non-smoothable Calabi–Yau threefolds | `paper3/` | `cy-mirror-pairs.pdf` |
| 4. A vanishing-cycle obstruction to smoothing Calabi–Yau threefolds | `paper4/` | `cy-vanishing-cycle.pdf` |

Each `paperN/` directory holds the LaTeX source (`main.tex`) and the
compiled PDF under the name above. The supporting code lives in `src/`
(papers 1–3 and shared exact toric modules) and
`paper4/certificates/` (paper 4). The three core Paper 4 checkers are
standalone; its general census, fixed-example, and fan programs reuse
modules in `src/`;
scan results in `output/`; the pinned Kreuzer–Skarke input manifest in
`manifests/`. Every paper builds from its directory with
`latexmk -pdf main.tex`, and every quantitative claim is recomputed by
a script listed below.

## The papers

1. **Non-smoothable Calabi–Yau threefolds from reflexive polytopes**
   (`paper1/`) — explicit compact Calabi–Yau threefolds that admit no
   smoothing: threefolds whose singular points are anticanonical cones over
   the Hirzebruch surface F₁ (down to a *single* such point), and — a new
   phenomenon — threefolds all of whose singular points deform nontrivially
   while no global smoothing exists. Includes a complete census of the
   isolated Gorenstein toric threefold germs with small edge data —
   exactly 217 classes, proved and not merely computed, sorted into the
   trichotomy type (R) rigid / type (D) deformable but non-smoothable /
   type (S) smoothable — and a sweep
   of **all 473,800,776 reflexive 4-polytopes** of the Kreuzer–Skarke
   classification: **8.27% (39,175,536)** carry a unit-edge non-smoothable
   2-face, so their generic anticanonical hypersurfaces admit no smoothing.

2. **Smoothing Calabi–Yau threefolds in Gorenstein toric Fano fourfolds**
   (`paper2/`) — the positive direction: an explicit deformation of the
   ambient toric pair, cut out by trinomials in Cox coordinates (a
   codimension-two transplant of Petracci's homogeneous deformations of
   toric pairs), smooths the points over a type-(S) face whenever the dual
   edge has lattice length ≥ 2; the threshold is sharp — no unit-edge
   Batyrev hypersurface whose singular locus is a single ordinary double
   point admits a smoothing, and a 7-vertex polytope realizes this;
   the ℓ = 1 del Pezzo-cone cases lie outside both Friedman's criterion and
   Gross's theorems and are open. Comparison with the Batyrev–Kreuzer
   all-conifold census: their criterion constrains only the dual-length-1
   faces, so of their 30,241 Namikawa-certified smoothable polytopes
   exactly 3,774 — one in eight — are explained per-face, and the
   remaining exactly 26,467 are cross-face rescues.

3. **Doubly isolated Batyrev mirror pairs and non-smoothable Calabi–Yau
   threefolds** (`paper3/`) — the mirror-symmetric capstone: the dual-edge
   length equals the transverse-cone multiplicity, so X and its Batyrev
   mirror X° are both isolated-singular exactly on "both-sides unit"
   polytopes. A complete scan of the classification shows there are
   **exactly 590** such polytopes; every 2-face is a triangle, a zonotope,
   or a reflexive polygon; the only non-smoothable germs occurring are the
   cyclic quotients ⅓(1,1,1) and ⅕(1,1,3) and — exactly once — the
   F₁-cone, giving a **unique mirror pair** (22 and 26 vertices,
   resolution Hodge numbers (20,26) and (26,20)) in which X admits no
   smoothing while every singular point of X° is locally smoothable. The
   germs that deform but admit no smoothing never occur at all, so a
   Calabi–Yau threefold carrying one never has a mirror with isolated
   singularities; at the opposite extreme exactly one polytope, the
   self-dual 24-cell, gives a mirror pair with both members smooth. The
   paper also records two refuted conjectures: the natural low-vertex
   guess holds for all 395,406,329 polytopes with at most 18 vertices and
   fails from 19 on. An exact row reduction of the 26 node relations of the
   mirror X° gives rank 18 with two coloops, so the nodes alone admit no
   relation with every coefficient non-zero: any smoothing of X° would have
   to couple them to its four del Pezzo-cone directions.

4. **A vanishing-cycle obstruction to smoothing Calabi–Yau threefolds**
   (`paper4/`) — the necessity theory for the mixed case: any one-parameter
   smoothing of a Calabi–Yau threefold whose singular points are nodes and
   anticanonical del Pezzo cones (degrees 6 and 7) produces a rational
   homology class on the singularity links that dies on a crepant
   resolution, lies in canonical root subspaces at the cone points, and is
   nonzero at every node, every dP₇ point, and every line-smoothed dP₆
   point. The engine is a vanishing-period identity (the period of the
   holomorphic 3-form over a vanishing cycle equals the local smoothing
   parameter times a unit); in the purely nodal case this recovers the
   necessity half of Friedman's criterion at all orders, without
   unobstructedness. The criterion is a finite linear-algebra test, and it
   obstructs: the X° member of paper 3's unique mirror pair (26 nodes, two
   dP₇ and two dP₆ cone points, every germ locally smoothable) admits **no
   smoothing**, so both members of the pair are non-smoothable, one for
   local and one for global reasons. On a second, 19-vertex example (14
   nodes and one dP₇ point) the test is silent, and its smoothability is
   posed as an open problem.

The classification-wide counts of papers 1–3 (the 8.27%, the
3,774 / 26,467 split, the 590 both-sides unit polytopes and the
uniqueness of the mirror pair) are stated in the papers under an explicit
hypothesis on the database copy scanned: that the per-vertex-count files
contain, without repetition, exactly the classification members in their
stated vertex range, the one 36-vertex member being supplied separately
(`missing_polytope.py`). The transverse identity, the local trichotomy and
census, and every assertion about an explicitly displayed polytope are
unconditional.

## Code (`src/`)

Exact-arithmetic Python (stdlib only for the core; `numpy` + `pyarrow`
for the database scanners). Every quantitative claim in the papers is
produced by one of these scripts, and the anchor examples are asserted on
every run:

- `toric_census.py` — the local classification engine: the type (R) /
  type (D) / type (S) trichotomy for cones over unit-edge lattice
  polygons, with the 217-class census.
- `batyrev_global.py` — reflexive-4-polytope toolkit (facets, 2-faces in
  induced lattices, dual-edge lengths) and the headline example polytopes.
- `hodge_numbers.py` — Batyrev Hodge numbers of the MPCP resolutions.
- `plant_search.py` — planting non-smoothable polygons as 2-faces.
- `ks_sweep.py` — the full Kreuzer–Skarke sweep (fast integer engine,
  selftested per file against the reference path).
- `missing_polytope.py` — the one polytope the per-vertex-count parquet
  files omit (the 36-vertex hexagon×hexagon product), identified and
  verified.
- `paper2_check.py`, `cascade_check.py` — machine checks for paper 2.
- `bk_check.py` — the Batyrev–Kreuzer all-conifold census.
- `both_sides_ks.py`, `both_sides_fast.py`, `both_sides_census.py` (the
  590-polytope census, fully asserted), `both_sides_search.py`,
  `both_sides_chain.sh` (the full-database driver: downloads each input
  from the pinned dataset revision, verifies its digest, validates an
  existing result before skipping it, writes new results atomically and
  records a transcript), `verify_both_sides_artifact.py` (the resumption
  gate: input digest, result schema and row count, and every recorded
  positive hit re-checked exactly), `b1_*.py`, `mirror_check.py` — the
  both-sides-unit scans for paper 3.
- `paper3_node_relations.py` — the node subsystem of the mirror X°: exact
  rank and coloops of its 26 diagonal relations.
- `face_data.py` — the named 2-faces printed in paper 3, each re-derived
  from the vertex list exactly as it appears in the text.

All predicates and invariants — including the cyclic ordering of the
vertices of a planar face, which uses an exact half-plane and
cross-product comparator — are computed in integer or rational arithmetic;
no floating-point operation enters any classification.

### Paper 4's certificates (`paper4/certificates/`)

The certificate programs named in paper 4's Appendix A are runnable in
the repository layout. The first three Python programs below are
standalone; `examples.py`, `fixed_example.py`, and the Sage fan programs
reuse exact modules in `src/`, so the arXiv ancillary layout preserves
both directories:

- `milnor_kernels.py` — the lattice package: the link lattices and root
  subspaces, the (−2)-enumerations with rigorous bounds, and the
  isometries to the standard del Pezzo markings (49 assertions).
- `global_kernel.py` — the matrix package: the 36×26 matrix rebuilt from
  the paper's printed tables; ranks, coloops, kernels, the four profiles
  and the forced-zero lists (27 assertions).
- `dp_periods.py` — the period package: the fan combinatorics and the
  constancy of the pulled-back holomorphic 3-form behind the 4π² period
  (37 assertions).
- `examples.py` — the finite test of the paper's Section 7, run on both
  the 22-vertex and the 19-vertex examples (8 assertions).

On Sage: `resolution_fan.sage`, `restriction_data.sage`,
`mixed_candidate.sage`, `global_section.sage` (loads `cox_data.sage`),
`local_models.sage` and `chart_characters.sage`, with the face-data
module `fixed_example.py`, certify the fixed crepant subdivision, the
restriction lattices, the 36×26 matrix, Δ-regularity, and the local
deformation bases. The `*.md` files in the same directory are the data
documents these programs derive and check.

## Data (`output/`)

JSON results of every scan, so all counts can be checked without redoing
the compute (the full sweep is ≈ 60 h): per-vertex-count sweep results
(`ks_v*.json`), the Batyrev–Kreuzer census (`bk_*.json`), planting results
(`plant_*.json`), and the both-sides census (`both_sides_*.json`).

The polytope data itself is the Kreuzer–Skarke classification, republished
as parquet at
[huggingface.co/datasets/calabi-yau-data/polytopes-4d](https://huggingface.co/datasets/calabi-yau-data/polytopes-4d)
(not redistributed here); the scanners download per-vertex-count files into
`data/ks/`. `manifests/ks_polytopes_4d_sha256.tsv` pins that dataset to the
immutable repository revision `60c0e119a03608418df538191f65da3f43b5b819` and
records the byte size and SHA-256 digest of every per-vertex-count file
(5–33 vertices and the separate 36-vertex file), so the finite input of every
scan is identified exactly; `verify_both_sides_artifact.py` checks a
downloaded file against it.

## Reproducing

```bash
python3 src/toric_census.py        # local census + self-tests   (~1 s)
python3 src/batyrev_global.py      # example polytopes, asserted  (~1 s)
python3 src/hodge_numbers.py       # Hodge numbers, asserted      (~1 s)
python3 src/missing_polytope.py    # the 36-vertex polytope       (~2 s)
python3 src/paper2_check.py        # paper-2 machine checks       (~2 s)
python3 src/cascade_check.py       # cascade bookkeeping          (~1 s)
python3 paper4/certificates/milnor_kernels.py  # paper 4: lattice package (~1 s)
python3 paper4/certificates/global_kernel.py   # paper 4: matrix package  (~1 s)
python3 paper4/certificates/dp_periods.py      # paper 4: period package  (~1 s)
python3 paper4/certificates/examples.py        # paper 4: finite test     (~1 s)

# database scans (need: pip install numpy pyarrow; and the parquet files)
./venv/bin/python src/ks_sweep.py data/ks/polytopes-4d-06-vertices.parquet
./venv/bin/python src/bk_check.py data/ks/polytopes-4d-0*-vertices.parquet
./venv/bin/python src/both_sides_fast.py data/ks/polytopes-4d-09-vertices.parquet --procs 8
./venv/bin/python src/verify_both_sides_artifact.py \
  --input data/ks/polytopes-4d-09-vertices.parquet \
  --manifest manifests/ks_polytopes_4d_sha256.tsv \
  --result output/both_sides_v0809_fast.json
./src/both_sides_chain.sh 8              # the full pinned, verifying, resumable scan
```

Each paper builds from its directory with `latexmk -pdf main.tex`.

## License

The code (`src/`), data files (`output/`) and manifest (`manifests/`) are
released under the MIT License (see `LICENSE`). The paper sources and PDFs (`paper*/`) are
© Bernd Johannes Wuebben; all rights reserved pending journal publication.
