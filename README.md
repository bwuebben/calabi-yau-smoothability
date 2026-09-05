# Smoothability of Calabi–Yau threefolds

Paper sources, exact computations, and census data for five papers on
smoothings of singular Calabi–Yau threefolds. Papers 1–3 study
anticanonical hypersurfaces in Gorenstein toric Fano fourfolds: local
singularity types, ambient smoothing constructions, and isolated mirror
pairs. Paper 4 develops a global obstruction using vanishing periods.
Paper 5 proves a necessary and sufficient mixed smoothing criterion and
convergent integrability for its specified deformation spaces.

**Updated 5 September 2026.** Paper 5 now covers nodes and exact
anticanonical cones over smooth del Pezzo surfaces of degrees 5, 6 and 7,
and over P¹×P¹. Integrability is proved under the hypotheses stated in
its theorems. The distinction between smooth deformation bases and the
existence of smooth fibres is essential: the paper also constructs a
nonsmoothable example with four smooth deformation components.
The degree-by-degree scope table is in the introduction.

**Author:** Bernd Johannes Wuebben (wuebben@gmail.com)

| paper | directory | compiled PDF |
|---|---|---|
| 1. Non-smoothable Calabi–Yau threefolds from reflexive polytopes | `paper1/` | [cy-non-smoothable.pdf](paper1/cy-non-smoothable.pdf) |
| 2. Deformations of toric pairs and the smoothing of Batyrev Calabi–Yau threefolds | `paper2/` | [cy-toric-pairs.pdf](paper2/cy-toric-pairs.pdf) |
| 3. Doubly isolated Batyrev mirror pairs and non-smoothable Calabi–Yau threefolds | `paper3/` | [cy-mirror-pairs.pdf](paper3/cy-mirror-pairs.pdf) |
| 4. A vanishing-cycle obstruction to smoothing Calabi–Yau threefolds | `paper4/` | [cy-vanishing-cycle.pdf](paper4/cy-vanishing-cycle.pdf) |
| 5. Smoothing Calabi–Yau threefolds with nodes and del Pezzo cone points | `paper5/` | [cy-mixed-smoothing.pdf](paper5/cy-mixed-smoothing.pdf) |

Each `paperN/` directory holds `main.tex` and its compiled PDF. Paper 5's
source also includes `sections/`. Supporting computations are in `src/`
(papers 1–3 and shared exact toric modules), `paper4/certificates/`,
`paper5/certificates/`, and `paper5/figures/code/`. Saved scan results are
in `output/`; the pinned Kreuzer–Skarke input manifest is in `manifests/`.
The reproduction commands below distinguish standard Python from SageMath.
Every paper builds from its directory with `latexmk -pdf main.tex`.

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

2. **Deformations of toric pairs and the smoothing of Batyrev Calabi–Yau
   threefolds** (`paper2/`) — the positive direction: an explicit deformation of the
   ambient toric pair, cut out by trinomials in Cox coordinates (a
   codimension-two transplant of Petracci's homogeneous deformations of
   toric pairs), smooths the points over a type-(S) face whenever the dual
   edge has lattice length ≥ 2; the threshold is sharp — no unit-edge
   Batyrev hypersurface whose singular locus is a single ordinary double
   point admits a smoothing, and a 7-vertex polytope realizes this;
   the single ℓ = 1 cone cases covered by Paper 5 (degree 6, degree 7,
   and P¹×P¹) are also nonsmoothable, because their relation kernel is zero.
   The general simultaneous ambient construction retains its stated
   irreducibility hypothesis on the deformation space. Comparison with the
   Batyrev–Kreuzer
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
   relation with every coefficient non-zero. Paper 4 further proves that
   coupling them to the four del Pezzo-cone directions cannot give a
   smoothing either.

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
   obstructs the X° member of paper 3's unique mirror pair (26 nodes, two
   dP₇ and two dP₆ cone points, every germ locally smoothable), which admits
   **no smoothing**. A version allowing rigid germs shows independently
   that no deformation of the mirror partner smooths either of two specified
   nodes; its rigid F₁-cone point is not the source of that global
   obstruction. The relation rows give a rational form of each toric germ's
   T¹, the ambient matrix kernel equals the topological kernel on the
   admissible class, and positive cone-point block rank certifies failure of
   ℚ-factoriality. The seven-vertex bottom of the pentagonal census has one
   dP₇-cone point and no smoothing. Among all 77 admissible polytopes with at
   most nine vertices and a pentagonal face, the criterion obstructs 76. On a
   separate hypersurface X₁₉ with 14 nodes and one dP₇ point the test is
   silent; Paper 5 now proves that it is smoothable and that its full
   deformation base is smooth of dimension 30.

5. **Smoothing Calabi–Yau threefolds with nodes and del Pezzo cone points**
   (`paper5/`) — a necessary and sufficient smoothing criterion for connected
   normal projective complex threefolds with trivial dualizing sheaf,
   H¹(O_X) = 0, and only nodes and the exact cone germs specified above.
   For a chosen profile of local deformation branches, a homological
   relation must be nonzero at each rank-one summand, have all three
   degree-six A₂ coefficients nonzero, and have all five conic-pencil
   pairings nonzero at each degree-five point. At a P¹×P¹ cone the
   coefficient is the difference of the two rulings. These conditions
   characterize actual analytic smoothings. The selected deformation
   base is smooth when the relation space projects nontrivially to every
   degree-six summand, and every tangent on that base integrates to a
   convergent curve. No additional projection hypothesis is imposed at
   the other covered germs.

   A threefold with **30 nodes and two degree-six cone points** disproves
   the weaker condition that merely asks for a nonzero vector at each
   germ. Its full analytic deformation ring is
   C{z₁,…,z₃₂,x,y,u,v,b}/(xy,xb,uv,ub), with four smooth components of
   dimensions 34, 34, 34 and 35, but no smoothing.

   The toric applications distinguish intrinsic smoothability from the
   specified ambient constructions. **X₁₉ is smoothable**, with full
   deformation base smooth of dimension 30, although each primitive
   single-degree ambient construction considered in the paper retains a
   singular curve. An explicit ambient family smooths **X₉**, whose two
   nodes lie in no purely nodal relation. Its general ambient fibre is
   toric: the paper prints a fan with seven rays, ten maximal cones and
   one singular fixed point. This corrects the former non-toricity claim.

   Among reflexive 4-polytopes with at most nine vertices in the stated
   framework (only isolated nodes and degree-six/seven cone points,
   **one singular point per singular two-face**), exactly **77** have a
   pentagonal face: X₉ is smoothable and the other **76** are not. The
   larger facet-rigidity census has **12,508** pentagonal faces:
   **9,466** locked, **37** additional unlocked cases with both containing
   lattice facets rigid, and **3,005** with a decomposable containing facet. These are
   counts of facet conditions, not an unrestricted smoothability census.

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
- `mirror_partner.py` — the mirror-partner package: its face inventory,
  crepant resolution data, reduced-rigidity test, relation matrix, and the
  two forced node coordinates (31 assertions).
- `examples.py` — the finite test of Section 9, including the 22-vertex,
  19-vertex, and 20-vertex examples (13 assertions).
- `fixed_example.py` — an independent reconstruction of the distinguished
  mirror pair and its node relations.
- `relation_class.py` — the ray-relation interpretation of the local rows,
  their rational T¹ forms, the block ranks, and an independent reconstruction
  of the rank-21 matrix (36 checks).
- `kernel_equality.py` — the vanishing of the Batyrev correction terms on the
  named examples and all 77 census members (26 checks).
- `lone_germ.py` — the seven-vertex one-germ example (8 checks).
- `classification.py` — the 77-member census, rebuilt from
  `framework_77.json` (10 checks).

On Sage: `resolution_fan.sage`, `restriction_data.sage`,
`mixed_candidate.sage`, `global_section.sage` (loads `cox_data.sage`),
`local_models.sage`, `mirror_partner_fan.sage`, and `chart_characters.sage`, with the face-data
module `fixed_example.py`, certify the fixed crepant subdivision, the
restriction lattices, the 36×26 matrix, Δ-regularity, and the local
deformation bases. The `*.md` files in the same directory are the data
documents these programs derive and check.

### Paper 5's certificates (`paper5/certificates/`)

The toric programs in Appendix B reuse exact modules in `src/` and lattice
data in `paper4/certificates/`, so both directories must be present. The
additional local and deformation-space computations are listed below.

- `v09_candidate.py` — the nine-vertex polytope Δ₉: reflexivity, its face
  inventory, and that it lies in the admissible framework.
- `sing_locus.py` — the singular locus of X₉ is exactly the three germ points,
  with both provisos of the Batyrev statement checked rather than assumed
  (6 checks).
- `locking.py` — the forcing rules and the locking closure,
  including the cube counterexample that shows why the adjacent-pair rule needs
  consecutive edges (13 checks).
- `one_facet.py` — the reachability lemma and its scope (7 checks).
- `slice_rigidity.py` and `slice_rigidity.sage` — cell rigidity by the forcing
  chain, implemented twice, in pure Python with an integer chain and in Sage
  with convex hulls over ℚ and ℚ(s) (49 and 24 checks).
- `gate1_admissible.sage` — a diagnostic showing why the retired area search
  did not certify Theorem A: its canonical candidate decompositions fail
  Ilten–Vollmert admissibility (D2).
- `def41_check.sage` — the passage from the cell to the germ (16 checks).
- `global_decomp.py` — Δ₁₉ and Δ₂₀ are Minkowski-indecomposable, so a global
  decomposition is not available either (12 checks).
- `refine.py` — why one cannot pass to a simplicial ambient (5 checks).
- `dp7_facet_sweep.py` — the census: the locking verdict for every pentagonal
  2-face of every reflexive 4-polytope with at most nine vertices. Results in
  `dp7_sweep_all.json` (the 3,005 decomposing candidates) and
  `dp7_framework_77.json` (the 77 framework pentagons).
- `slice_admissible.sage`, `general_fibre.sage` — the full two-sided
  Ilten–Vollmert axioms, including the arbitrary-collection face condition,
  and completeness of the general fibre (36 and 6 checks).
- `hzero.sage` — an independent invariant-divisor calculation of the
  anticanonical sections.
- `charts.sage` — smoothness of the proper faces and lower-dimensional
  full chart cones (3 checks).
- `toric_fibre.sage` — the complete marked slices of the general ambient
  fibre, its toric fan, its unique singular fixed point, and its 162
  anticanonical lattice points (20 checks).
- `final_check.sage` — the relative canonical and local smoothing
  comparisons (3 checks).
- `branch_locus.py`, `smoothing_locus.py` — the smoothing locus in T¹ agrees
  with paper 4's branch subspace at both germ types, and why the cyclic
  criterion does not apply here; the first script also checks the local Hodge
  and equivariance identifications used in the global comparison (17 and
  13 checks).
- `theoremB_class.py` — the relation class of the explicit X₉ family,
  computed from its local deformation parameters (10 checks).
- `candidates.sage` — what the census survivors induce on their pentagons.

### Paper 5's local and deformation-space computations

The six principal programs verify **507 exact assertions**:

| Program in `paper5/certificates/` | Finite calculation | Assertions |
|---|---|---:|
| `a2_simultaneous_partial_resolution.py` | Incidence charts, nodal Hessians and inverse maps | 33 |
| `dp6_nodal_partial_resolutions.py` | Degree-six subdivisions, root classes and parameter maps | 43 |
| `a2_branch_selection_counterexample.py` | Counterexample polytope, divisor identities and four relation kernels | 272 |
| `branch_smoothness_inputs.py` | Hodge data, local projections and nodal replacements | 52 |
| `counterexample_deformation_germ.py` | Component tangent ideals and their intersection | 57 |
| `a2_counterexample_rational_replay.py` | Independent rational linear algebra and Hodge calculation | 50 |

The first five run with `sage -python`; the rational replay uses Python's
standard library. Their saved JSON results are included, so the replay
has its required input immediately after cloning.

`d19_quadric_partial_resolution.py` additionally checks the X₁₉ fan,
quadric-cone cohomology, Hodge data and relation matrices. The two standard
Python scripts in `paper5/figures/code/` check the degree-five pencil
lattice, permutations, intersection numbers and exterior-product ranks,
and the degree-eight Cayley and ruling-difference lattices. The analytic
lifting and period arguments are proved in the paper.

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
python3 paper4/certificates/mirror_partner.py  # paper 4: mirror partner  (~1 s)
python3 paper4/certificates/examples.py        # paper 4: finite test     (~1 s)
python3 paper4/certificates/fixed_example.py   # paper 4: fixed geometry  (~2 s)
python3 paper4/certificates/relation_class.py  # paper 4: ray relations   (~1 s)
python3 paper4/certificates/kernel_equality.py # paper 4: kernel equality (~2 s)
python3 paper4/certificates/lone_germ.py        # paper 4: one-germ case   (~2 s)
python3 paper4/certificates/classification.py  # paper 4: 77-member census (~3 s)
python3 paper5/certificates/sing_locus.py      # paper 5: Sing(X_9)       (~1 s)
python3 paper5/certificates/locking.py         # paper 5: the closure rule (~1 s)
python3 paper5/certificates/slice_rigidity.py  # paper 5: cell rigidity   (~2 s)
python3 paper5/certificates/global_decomp.py   # paper 5: indecomposability (~2 s)
python3 paper5/certificates/branch_locus.py    # paper 5: the smoothing locus (~2 s)
python3 paper5/certificates/theoremB_class.py  # paper 5: the relation class (~1 s)
sage paper5/certificates/slice_admissible.sage # paper 5: full two-sided axioms
sage paper5/certificates/toric_fibre.sage      # general X9 fan and marked slices
sage paper5/certificates/charts.sage           # chart smoothness
sage -python paper5/certificates/a2_simultaneous_partial_resolution.py
sage -python paper5/certificates/dp6_nodal_partial_resolutions.py
sage -python paper5/certificates/a2_branch_selection_counterexample.py
sage -python paper5/certificates/branch_smoothness_inputs.py
sage -python paper5/certificates/counterexample_deformation_germ.py
python3 paper5/certificates/a2_counterexample_rational_replay.py
sage -python paper5/certificates/d19_quadric_partial_resolution.py
python3 paper5/figures/code/compute01_degree5_local.py
python3 paper5/figures/code/compute01_degree8_local.py

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
