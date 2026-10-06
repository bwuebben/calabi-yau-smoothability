# Smoothability of Calabi–Yau threefolds

Paper sources, exact computations, and census data for six papers on
smoothings of singular Calabi–Yau threefolds. Papers 1–3 study
anticanonical hypersurfaces in Gorenstein toric Fano fourfolds: local
singularity types, ambient smoothing constructions, and isolated mirror
pairs. Paper 4 develops a global obstruction using vanishing periods.
Paper 5 proves a necessary and sufficient mixed smoothing criterion and
convergent integrability for its specified deformation spaces. Paper 6 proves
the mirror counterpart: the vanishing spheres of the Batyrev mirror satisfy
exactly the relations of the exceptional curves.

Papers 1–3 give explicit examples, local classifications and computer-assisted
database-wide counts using the verified input copy described below. Paper 4
specifies the singularity and resolution hypotheses for its homological
matrix, and Paper 5 gives the precise slice hypotheses and threshold count
in its ambient obstruction argument.

Paper 5 covers nodes and exact anticanonical cones over smooth del Pezzo
surfaces of degrees 5, 6 and 7, and over P¹×P¹. Its selected deformation
bases are smooth under the stated nonzero degree-six projection hypothesis;
this is distinct from the existence of smooth fibres. The degree-by-degree
scope table is in the introduction. The paper also explains why the
full-tangent/link comparison does not extend directly to degrees 1–4: for a
cubic cone those dimensions are 16 and 6. This does not rule out another
mixed criterion for those degrees.

**Author:** Bernd Johannes Wuebben (wuebben@gmail.com)

| paper | guide | compiled PDF | pages |
|---|---|---|---:|
| 1. Non-smoothable Calabi–Yau threefolds from reflexive polytopes | [paper1/](paper1/README.md) | [cy-non-smoothable.pdf](paper1/cy-non-smoothable.pdf) | 18 |
| 2. Deformations of toric pairs and the smoothing of Batyrev Calabi–Yau threefolds | [paper2/](paper2/README.md) | [cy-toric-pairs.pdf](paper2/cy-toric-pairs.pdf) | 18 |
| 3. Doubly isolated Batyrev mirror pairs and non-smoothable Calabi–Yau threefolds | [paper3/](paper3/README.md) | [cy-mirror-pairs.pdf](paper3/cy-mirror-pairs.pdf) | 11 |
| 4. A vanishing-cycle obstruction to smoothing Calabi–Yau threefolds | [paper4/](paper4/README.md) | [cy-vanishing-cycle.pdf](paper4/cy-vanishing-cycle.pdf) | 35 |
| 5. Smoothing Calabi–Yau threefolds with nodes and del Pezzo cone points | [paper5/](paper5/README.md) | [cy-mixed-smoothing.pdf](paper5/cy-mixed-smoothing.pdf) | 61 |
| 6. Mirror Symmetry for Conifold Relations in Calabi–Yau Threefolds | [paper6/](paper6/README.md) | [cy-vanishing-spheres.pdf](paper6/cy-vanishing-spheres.pdf) | 206 |

Each `paperN/` directory holds `main.tex` and its compiled PDF. Paper 5's
source also includes `sections/`. Supporting computations are in `src/`
(papers 1–3 and shared exact toric modules), `paper4/certificates/`,
`paper5/certificates/`, and `paper5/figures/code/`. Saved scan results are
in `output/`; the pinned Kreuzer–Skarke input manifest is in `manifests/`.
The reproduction commands below distinguish standard Python from SageMath.
Every paper builds from its directory with `latexmk -pdf main.tex`.

The full database-copy verification has passed for all **473,800,776 entries**
in the **30 pinned files**, including reflexivity, lattice-equivalence
uniqueness and final checkpoint validation. The
[complete report](output/ks_dataset_check_2026-09-09.json) records the result;
[checking the dataset copy](DATASET_CHECK.md) gives its scope and reproduction
commands. This check validates the input copy; it does not recompute the
papers' singularity-subset counts. Papers 1–3 use it to discharge their
database-completeness and no-repetition hypotheses; Papers 4–5 identify
the verified inputs to their bounded censuses.

## The papers

1. **Non-smoothable Calabi–Yau threefolds from reflexive polytopes**
   ([paper1/](paper1/README.md)) — explicit projective examples with a
   single F₁-cone point or a single locally deformable nonsmoothable
   pentagon-cone point. Gross's earlier compact constructions are
   distinguished from these explicit reflexive-polytope realizations.
   The local census has **217** isolated singularity classes whose primitive
   polygon edge vectors have coordinates in [−2,2] in some lattice basis:
   **8 reduced-rigid**, **79 deformable nonsmoothable**, and **130 smoothable**,
   including the ordinary double point. Reduced-rigid means that the reduced
   miniversal base is a point; it does not assert vanishing of the first-order
   deformation space. The smooth unimodular triangle is excluded.
   The scan finds that **39,175,536** of the **473,800,776**
   classification polytopes carry a nonsmoothable unit-edge two-face,
   approximately **8.27%**. This is a sufficient local obstruction and gives
   a lower bound for global nonsmoothability, not a complete count of it.

2. **Deformations of toric pairs and the smoothing of Batyrev Calabi–Yau
   threefolds** ([paper2/](paper2/README.md)) — an explicit deformation of
   the ambient toric pair, given by trinomials in Cox coordinates, smooths
   the singularities over a chosen smoothable face when its dual edge has
   an interior lattice point. It gives a global smoothing when this is the
   only singular face. Simultaneous smoothing is unconditional for nodes
   and otherwise assumes irreducibility of the semiuniversal deformation
   space. A seven-vertex one-node example has no smoothing, showing the
   dual-edge threshold is sharp. The lone length-one del Pezzo cases
   (degrees 6 and 7, and P¹×P¹) are nonsmoothable by the surface-homology
   injection and the explicitly recalled necessary condition from Paper 5.
   Exactly **3,774** of the **30,241** smoothings in
   the Batyrev–Kreuzer census meet the criterion face by face; the remaining
   **26,467** require relations among exceptional curves from different faces.

3. **Doubly isolated Batyrev mirror pairs and non-smoothable Calabi–Yau
   threefolds** ([paper3/](paper3/README.md)) — exactly **590** reflexive
   four-polytopes have both associated Batyrev
   hypersurfaces with at most isolated singularities. Their two-faces belong
   to **twelve** lattice-isomorphism classes: triangles, zonotopes, and
   reflexive polygons. The only nonsmoothable germs are ⅓(1,1,1), ⅕(1,1,3)
   and the F₁ cone; type-(D) germs do not occur. The F₁ cone occurs in a
   unique mirror pair in this classification, with **22 and 26 vertices**
   and resolution Hodge numbers **(20,26)** and **(26,20)**. Its displayed
   geometry is verified independently of database completeness. One member
   is locally nonsmoothable; all germs of its mirror are smoothable, but
   Paper 4 proves that the mirror has no global smoothing. The self-dual
   **24-cell** is the unique smooth/smooth pair within the classification.
   The transverse-multiplicity identity supplies the elementary geometric
   interpretation of the enumeration predicate.

4. **A vanishing-cycle obstruction to smoothing Calabi–Yau threefolds**
   ([paper4/](paper4/README.md)) — any smoothing with the stated nodes and
   degree-six/seven anticanonical cone germs forces a relation among link
   homology classes on a crepant resolution. The local classes must bound
   in the chosen Milnor fibres; the relation is nonzero at every node,
   degree-seven point and degree-six point using its one-dimensional local
   component. The period argument applies to arbitrary contact order.
   It obstructs the **26-node/four-cone** mirror member in Paper 3 even
   though every germ is smoothable. A version allowing reduced-rigid germs
   proves that two specified nodes on the other member cannot be smoothed.
   The equality between the ambient matrix kernel and the topological
   kernel has its precise isolated-germ and resolution hypotheses; the
   general ambient-divisor spanning assertion is stated separately.
   Among the **77** polytopes in the stated pentagonal framework with at
   most nine vertices, the criterion obstructs **76**. Paper 5 constructs
   the remaining smoothing and determines the smooth 30-dimensional
   deformation base of the separate example X₁₉.

5. **Smoothing Calabi–Yau threefolds with nodes and del Pezzo cone points**
   ([paper5/](paper5/README.md)) — a necessary and sufficient smoothing
   criterion for connected normal projective complex threefolds with trivial dualizing sheaf,
   H¹(O_X) = 0, and only nodes and the exact cone germs specified above.
   For a chosen profile of local deformation branches, a homological
   relation must be nonzero at each rank-one summand, have all three
   degree-six A₂ coefficients nonzero, and have all five conic-pencil
   pairings nonzero at each degree-five point. At a P¹×P¹ cone the
   rank-one class is the difference of the two rulings. These conditions
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
   one singular fixed point.

   Among reflexive 4-polytopes with at most nine vertices in the stated
   framework (only isolated nodes and degree-six/seven cone points,
   **one singular point per singular two-face**), exactly **77** have a
   pentagonal face: X₉ is smoothable and the other **76** are not. The
   larger facet-rigidity census has **12,508** pentagonal faces:
   **9,466** locked, **37** additional unlocked cases with both containing
   lattice facets rigid, and **3,005** with a decomposable containing facet. These are
   counts of facet conditions, not an unrestricted smoothability census.

6. **Mirror Symmetry for Conifold Relations in Calabi–Yau Threefolds**
   ([paper6/](paper6/README.md)) — establishes and strengthens Morrison's
   numerical prediction for the specified class of toric hypersurfaces with
   nodes and degree-six/seven del Pezzo cone points, assuming a projective
   toric nodal partial resolution. For suitable degenerations of the Batyrev
   mirror family, the exceptional curves and mirror vanishing spheres have
   the same integral relation lattice, and the spheres generate a torsion-free
   subgroup. Explicit integral four-chains realise every relation. The nearby
   nodal locus has the predicted codimension, and smoothability of the nodal
   partial resolution is equivalent to a mirror relation with every coefficient
   nonzero; when a smoothing exists, the mirror admits a symplectic conifold
   transition. When no two-face is a hexagon, the same criterion decides the
   smoothability of the hypersurface itself, with its degree-seven cone points. The required triangulations are verified for all **51,827**
   admissible reflexive polytopes. The [computation bundle](paper6/COMPUTATIONS.md) supplies the
   certificates and checkers, five worked examples, forced-node data and
   explicitly labelled census summaries.

The classification-wide counts use a complete set of distinct reflexive
polytopes up to lattice equivalence. PALP verifies these properties for
all 473,800,776 entries; the published Kreuzer–Skarke total then establishes
completeness. The original scan range covered 5–33 vertices; the 36-vertex
file belongs to the full dataset and was checked separately.

The [input-provenance record](DATASET_CHECK.md#connection-to-the-original-scans)
connects the saved scans to the verified inventory: the upstream data objects
have been unchanged since February 2024, and the recorded filenames and row
counts agree. Legacy scan results did not record execution-time input hashes.
The dataset verification does not independently repeat the singularity scans.

## Code (`src/`)

Exact-arithmetic Python (stdlib only for the core; `numpy` + `pyarrow`
for the database scanners). The scripts reconstruct the finite numerical
data and test the named examples; the geometric and analytic arguments are given in the papers:

- `toric_census.py` — the local classification engine: the type (R) /
  type (D) / type (S) trichotomy for cones over unit-edge lattice
  polygons, with the 217-class census.
- `batyrev_global.py` — reflexive-4-polytope toolkit (facets, 2-faces in
  induced lattices, dual-edge lengths) and the headline example polytopes.
- `hodge_numbers.py` — Batyrev Hodge numbers of the MPCP resolutions.
- `plant_search.py` — planting non-smoothable polygons as 2-faces.
- `reconcile_ks_inputs.py` — compare the verified inventory, pinned upstream
  file history and saved scan coverage without rerunning singularity tests.
- `ks_sweep.py` — the full Kreuzer–Skarke sweep (fast integer engine,
  selftested per file against the reference path).
- `missing_polytope.py` — the polytope outside the original scan
  range (the 36-vertex hexagon×hexagon product), identified and
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
- `gate1_admissible.sage` — a diagnostic showing that the candidate
  decompositions in the area search fail Ilten–Vollmert admissibility (D2).
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

Saved JSON scan results allow the aggregate counts to be reconstructed
without repeating the full scan. This reaggregation does not independently
verify every discarded polytope. The files include per-vertex-count sweep results
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
python3 src/reconcile_ks_inputs.py # input provenance and saved row coverage
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
python3 paper5/figures/code/compute01_cubic_scope.py

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
