# Smoothability of Batyrev Calabi–Yau threefolds

Code, data, and paper sources for a series of three papers on the
smoothability of Calabi–Yau threefolds arising as anticanonical
hypersurfaces in Gorenstein toric Fano fourfolds, and on what
smoothability means for their Batyrev mirrors.

**Author:** Bernd J. Wuebben (wuebben@gmail.com)

## The papers

1. **Non-smoothable Calabi–Yau threefolds from reflexive polytopes**
   (`paper1/`) — explicit compact Calabi–Yau threefolds that admit no
   smoothing: threefolds whose singular points are anticanonical cones over
   the Hirzebruch surface F₁ (down to a *single* such point), and — a new
   phenomenon — threefolds all of whose singular points deform nontrivially
   while no global smoothing exists. Includes a complete census of the
   isolated Gorenstein toric threefold germs with small edge data (217
   classes: rigid / deformable-but-non-smoothable / smoothable) and a sweep
   of **all 473,800,776 reflexive 4-polytopes** of the Kreuzer–Skarke
   classification: **8.27% (39,175,536)** carry a unit-edge non-smoothable
   2-face, so their generic anticanonical hypersurfaces admit no smoothing.

2. **Smoothing Calabi–Yau threefolds in Gorenstein toric Fano fourfolds**
   (`paper2/`) — the positive direction: an explicit binomial deformation
   (a codimension-two transplant of Petracci's homogeneous deformations of
   toric pairs) smooths the points over a type-(S) face whenever the dual
   edge has lattice length ≥ 2; the threshold is sharp (a 7-vertex polytope
   whose hypersurface has a single ordinary double point and no smoothing);
   the ℓ = 1 del Pezzo-cone cases lie outside both Friedman's criterion and
   Gross's theorems and are open. Comparison with the Batyrev–Kreuzer
   all-conifold census: of their 30,241 Namikawa-certified smoothable
   polytopes, at most 3,774 are explained per-face — at least 26,467 are
   cross-face rescues.

3. **Non-smoothability and mirror non-isolatedness of Batyrev Calabi–Yau
   threefolds** (`paper3/`) — the mirror-symmetric meaning: the dual-edge
   length equals the transverse cone multiplicity, so X and its Batyrev
   mirror X° are both isolated-singular exactly on "both-sides unit"
   polytopes; the paper studies which local germs occur there, with a
   census of both-sides polytopes across the entire classification.
   *Draft under active revision as the classification-wide census
   completes.*

## Code (`src/`)

Exact-arithmetic Python (stdlib only for the core; `numpy` + `pyarrow`
for the database scanners). Every quantitative claim in the papers is
produced by one of these scripts, and the anchor examples are asserted on
every run:

- `toric_census.py` — the local classification engine: rigid /
  deformable-but-non-smoothable / smoothable trichotomy for cones over
  unit-edge lattice polygons, with the 217-class census.
- `batyrev_global.py` — reflexive-4-polytope toolkit (facets, 2-faces in
  induced lattices, dual-edge lengths) and the headline example polytopes.
- `hodge_numbers.py` — Batyrev Hodge numbers of the MPCP resolutions.
- `plant_search.py` — planting non-smoothable polygons as 2-faces.
- `ks_sweep.py` — the full Kreuzer–Skarke sweep (fast integer engine,
  selftested per file against the reference path).
- `missing_polytope.py` — the one polytope the per-vertex-count mirror
  omits (the 36-vertex hexagon×hexagon product), identified and verified.
- `paper2_check.py`, `cascade_check.py` — machine checks for paper 2.
- `bk_check.py` — the Batyrev–Kreuzer all-conifold census.
- `both_sides_ks.py`, `both_sides_fast.py`, `both_sides_search.py`,
  `b1_*.py`, `mirror_check.py` — the both-sides-unit scans for paper 3.

## Data (`output/`)

JSON results of every scan, so all counts can be checked without redoing
the compute (the full sweep is ≈ 60 h): per-vertex-count sweep results
(`ks_v*.json`), the Batyrev–Kreuzer census (`bk_*.json`), planting results
(`plant_*.json`), and the both-sides census (`both_sides_*.json`).

The polytope data itself is the Kreuzer–Skarke classification, mirrored as
parquet at
[huggingface.co/datasets/calabi-yau-data/polytopes-4d](https://huggingface.co/datasets/calabi-yau-data/polytopes-4d)
(not redistributed here); the scanners download per-vertex-count files into
`data/ks/`.

## Reproducing

```bash
python3 src/toric_census.py        # local census + self-tests   (~1 s)
python3 src/batyrev_global.py      # example polytopes, asserted  (~1 s)
python3 src/hodge_numbers.py       # Hodge numbers, asserted      (~1 s)
python3 src/missing_polytope.py    # the 36-vertex polytope       (~2 s)
python3 src/paper2_check.py        # paper-2 machine checks       (~2 s)
python3 src/cascade_check.py       # cascade bookkeeping          (~1 s)

# database scans (need: pip install numpy pyarrow; and the parquet files)
./venv/bin/python src/ks_sweep.py data/ks/polytopes-4d-06-vertices.parquet
./venv/bin/python src/bk_check.py data/ks/polytopes-4d-0*-vertices.parquet
./venv/bin/python src/both_sides_fast.py data/ks/polytopes-4d-09-vertices.parquet --procs 8
```

Each paper builds from its directory with `latexmk -pdf main.tex`.

## License

The code (`src/`) and data files (`output/`) are released under the MIT
License (see `LICENSE`). The paper sources and PDFs (`paper*/`) are
© Bernd J. Wuebben; all rights reserved pending journal publication.
