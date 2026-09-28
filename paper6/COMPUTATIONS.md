# Paper 6 computation bundle

This bundle accompanies *Mirror Symmetry for Conifold Relations in
Calabi–Yau Threefolds*. It contains exact data and checkers for the five
worked examples, forced-node and projectivity calculations, and unimodular
triangulations, together with explicitly identified summaries of the larger
tropical-cycle computations. See `census/README.md` for the distinction between
rechecking geometry and tabulating retained results.

Run commands below from this directory (eventually `paper6/` in the public
repository). Python 3.10+ suffices for the integer, example and forced-node
checks. The rational checker and certificate generator require SageMath;
10.9 was used for the acceptance checks. The classification filter additionally
requires NumPy, PyArrow and PALP's `poly-4d.x`.

## Positive certificates

```sh
python3 unimodular/check_unimodular_certificates.py --table
python3 unimodular/check_unimodular_certificates.py --quick 1 --procs 2
sage -python unimodular/check_certificates_rational.py --sample 1 --procs 2
sage -python unimodular/corrupted_certificates_test.py
```

All 51,827 positive certificates were previously checked in integer and exact
rational arithmetic (12,001,172 tetrahedra). Full rerun commands, for a suitable
machine, are the integer and rational commands without `--quick`/`--sample`.
The accompanying corpus and corrected mathematical checkers reuse the frozen
completed run, rather than claim that today's short checks constitute a full
rerun. File 32 is legitimately empty. Missing input and wholly empty runs fail.

`--table` derives positive counts from the certificates and reads total file
sizes from `provenance/classification_inventory.json`. These latter counts are
provenance, not a new verification of the excluded polytopes. Missing inventory
is an error, never a zero count.

## Worked examples

```sh
python3 examples/examples_check.py
python3 examples/corrupted_examples_test.py
```

The five JSON files include the actual reconstructed fans, polar
triangulations, heights, polyhedral bases and saved cycles. Rays and node
parallelograms use the manuscript's order. In `nodes_in_pentagon_or_hexagon`,
node numbers are one-based; array indices elsewhere are zero-based.
Triangle vertices are flags `[[dimension, cell_index], ...]`, with dimensions
0/1/2/3 for vertices/edges/two-cells/three-cells. Exact rational geometric
coordinates are serialized as strings. Each cycle's boundary leg is
`[edge_index, two_cell_index, integer_coefficient]`.

The checker verifies faces, node circuits, saturated relation lattices,
triangulations and convexity, fan subdivisions, discriminant boundary
networks and their integral coefficients, including equality of the oriented
triangle-chain boundary with the recorded network in the tangent charts. For the designated belt cycles it
also checks the disc, constant field and absence of interior discriminant
vertices. A closed boundary alone does not imply a belt cycle: X88's tenth
cycle has a closed boundary and a more complicated interior.

The public checker does **not** reconstruct the coupled base from the two
fans, or verify all CBM interior conditions for non-belt cycles. Its docstring
states this limit. Those constructions and broader proof checks remain in the
research record; the exporter is not a substitute for final correctness review.
The printed boundary paths, coefficients, cycle coordinates and specified Euler
characteristics are compared with the exported data. Interior intersections
with the discriminant are counted and their circular links checked; this alone
does not establish affine transversality. The X88 tenth-cycle count is 18.

## Tropical-cycle summaries

```sh
python3 census/check_census.py
```

This recomputes the reported totals, including 340 strict realisations among
364 cases and the corroborating 367 and 72 bases, from retained records. It is
not a replay of those geometric constructions. Two entries rely on written or
aggregate reports; `census/README.md` identifies them.

## Forced nodes and projectivity

```sh
python3 forced_nodes/forced_nodes_check.py
python3 forced_nodes/forced_nodes_check.py --sources
```

`list_polytopes.json.gz` contains 613 admissible vertex sets, an explicit
unimodular equivalence for the recorded duplicate, node profiles, positive
projectivity witnesses and Farkas certificates covering negative assertions.
The checker recomputes the faces, integral kernels, forced nodes, witness
inequalities, certificates and profile coverage. It reproduces the manuscript's
612-polytope count, 447 cases in scope and 161 satisfying Corollary B's
no-forced-node condition. Its `stored` fields are generated summaries checked
against these recomputations, not independent proofs.

`source_vertices.json` retains the 670 original labelled input entries (668
distinct vertex sets). `--sources` checks that the 613 exported sets are exactly
the admissible ones among these inputs. This is a finite-list check, not the
full Kreuzer–Skarke classification. Labels identify the committed source files;
source provenance and comparison to earlier per-polytope results are retained
in the private export record.

The optional generator is not needed to verify the results. A fresh generation
may be long. On the research machine:

```sh
sage -python forced_nodes/proj_certificates_generate.py \
  --hints forced_nodes/projectivity_hints.json --out /tmp/forced.json.gz
python3 forced_nodes/forced_nodes_check.py --data /tmp/forced.json.gz --sources
```

The generator checkpoints after each polytope and never labels a partial run
complete. The two large-case hints are previously saved positive witnesses,
translated by vertex coordinates and checked exactly; they do not bypass any
verification. Preserve matching `.prepared.json` and `.checkpoint.json` files
when resuming. Use a fresh output prefix if the input corpus changes.

## Classification filter

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python unimodular/admissible_filter.py /path/to/polytopes-4d-05-vertices.parquet \
  --out /tmp/admissible_v05.jsonl --compare
```

The original helper modules are included in `unimodular/_lib/`; no private
repository imports are required. The raw parquet files are external. The pinned
revision, byte sizes and SHA-256 hashes are in
`provenance/ks_polytopes_4d_sha256.tsv`. Retrieve from that immutable revision
of `calabi-yau-data/polytopes-4d`, then verify each file's size and hash before
running the filter. `provenance/DATASET_CHECK.md` describes the input provenance and verification
scope. NumPy and PyArrow must be installed in the interpreter that runs this
command; a SageMath installation does not necessarily include PyArrow. PALP
may be on `PATH` or supplied by an installed SageMath. The five-vertex file is a
small end-to-end check: 1,561 rows and five admissible polytopes, all matching the
certificate vertex sets and face types.

The original exhaustive filter processed 473,800,776 rows. The completed
independent replay covered all files 05–09 and 24–33/36, all d≤7 rows in files
10–23 and prescribed seeded samples. Here d is the number of nonzero lattice
points that are not vertices. Rows with d>7 in files 10–23 remain sampled.
The B/C extension counts overlap and must not be added as a distinct-row count.
Checking positive certificates alone does not establish classification completeness.

## File integrity

```sh
python3 verify_manifest.py
```

This checks every bundled file against `MANIFEST.json`. It checks file integrity,
not mathematical correctness. In the public repository these instructions are
named `COMPUTATIONS.md`; its existing `README.md` is preserved.
