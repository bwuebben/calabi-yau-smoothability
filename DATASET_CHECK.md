# Checking the Kreuzer–Skarke dataset copy

`src/check_ks_dataset.py` checks the parquet files identified by
`manifests/ks_polytopes_4d_sha256.tsv`. The manifest pins revision
`60c0e119a03608418df538191f65da3f43b5b819` of
[calabi-yau-data/polytopes-4d](https://huggingface.co/datasets/calabi-yau-data/polytopes-4d/tree/60c0e119a03608418df538191f65da3f43b5b819).
The complete copy consists of 30 files, for 5–33 and 36 vertices,
totalling 15,773,290,651 bytes (15.77 GB).

The checker separates three questions:

1. **File integrity and coverage:** are all expected files present, do their
   byte sizes and SHA-256 digests match the manifest, and do their parquet
   row counts sum to 473,800,776?
2. **Repeated coordinate sets:** after sorting each row's vertices, is any
   coordinate set repeated? This catches vertex reorderings but cannot
   detect equivalent polytopes written in different lattice bases.
3. **Reflexivity and lattice equivalence:** recompute facet equations and
   canonical normal forms using PALP, check dimension four, reflexivity,
   the number of actual vertices and the recorded facet count, and reject
   repeated normal forms. This checks the mathematical representatives,
   independently of whether the supplied coordinates are already in normal
   form.

PALP's `poly.x -fneN` processes batches, without enumerating all lattice
points or calculating Hodge numbers. The facet output distinguishes
reflexive polytopes, and `-N` computes a normal form under `GL(4,Z)`;
see the [PALP documentation](https://hep.itp.tuwien.ac.at/~kreuzer/CY/palpwiki/Documentation_on_poly.x)
and [normal-form reference](https://arxiv.org/abs/1301.6641).
The checker relies on those PALP routines; it is not a separately formalized
verification of PALP.

After **all** 473,800,776 rows pass the lattice check, they are distinct
four-dimensional reflexive polytopes. The
[published Kreuzer–Skarke classification](https://arxiv.org/abs/hep-th/0002240)
then implies that the copy is complete. A checksum pass, coordinate-only
deduplication, or a successful sample does not establish this conclusion.
The script reports `classification_verified: true` only for the complete
lattice check. It does not rerun the papers' singularity classifications,
verify their counts, or validate the supplied Hodge-number columns.

## Commands

Use Python with `pyarrow` installed. The project environment is
`./venv/bin/python`; replace that interpreter path as needed. The lattice
check additionally requires PALP. The script looks on `PATH` and then asks
an installed SageMath for its `poly.x`; alternatively pass
`--palp /path/to/poly.x`.

```bash
# Always verifies bytes against the pinned manifest; lists missing files.
./venv/bin/python src/check_ks_dataset.py inventory

# Check every available row for repeated coordinate sets.
./venv/bin/python src/check_ks_dataset.py check \
  --mode stored --allow-partial --batch-size 8192

# Check every available row for reflexivity and GL(4,Z) uniqueness.
./venv/bin/python src/check_ks_dataset.py check \
  --mode lattice --allow-partial --batch-size 2048

# Measure both modes on deterministic uniform samples of available files.
./venv/bin/python src/check_ks_dataset.py benchmark --sample-per-file 500

# Full mathematical check: refuses to start if any input file is missing.
./venv/bin/python src/check_ks_dataset.py check --mode lattice --batch-size 2048

# Adversarial checks, including a changed lattice basis and invalid inputs.
./venv/bin/python src/test_check_ks_dataset.py
```

Use `--vertices 5 6`, for example, to select available vertex-count files.
A partial check requires `--allow-partial` and never reports full
classification verification. Reports are JSON; use `--report PATH` to
choose a destination. The script never downloads files automatically.

SQLite files under `data/ks-check/` store complete, injectively encoded
normal forms, not just hashes. There is one file per vertex count and mode;
equivalent polytopes necessarily have the same vertex count. Each completed
batch commits its keys and row offset in a single transaction. Repeating
the command resumes it. Input-file, checker-source, or PALP-binary changes
invalidate the checkpoint; use a new `--work-dir` in that case. Do not edit
the checkpoint files. A duplicate, invalid polytope, malformed PALP output,
or timeout produces a failing exit status and report. The batch timeout
defaults to 120 seconds; smaller batches or `--palp-timeout` can be used
for unusually expensive files.

The checks run serially and stream the parquet data in bounded batches.
SQLite uses disk for its index. Larger data and indexes can slow the check;
sample timings are not guaranteed full-run timings. Benchmark sampling time
is reported separately from normal-form and index-insertion time, so the
sample's linear extrapolation omits full sequential decoding costs and
large-index effects. In particular, samples from files with at most nine
vertices do not measure the cost of the larger polytopes that dominate the
complete classification.

## Recorded partial verification

The verification run of 8 September 2026 used only the 5–9 vertex files:

| Vertices | Rows |
|---:|---:|
| 5 | 1,561 |
| 6 | 24,189 |
| 7 | 177,446 |
| 8 | 834,638 |
| 9 | 2,867,955 |
| **Total** | **3,905,789** |

These files occupy 83,732,464 bytes (79.9 MiB), cover 0.82435% of the
classification, and all match the pinned SHA-256 digests. The other 25 files
were outside this run and total 15,689,558,187 bytes (15.69 GB). All tested
coordinate sets were distinct; the stored-coordinate check took 66.99
seconds and used 200.6 MB of SQLite files. This is a partial dataset result.

The full lattice check on those same **3,905,789 rows passed**, finding no
repetitions up to `GL(4,Z)` and no invalid reflexive polytopes, vertex counts,
or facet counts. It took **196.48 seconds (3 minutes 16 seconds)** and used
199.7 MB of SQLite files on an Apple M5. Eleven adversarial tests passed,
including detection of a duplicate in a different lattice basis,
non-reflexive and lower-dimensional input, and checkpoint rollback.
Both modes together used about 400 MB of checkpoint files in this run.

Scaling these complete local runs solely by row count gives **2.26 hours**
for coordinate deduplication or **6.62 hours** for the lattice check across
the full classification. These are provisional extrapolations, **not
measurements of the complete dataset**: this run did not time the 10–36
vertex inputs, whose larger polytopes and disk indexes may take longer.
Allow an overnight run initially, with additional time if those files
benchmark more slowly; rerun the benchmark after downloading them. Time
for any required downloads is additional and depends on the connection.
The complete copy plus one lattice-check index should be budgeted at tens
of GB, with additional disk headroom; current local index sizes cannot
predict its exact size.

[Recorded inputs, per-file results and timings](output/ks_dataset_check_2026-09-08.json)
contain the reproducible evidence and sample benchmark. The report covers
only the stated 3,905,789 rows. The abstracts report the database computations;
the body of each paper retains its database hypothesis until the complete
check has actually passed.
