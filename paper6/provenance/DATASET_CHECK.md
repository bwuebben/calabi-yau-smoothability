# Classification input provenance

The external input is the 30 parquet files for 5–33 and 36 vertices at revision
`60c0e119a03608418df538191f65da3f43b5b819` of
[calabi-yau-data/polytopes-4d](https://huggingface.co/datasets/calabi-yau-data/polytopes-4d/tree/60c0e119a03608418df538191f65da3f43b5b819).
`ks_polytopes_4d_sha256.tsv` pins their SHA-256 digests and byte sizes.
`classification_inventory.json` records all 473,800,776 rows in those files,
totalling 15,773,290,651 bytes. These inventory records are not themselves a
proof of the classification or a new test of the excluded polytopes.

The inventory intentionally has `classification_verified: false`: the
inventory operation verifies input integrity and coverage only. A separate
historical full check of the same pinned corpus tested reflexivity and
uniqueness of PALP normal forms for all rows. That computation and the
published Kreuzer–Skarke classification underlie the paper's statement about
the database copy. The present bundle supplies the singularity filter and
triangulation checks; it does not purport to repeat the original
Kreuzer–Skarke classification.

The original exhaustive admissibility filter processed every row. The
independent implementation covered every row in files 05–09 and 24–33/36,
and every row with d ≤ 7 in files 10–23, together with seeded samples. Here d
is the number of nonzero lattice points that are not vertices. Rows with
d > 7 in files 10–23 remain sampled. The completed extension runs tested
15,637,714 and 18,627,401 rows, with overlap; their sum is not a distinct-row
count. All 51,827 positive candidates are accounted for by row identity,
vertices and face types.

For a bounded reproduction, retrieve `polytopes-4d-05-vertices.parquet` from
the pinned revision, check its size (27,519 bytes) and its digest in the TSV,
and run the command in the bundle README. It must return five admissible
polytopes among 1,561 rows, all matching the positive corpus.

The raw parquet files are deliberately external. A successful positive
certificate check, checksum comparison or sampled filter run does not prove
exhaustiveness of the original admissibility scan.
