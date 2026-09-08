#!/usr/bin/env python3
"""Adversarial checks for the dataset verifier; requires pyarrow and PALP."""
import copy
from pathlib import Path
import tempfile
import unittest

import check_ks_dataset as check


SIMPLEX = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0],
           [0, 0, 0, 1], [-1, -1, -1, -1]]


def row(vertices, facets=5):
    return {"vertices": vertices, "vertex_count": len(vertices), "facet_count": facets}


class DatasetCheckerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.palp = check.find_palp()

    def test_unimodular_change_and_reordering(self):
        original = row(SIMPLEX)
        # Composition of integral shears, permutation, and sign change.
        transformed = row([[-v[3], v[0] + 7*v[1], v[1], v[2] - 3*v[3]]
                           for v in reversed(SIMPLEX)])
        keys = check.palp_keys([original, transformed], 5, self.palp)
        self.assertEqual(keys[0], keys[1])
        # Demonstrates why coordinate dedup alone cannot prove GL uniqueness.
        raw = check.batch_keys([original, transformed], 5, "stored", None, 120)
        self.assertNotEqual(raw[0], raw[1])

    def test_vertex_permutation_is_stored_duplicate(self):
        keys = check.batch_keys([row(SIMPLEX), row(list(reversed(SIMPLEX)))], 5, "stored", None, 120)
        self.assertEqual(keys[0], keys[1])

    def test_distinct_reflexive_polytopes(self):
        weighted = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0],
                    [0, 0, 0, 1], [-1, -1, -1, -2]]
        a, b = check.palp_keys([row(SIMPLEX), row(weighted)], 5, self.palp)
        self.assertNotEqual(a, b)

    def test_scaled_nonreflexive_input_rejected(self):
        with self.assertRaisesRegex(check.CheckError, "reflexivity"):
            check.palp_keys([row([[2*x for x in v] for v in SIMPLEX])], 5, self.palp)

    def test_origin_outside_input_rejected(self):
        with self.assertRaisesRegex(check.CheckError, "reflexivity"):
            check.palp_keys([row([[v[0]+3, *v[1:]] for v in SIMPLEX])], 5, self.palp)

    def test_redundant_nonvertex_rejected(self):
        with self.assertRaisesRegex(check.CheckError, "vertex count"):
            check.palp_keys([row(SIMPLEX + [[0, 0, 0, 0]])], 6, self.palp)

    def test_lower_dimensional_input_rejected(self):
        flat = [[1, 0, 0, 0], [-1, 0, 0, 0], [0, 1, 0, 0],
                [0, -1, 0, 0], [1, 1, 0, 0]]
        with self.assertRaises(check.CheckError):
            check.palp_keys([row(flat)], 5, self.palp)

    def test_metadata_mismatch_rejected(self):
        with self.assertRaisesRegex(check.CheckError, "facet count"):
            check.palp_keys([row(SIMPLEX, facets=6)], 5, self.palp)
        with self.assertRaisesRegex(check.CheckError, "filename"):
            check.validate_row(row(SIMPLEX), 6)

    def test_nonintegral_coordinate_and_repeated_vertex_rejected(self):
        bad = copy.deepcopy(SIMPLEX)
        bad[0][0] = 1.5
        with self.assertRaises(check.CheckError):
            check.validate_row(row(bad), 5)
        with self.assertRaisesRegex(check.CheckError, "repeated vertex"):
            check.validate_row(row(SIMPLEX + [SIMPLEX[0]]), 6)

    def test_checkpoint_resume_and_duplicate_rollback(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "keys.sqlite"
            db, offset = check.open_database(path, {"test": 1})
            self.assertEqual(offset, 0)
            check.insert_keys(db, [b"first", b"second"], 0)
            db.close()
            db, offset = check.open_database(path, {"test": 1})
            self.assertEqual(offset, 2)
            with self.assertRaisesRegex(check.CheckError, "duplicate rows 0 and 3"):
                check.insert_keys(db, [b"third", b"first"], offset)
            self.assertEqual(db.execute("SELECT count(*) FROM keys").fetchone()[0], 2)
            self.assertEqual(db.execute("SELECT next_row FROM state").fetchone()[0], 2)
            db.close()
            with self.assertRaisesRegex(check.CheckError, "identity changed"):
                check.open_database(path, {"test": 2})

    def test_missing_coverage_and_corrupted_download(self):
        import pyarrow as pa
        import pyarrow.parquet as pq

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = root / "polytopes-4d-05-vertices.parquet"
            pq.write_table(pa.Table.from_pylist([row(SIMPLEX)]), data)
            manifest = root / "manifest.tsv"
            lines = []
            for n in [*range(5, 34), 36]:
                digest, size = (check.sha256(data), data.stat().st_size) if n == 5 else ("0"*64, 1)
                lines.append(f"{digest}\t{size}\tpolytopes-4d-{n:02d}-vertices.parquet")
            manifest.write_text("\n".join(lines) + "\n")
            inv = check.inventory(root, manifest)
            self.assertEqual(inv["present_rows"], 1)
            self.assertEqual(len(inv["missing"]), 29)
            self.assertFalse(inv["coverage_complete"])
            self.assertFalse(inv["classification_verified"])
            contents = data.read_bytes()
            data.write_bytes(contents[:-1] + bytes([contents[-1] ^ 1]))
            with self.assertRaisesRegex(check.CheckError, "SHA-256 mismatch"):
                check.inventory(root, manifest)


if __name__ == "__main__":
    unittest.main()
