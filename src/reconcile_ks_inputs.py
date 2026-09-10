#!/usr/bin/env python3
"""Reconcile saved census inputs with the verified Kreuzer--Skarke inventory.

This checks recorded provenance and row coverage; it does not rerun the
singularity predicates or retroactively attest to historical input bytes.
Run from any directory. Requires only the Python standard library.
"""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REVISION = "60c0e119a03608418df538191f65da3f43b5b819"
REPORT = "ks_dataset_check_2026-09-09.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reconcile():
    report_path = ROOT / "output" / REPORT
    report = json.loads(report_path.read_text())
    assert report["classification_verified"] is True
    assert report["status"] == "full_lattice_check_passed"
    assert report["checker_sha256"] == digest(ROOT / "src/check_ks_dataset.py")
    inventory = {r["name"]: r for r in report["inventory"]["files"]}
    assert len(inventory) == 30
    assert sum(r["rows"] for r in inventory.values()) == 473800776
    manifest = {}
    for line in (ROOT / "manifests/ks_polytopes_4d_sha256.tsv").read_text().splitlines():
        if line and not line.startswith("#"):
            sha, size, name = line.split()
            manifest[name] = {"sha256": sha, "bytes": int(size)}
    assert manifest == {n: {k: r[k] for k in ("sha256", "bytes")}
                        for n, r in inventory.items()}
    upstream = json.loads((ROOT / "manifests/ks_upstream_history.json").read_text())
    assert upstream["history"][0]["id"] == REVISION
    assert upstream["history"][0]["date"] < "2026-07-10"
    assert len(upstream["file_objects"]) == 2
    assert all(objects == manifest for objects in upstream["file_objects"].values())
    result = {
        "verification_report": REPORT,
        "verification_report_sha256": digest(report_path),
        "dataset_revision": REVISION,
        "upstream_data_objects_unchanged_since": "2024-02-18",
        "verified_files": len(inventory),
        "verified_rows": sum(r["rows"] for r in inventory.values()),
        "scope": "Recorded source identity and row coverage; no singularity rescan or retrospective input-byte attestation.",
        "scan_families": {},
    }
    for family, pattern, expected_count in [
        ("local_obstruction", "ks_v*.json", 30),
        ("both_sides", "both_sides_v*.json", 29),
    ]:
        paths = sorted((ROOT / "output").glob(pattern))
        if not paths:  # Each ancillary package contains its own census results.
            continue
        records = {}
        artifacts = []
        for path in paths:
            data = json.loads(path.read_text())
            assert not (records.keys() & data.keys()), path.name
            records.update(data)
            artifacts.append({"name": path.name, "sha256": digest(path)})
        assert len(records) == expected_count
        expected_names = set(inventory)
        if family == "both_sides":
            expected_names.remove("polytopes-4d-36-vertices.parquet")
        assert set(records) == expected_names
        for name, r in records.items():
            assert r["n"] == inventory[name]["rows"], name
            if "input_sha256" in r:
                assert r["input_sha256"] == inventory[name]["sha256"], name
        result["scan_families"][family] = {
            "files": len(records), "rows": sum(r["n"] for r in records.values()),
            "records_with_execution_input_digest": sum("input_sha256" in r for r in records.values()),
            "artifacts": artifacts,
            "separate_36_vertex_check": family == "both_sides",
        }
    groups = [("bk_v0507.json", 5, 7), ("bk_v0809.json", 8, 9),
              ("bk_v10_13.json", 10, 13), ("bk_v14_33.json", 14, 33)]
    if any((ROOT / "output" / name).exists() for name, _, _ in groups):
        records = []
        for name, lo, hi in groups:
            path = ROOT / "output" / name
            r = json.loads(path.read_text())
            expected = sum(v["rows"] for v in inventory.values() if lo <= v["vertices"] <= hi)
            assert r["n"] == expected, name
            records.append({"name": name, "sha256": digest(path), "rows": r["n"]})
        result["scan_families"]["batyrev_kreuzer"] = {
            "groups": records, "rows": sum(r["rows"] for r in records),
            "separate_36_vertex_check": True,
        }
    return result


if __name__ == "__main__":
    print(json.dumps(reconcile(), indent=2, sort_keys=True))
