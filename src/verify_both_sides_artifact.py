#!/usr/bin/env python3
"""Verify a pinned KS input and, optionally, its both-sides result artifact.

This is the resumability gate for ``both_sides_chain.sh``.  It checks the
input byte size and SHA-256 digest against the committed manifest.  For a
saved result it additionally checks JSON structure, the full parquet row
count, aggregate consistency, and every recorded positive hit with the exact
reference toolkit.  The latter does not independently certify negative rows;
that separate limitation is stated explicitly in paper 3.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path


PINNED_REVISION = "60c0e119a03608418df538191f65da3f43b5b819"


class VerificationError(RuntimeError):
    pass


def read_manifest(path: Path) -> dict[str, tuple[str, int]]:
    entries: dict[str, tuple[str, int]] = {}
    for lineno, raw in enumerate(path.read_text().splitlines(), 1):
        if not raw or raw.startswith("#"):
            continue
        fields = raw.split("\t")
        if len(fields) != 3:
            raise VerificationError(f"{path}:{lineno}: expected three TSV fields")
        digest, size, name = fields
        if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            raise VerificationError(f"{path}:{lineno}: malformed SHA-256")
        entries[name] = (digest, int(size))
    return entries


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_input(path: Path, manifest: Path) -> str:
    if not path.is_file():
        raise VerificationError(f"missing input: {path}")
    entries = read_manifest(manifest)
    try:
        expected_digest, expected_size = entries[path.name]
    except KeyError as exc:
        raise VerificationError(f"{path.name} is absent from {manifest}") from exc
    actual_size = path.stat().st_size
    if actual_size != expected_size:
        raise VerificationError(
            f"{path}: size {actual_size}, expected {expected_size}")
    actual_digest = sha256(path)
    if actual_digest != expected_digest:
        raise VerificationError(
            f"{path}: SHA-256 {actual_digest}, expected {expected_digest}")
    return actual_digest


def jsonable(value):
    if isinstance(value, dict):
        return {key: jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    return value


def verify_result(input_path: Path, result_path: Path, input_digest: str) -> None:
    try:
        result = json.loads(result_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise VerificationError(f"cannot parse {result_path}: {exc}") from exc
    if not isinstance(result, dict) or input_path.name not in result:
        raise VerificationError(
            f"{result_path}: missing record for {input_path.name}")
    record = result[input_path.name]
    required = {"n", "both_sides_any", "asymmetric_hits", "hits", "secs"}
    missing = required - set(record) if isinstance(record, dict) else required
    if missing:
        raise VerificationError(f"{result_path}: missing fields {sorted(missing)}")

    import pyarrow.parquet as pq

    rows = pq.ParquetFile(input_path).metadata.num_rows
    if record["n"] != rows:
        raise VerificationError(
            f"{result_path}: n={record['n']}, parquet rows={rows}")
    hits = record["hits"]
    if not isinstance(hits, list):
        raise VerificationError(f"{result_path}: hits is not a list")
    if record["both_sides_any"] != len(hits):
        raise VerificationError(
            f"{result_path}: both_sides_any does not equal len(hits)")
    asymmetric = sum(bool(hit.get("asymmetric")) for hit in hits
                     if isinstance(hit, dict))
    if asymmetric != record["asymmetric_hits"]:
        raise VerificationError(
            f"{result_path}: asymmetric_hits={record['asymmetric_hits']}, "
            f"recount={asymmetric}")
    if not isinstance(record["secs"], (int, float)) or record["secs"] < 0:
        raise VerificationError(f"{result_path}: invalid elapsed time")
    if "input_sha256" in record and record["input_sha256"] != input_digest:
        raise VerificationError(f"{result_path}: embedded input digest mismatch")
    if "input_bytes" in record and record["input_bytes"] != input_path.stat().st_size:
        raise VerificationError(f"{result_path}: embedded input size mismatch")
    if (record.get("dataset_revision") is not None
            and record["dataset_revision"] != PINNED_REVISION):
        raise VerificationError(f"{result_path}: embedded dataset revision mismatch")
    if record.get("complete") is False:
        raise VerificationError(f"{result_path}: segmented result is not complete")
    if "row_start" in record and record["row_start"] != 0:
        raise VerificationError(f"{result_path}: result does not start at row zero")
    if "row_stop" in record and record["row_stop"] != rows:
        raise VerificationError(f"{result_path}: result does not end at the final row")

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from both_sides_fast import verify_hit

    expected_vertices = int(input_path.name.split("-")[2])
    for index, saved in enumerate(hits):
        if not isinstance(saved, dict) or not {"verts", "facekinds", "asymmetric"} <= set(saved):
            raise VerificationError(f"{result_path}: malformed hit {index}")
        vertices = [tuple(int(x) for x in vertex) for vertex in saved["verts"]]
        if len(vertices) != expected_vertices or any(len(vertex) != 4 for vertex in vertices):
            raise VerificationError(f"{result_path}: malformed vertices in hit {index}")
        fresh = jsonable(verify_hit(vertices))
        if fresh != saved:
            raise VerificationError(
                f"{result_path}: exact re-verification differs at hit {index}")

    print(f"verified result: {result_path} ({rows:,} rows, {len(hits)} exact hits)")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--result", type=Path)
    args = parser.parse_args()
    try:
        digest = verify_input(args.input, args.manifest)
        print(f"verified input: {args.input} sha256={digest}")
        if args.result is not None:
            verify_result(args.input, args.result, digest)
    except VerificationError as exc:
        print(f"verification failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
