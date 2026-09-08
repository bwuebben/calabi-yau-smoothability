#!/usr/bin/env python3
"""Check the pinned KS parquet copy; see ../DATASET_CHECK.md for scope.

inventory: SHA-256, byte sizes, coverage, parquet row counts.
check --mode stored: exact coordinate-set duplicates (not GL(4,Z) duplicates).
check --mode lattice: PALP reflexivity and canonical GL(4,Z) normal forms.
benchmark: deterministic samples, including SQLite insertion, in both modes.

Full keys, not hashes of polytopes, are stored in SQLite. The lattice check
discharges the database-copy hypothesis only after all 473,800,776 rows pass.
It relies on PALP and the published KS cardinality, and does not check the
singularity predicates or the papers' census results.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import time

TOTAL = 473_800_776
REVISION = "60c0e119a03608418df538191f65da3f43b5b819"
ROOT = Path(__file__).resolve().parent.parent
NAME = re.compile(r"polytopes-4d-(\d{2})-vertices\.parquet\Z")


class CheckError(RuntimeError):
    pass


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def manifest_entries(path):
    entries = {}
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        digest, size, name = line.split()
        if not NAME.fullmatch(name) or name in entries:
            raise CheckError(f"invalid or repeated manifest filename: {name}")
        if not re.fullmatch(r"[0-9a-f]{64}", digest) or int(size) <= 0:
            raise CheckError(f"invalid manifest entry: {name}")
        entries[name] = {"sha256": digest, "bytes": int(size)}
    expected = {f"polytopes-4d-{n:02d}-vertices.parquet" for n in [*range(5, 34), 36]}
    if entries.keys() != expected:
        raise CheckError("manifest must list each 5--33 and 36 vertex file exactly once")
    return entries


def inventory(data_dir, manifest):
    import pyarrow.parquet as pq

    entries = manifest_entries(manifest)
    present, missing = [], []
    unexpected = sorted(p.name for p in data_dir.glob("*.parquet") if p.name not in entries)
    if unexpected:
        raise CheckError(f"unexpected parquet files in input directory: {unexpected}")
    for name, expected in sorted(entries.items()):
        path = data_dir / name
        if not path.is_file():
            missing.append(name)
            continue
        if path.stat().st_size != expected["bytes"] or sha256(path) != expected["sha256"]:
            raise CheckError(f"size or SHA-256 mismatch: {path}")
        f = pq.ParquetFile(path)
        if not {"vertices", "vertex_count", "facet_count"} <= set(f.schema_arrow.names):
            raise CheckError(f"missing required columns: {path}")
        if f.metadata.num_rows <= 0:
            raise CheckError(f"empty parquet: {path}")
        present.append({"name": name, **expected, "rows": f.metadata.num_rows,
                        "vertices": int(NAME.fullmatch(name)[1])})
    rows = sum(p["rows"] for p in present)
    if rows > TOTAL or (not missing and rows != TOTAL):
        raise CheckError(f"row count {rows:,} conflicts with KS total {TOTAL:,}")
    return {"dataset_revision": REVISION, "manifest_sha256": sha256(manifest),
            "files": present, "missing": missing, "present_rows": rows,
            "expected_rows": TOTAL, "present_bytes": sum(p["bytes"] for p in present),
            "expected_bytes": sum(p["bytes"] for p in entries.values()),
            "coverage_complete": not missing and rows == TOTAL,
            "classification_verified": False}


def validate_row(row, vertex_count):
    vertices = row["vertices"]
    if not isinstance(vertices, list) or len(vertices) != vertex_count:
        raise CheckError("vertex list length disagrees with filename")
    if type(row["vertex_count"]) is not int or row["vertex_count"] != vertex_count:
        raise CheckError("vertex_count column disagrees with filename")
    if type(row["facet_count"]) is not int or row["facet_count"] < 5:
        raise CheckError("invalid facet_count")
    if any(not isinstance(v, list) or len(v) != 4 or
           any(type(x) is not int for x in v) for v in vertices):
        raise CheckError("vertices must be integer vectors of dimension four")
    if len({tuple(v) for v in vertices}) != vertex_count:
        raise CheckError("repeated vertex inside a row")
    return vertices


def encode_key(vertices):
    """Injective variable-length encoding of an ordered integer matrix."""
    out = bytearray()
    for x in [len(vertices), *(x for v in vertices for x in v)]:
        n = 2 * x if x >= 0 else -2 * x - 1
        while n >= 128:
            out.append((n & 127) | 128)
            n >>= 7
        out.append(n)
    return bytes(out)


def find_palp(explicit=None):
    if explicit:
        result = shutil.which(explicit)
        if not result:
            raise CheckError(f"PALP executable not found: {explicit}")
        return result
    for name in ("poly.x", "poly-6d.x", "poly-4d.x"):
        if result := shutil.which(name):
            return result
    if shutil.which("sage"):
        p = subprocess.run(["sage", "-python", "-c",
                            "import shutil; print(shutil.which('poly.x') or '')"],
                           capture_output=True, text=True, timeout=60)
        candidate = p.stdout.strip()
        if p.returncode == 0 and Path(candidate).is_file():
            return candidate
    raise CheckError("PALP unavailable; install PALP or pass --palp /path/to/poly.x")


def palp_keys(rows, vertex_count, executable, timeout=120):
    """Batch PALP calls; fail closed on non-reflexive or unexpected output."""
    inputs = []
    for row in rows:
        vertices = validate_row(row, vertex_count)
        inputs.append(f"4 {vertex_count}\n" + "\n".join(
            " ".join(str(v[j]) for v in vertices) for j in range(4)))
    p = subprocess.run([executable, "-fneN"], input="\n".join(inputs) + "\n",
                       text=True, capture_output=True, timeout=timeout)
    if p.returncode or p.stderr.strip():
        raise CheckError(f"PALP failed: {p.stderr[-1000:]} {p.stdout[-1000:]}")
    lines = iter(line.strip() for line in p.stdout.splitlines() if line.strip())
    keys = []
    try:
        for row in rows:
            header = next(lines)
            match = re.fullmatch(r"(\d+) 4\s+Vertices of P-dual <-> Equations of P", header)
            if not match:
                raise CheckError(f"PALP did not certify reflexivity: {header}")
            facets = int(match[1])
            if facets != row["facet_count"]:
                raise CheckError("recomputed facet count disagrees with dataset")
            for _ in range(facets):
                if len([int(x) for x in next(lines).split()]) != 4:
                    raise CheckError("malformed PALP facet equation")
            header = next(lines)
            if not re.fullmatch(rf"4 {vertex_count}\s+Normal form of vertices of P.*", header):
                raise CheckError(f"dimension or actual vertex count mismatch: {header}")
            matrix = [[int(x) for x in next(lines).split()] for _ in range(4)]
            if any(len(v) != vertex_count for v in matrix):
                raise CheckError("malformed PALP normal form")
            keys.append(encode_key(list(zip(*matrix))))
        if next(lines, None) is not None:
            raise CheckError("extra PALP output")
    except (StopIteration, ValueError) as exc:
        raise CheckError("truncated or malformed PALP output") from exc
    return keys


def batch_keys(rows, vertex_count, mode, palp, timeout):
    if mode == "lattice":
        return palp_keys(rows, vertex_count, palp, timeout)
    return [encode_key(sorted(validate_row(row, vertex_count))) for row in rows]


def open_database(path, identity):
    db = sqlite3.connect(path)
    db.execute("PRAGMA journal_mode=WAL")
    db.execute("PRAGMA synchronous=FULL")
    db.execute("PRAGMA cache_size=-65536")
    db.execute("CREATE TABLE IF NOT EXISTS keys (key BLOB PRIMARY KEY, row_number INTEGER NOT NULL) WITHOUT ROWID")
    db.execute("CREATE TABLE IF NOT EXISTS state (id INTEGER PRIMARY KEY CHECK(id=1), identity TEXT NOT NULL, next_row INTEGER NOT NULL)")
    saved = db.execute("SELECT identity, next_row FROM state WHERE id=1").fetchone()
    identity_text = json.dumps(identity, sort_keys=True)
    if saved:
        if saved[0] != identity_text:
            db.close()
            raise CheckError("checkpoint identity changed; use a new --work-dir")
        offset = saved[1]
    else:
        if db.execute("SELECT 1 FROM keys LIMIT 1").fetchone():
            db.close()
            raise CheckError("checkpoint has keys but no state")
        db.execute("INSERT INTO state VALUES (1, ?, 0)", (identity_text,))
        db.commit()
        offset = 0
    return db, offset


def insert_keys(db, keys, offset):
    try:
        with db:
            for i, key in enumerate(keys, offset):
                try:
                    db.execute("INSERT INTO keys VALUES (?, ?)", (key, i))
                except sqlite3.IntegrityError as exc:
                    previous = db.execute("SELECT row_number FROM keys WHERE key=?", (key,)).fetchone()[0]
                    raise CheckError(f"duplicate rows {previous} and {i} (zero based)") from exc
            db.execute("UPDATE state SET next_row=? WHERE id=1", (offset + len(keys),))
    except CheckError:
        raise


def iter_rows(path, batch_size, start=0):
    import pyarrow.parquet as pq

    f = pq.ParquetFile(path)
    position = 0
    for group in range(f.num_row_groups):
        count = f.metadata.row_group(group).num_rows
        if position + count <= start:
            position += count
            continue
        for batch in f.iter_batches(batch_size=batch_size, row_groups=[group],
                                    columns=["vertices", "vertex_count", "facet_count"]):
            stop = position + len(batch)
            if stop > start:
                yield batch.slice(max(0, start - position)).to_pylist()
            position = stop


def sample_rows(path, count, seed):
    """Uniform deterministic row sample, materializing only sampled Python rows."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    f = pq.ParquetFile(path)
    indices = sorted(random.Random(seed).sample(range(f.metadata.num_rows), min(count, f.metadata.num_rows)))
    rows, offset, cursor = [], 0, 0
    for batch in f.iter_batches(batch_size=8192, columns=["vertices", "vertex_count", "facet_count"]):
        selected = []
        while cursor < len(indices) and indices[cursor] < offset + len(batch):
            selected.append(indices[cursor] - offset)
            cursor += 1
        if selected:
            rows.extend(batch.take(pa.array(selected, type=pa.int64())).to_pylist())
        offset += len(batch)
        if cursor == len(indices):
            break
    return rows


def check_file(info, args, palp):
    path = args.data_dir / info["name"]
    cache = args.work_dir / f"{path.stem}-{args.mode}.sqlite"
    identity = {"input_sha256": info["sha256"], "rows": info["rows"], "mode": args.mode,
                "checker_sha256": sha256(__file__), "palp_sha256": sha256(palp) if palp else None}
    db, offset = open_database(cache, identity)
    start, initial, last = time.perf_counter(), offset, time.perf_counter()
    try:
        if not 0 <= offset <= info["rows"]:
            raise CheckError("checkpoint row offset outside input")
        for rows in iter_rows(path, args.batch_size, offset):
            keys = batch_keys(rows, info["vertices"], args.mode, palp, args.palp_timeout)
            insert_keys(db, keys, offset)
            offset += len(rows)
            if time.perf_counter() - last >= 10:
                print(f"{path.name}: {offset:,}/{info['rows']:,}", flush=True)
                last = time.perf_counter()
        if offset != info["rows"] or db.execute("SELECT count(*) FROM keys").fetchone()[0] != offset:
            raise CheckError("decoded row count or checkpoint key count mismatch")
        if db.execute("PRAGMA quick_check").fetchone()[0] != "ok":
            raise CheckError("SQLite integrity check failed")
    finally:
        db.close()
    return {"name": info["name"], "mode": args.mode, "rows_checked": offset,
            "new_rows": offset - initial, "seconds": time.perf_counter() - start,
            "database_bytes": cache.stat().st_size, "input_sha256": info["sha256"]}


def benchmark(info, args, palp):
    sample_start = time.perf_counter()
    rows = sample_rows(args.data_dir / info["name"], args.sample_per_file, args.seed + info["vertices"])
    sample_seconds = time.perf_counter() - sample_start
    result = {"name": info["name"], "samples": len(rows), "sample_read_seconds": sample_seconds}
    for mode in ("stored", "lattice"):
        with tempfile.TemporaryDirectory(prefix="ks-benchmark-", dir=args.work_dir) as directory:
            start = time.perf_counter()
            path = Path(directory) / "keys.sqlite"
            db, _ = open_database(path, {"benchmark": mode})
            try:
                for offset in range(0, len(rows), args.batch_size):
                    batch = rows[offset:offset + args.batch_size]
                    keys = batch_keys(batch, info["vertices"], mode, palp, args.palp_timeout)
                    insert_keys(db, keys, offset)
            finally:
                db.close()
            seconds = time.perf_counter() - start
            result[mode] = {"seconds": seconds, "rows_per_second": len(rows) / seconds,
                            "linear_full_database_days": seconds / len(rows) * TOTAL / 86400,
                            "sqlite_bytes_per_sample_row": path.stat().st_size / len(rows)}
    return result


def write_report(path, report):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=["inventory", "check", "benchmark"])
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data/ks")
    parser.add_argument("--manifest", type=Path, default=ROOT / "manifests/ks_polytopes_4d_sha256.tsv")
    parser.add_argument("--work-dir", type=Path, default=ROOT / "data/ks-check")
    parser.add_argument("--report", type=Path)
    parser.add_argument("--mode", choices=["stored", "lattice"], default="lattice")
    parser.add_argument("--allow-partial", action="store_true", help="check available files; never claims full coverage")
    parser.add_argument("--vertices", nargs="+", type=int, help="select available vertex-count files")
    parser.add_argument("--palp", help="path to poly.x; otherwise find on PATH or through Sage")
    parser.add_argument("--palp-timeout", type=float, default=120, help="seconds per batch; failure is checkpointed")
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--sample-per-file", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260908)
    args = parser.parse_args(argv)
    if args.batch_size <= 0 or args.sample_per_file <= 0 or args.palp_timeout <= 0:
        parser.error("batch size, sample count, and timeout must be positive")
    if not args.data_dir.is_dir():
        parser.error("data directory does not exist")
    report = {"command": args.command, "checker_sha256": sha256(__file__), "status": "incomplete",
              "classification_verified": False}
    report_path = args.report or args.work_dir / f"{args.command}-{args.mode}.json"
    start = time.perf_counter()
    try:
        inv = inventory(args.data_dir, args.manifest)
        report["inventory"] = inv
        print(f"Verified hashes: {len(inv['files'])}/30 files, {inv['present_rows']:,}/{TOTAL:,} rows; "
              f"{len(inv['missing'])} files missing.", flush=True)
        files = inv["files"]
        if args.vertices:
            if len(set(args.vertices)) != len(args.vertices):
                raise CheckError("repeated --vertices argument")
            files = [f for f in files if f["vertices"] in args.vertices]
            if {f["vertices"] for f in files} != set(args.vertices):
                raise CheckError("a selected vertex-count file is absent")
        report["results"] = []
        if args.command == "inventory":
            report["status"] = "integrity_passed"
        else:
            if not files:
                raise CheckError("no available inputs")
            if args.command == "check" and not args.allow_partial and (inv["missing"] or len(files) != 30):
                raise CheckError("full check requires all 30 files; use --allow-partial for an explicitly partial check")
            args.work_dir.mkdir(parents=True, exist_ok=True)
            palp = find_palp(args.palp) if args.mode == "lattice" or args.command == "benchmark" else None
            if palp:
                report["palp"] = {"path": str(Path(palp).resolve()), "sha256": sha256(palp)}
            for info in files:
                result = benchmark(info, args, palp) if args.command == "benchmark" else check_file(info, args, palp)
                report["results"].append(result)
                print(json.dumps(result, sort_keys=True), flush=True)
                write_report(report_path, report)
            report["classification_verified"] = (args.command == "check" and args.mode == "lattice"
                                                   and inv["coverage_complete"] and len(files) == 30
                                                   and sum(r["rows_checked"] for r in report["results"]) == TOTAL)
            report["status"] = "sample_benchmark_passed" if args.command == "benchmark" else "selected_files_passed"
            if report["classification_verified"]:
                report["status"] = "full_lattice_check_passed"
        code = 0
    except (CheckError, OSError, ValueError, sqlite3.DatabaseError, subprocess.TimeoutExpired) as exc:
        report["status"], report["error"] = "failed", str(exc)
        print(f"FAILED: {exc}", file=sys.stderr, flush=True)
        code = 1
    except KeyboardInterrupt:
        report["status"] = "interrupted"
        print("Interrupted; completed check batches can be resumed.", file=sys.stderr)
        code = 130
    report["seconds"] = time.perf_counter() - start
    write_report(report_path, report)
    print(f"Report: {report_path}; classification_verified={report['classification_verified']}", flush=True)
    return code


if __name__ == "__main__":
    sys.exit(main())
