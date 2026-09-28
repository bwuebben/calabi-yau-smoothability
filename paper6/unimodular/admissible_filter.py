#!/usr/bin/env python3
"""The admissible reflexive 4-polytopes of the Kreuzer-Skarke list (Appendix A, Step 2).

Run from the repository root, on one per-vertex-count file of the classification:

    ./venv/bin/python paper6/unimodular/admissible_filter.py data/ks/polytopes-4d-07-vertices.parquet \
        --out admissible_v07.jsonl --compare

Requirements: numpy and pyarrow (as for src/ks_sweep.py), the exact toric modules _lib/ks_sweep.py and
_lib/batyrev_global.py included beside this program, and the PALP program poly-4d.x, which ships with SageMath.  Its path
is taken from --palp, else from the environment variable PALP_POLY4D, else from PATH, else from
`sage -sh -c 'command -v poly-4d.x'`.  The parquet files are those of the pinned dataset revision in
paper6/provenance/ks_polytopes_4d_sha256.tsv (see the bundle README); they are not part of the repository.

Every row of the file is treated as the polytope Delta.  PALP (`poly-4d.x -vem`) gives the vertices, the facet
equations <u,x> = 1 and the vertex-facet incidences; the parse is checked against the pairing matrix on every
50th polytope and on every admissible one.  The two-faces of Delta are the maximal vertex sets (of size >= 3)
shared by two facets; each is put into the lattice of its plane by src/ks_sweep.polygon_of_face_int, which
returns its primitive edge vectors in cyclic order and their lattice lengths.  Delta is admissible when every
two-face has primitive edges and, with k vertices and twice its normalised area 2A (by the shoelace formula),
(k, 2A) is (3, 1), (4, 2), (5, 5) or (6, 6): by Pick's theorem these are the unimodular triangle, the unit
square, the dP7 pentagon and the dP6 hexagon; and every non-triangular two-face, contained in the facets with
normals u1, u2, has a dual edge [u1, u2] of lattice length one (gcd of the entries of u1 - u2 equal to 1).
(Primitivity of the edges of Delta is implied, since every edge lies in a two-face with primitive edges.)

Output: one JSON line per admissible polytope, {"row", "V", "two_faces"}, in the conventions of the
certificate files paper6/unimodular/certificates/heights_vNN.jsonl.gz.  With --compare the admissible rows and
their vertex lists are compared with the certificate file of the same number of vertices.
"""
import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(HERE, "_lib"))
from ks_sweep import polygon_of_face_int          # noqa: E402
from batyrev_global import vsub, vgcd              # noqa: E402

POLY = None


def find_palp(opt):
    for cand in (opt, os.environ.get("PALP_POLY4D"), shutil.which("poly-4d.x")):
        if cand and os.path.exists(cand):
            return cand
    try:
        out = subprocess.run(["sage", "-sh", "-c", "command -v poly-4d.x"], capture_output=True, text=True,
                             timeout=120).stdout.strip().splitlines()
        if out and os.path.exists(out[-1]):
            return out[-1]
    except (OSError, subprocess.SubprocessError):
        pass
    sys.exit("PALP's poly-4d.x not found: pass --palp PATH or set PALP_POLY4D")


def area2(evs, lens):
    """Twice the area of the polygon with the given edge vectors and lengths (shoelace formula)."""
    A2, x, y = 0, 0, 0
    for e, l in zip(evs, lens):
        A2 += x * e[1] * l - y * e[0] * l
        x += e[0] * l
        y += e[1] * l
    return abs(A2)


def classify(V, U, inc):
    """V: vertices; U: facet normals with <u,x> = 1 on the facet; inc[k]: vertex indices on facet k.
    Returns the numbers of two-faces of each kind, or None if Delta is not admissible."""
    kinds = {"triangle": 0, "square": 0, "pentagon": 0, "hexagon": 0}
    seen = set()
    nf = len(U)
    for a in range(nf):
        for b in range(a + 1, nf):
            I = inc[a] & inc[b]
            if len(I) < 3 or I in seen:
                continue
            seen.add(I)
            evs, lens = polygon_of_face_int(V, sorted(I), U[a], U[b])
            if any(l != 1 for l in lens):
                return None
            k, A2 = len(evs), area2(evs, lens)
            if (k, A2) == (3, 1):
                kinds["triangle"] += 1
                continue
            kd = {(4, 2): "square", (5, 5): "pentagon", (6, 6): "hexagon"}.get((k, A2))
            if kd is None or vgcd(vsub(U[a], U[b])) != 1:
                return None
            kinds[kd] += 1
    return kinds


def parse_matrix(lines, pos):
    r, c = (int(x) for x in lines[pos].split()[:2])
    M = np.array([[int(x) for x in lines[pos + 1 + i].split()] for i in range(r)], dtype=np.int64)
    assert M.shape == (r, c)
    return M, pos + 1 + r


def work(args):
    start, polys = args
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write("".join(f"{len(V)} 4\n" + "".join(f"{v[0]} {v[1]} {v[2]} {v[3]}\n" for v in V) for V in polys))
        fin = f.name
    try:
        out = subprocess.run([POLY, "-vem", fin], capture_output=True, text=True, check=True).stdout
    finally:
        os.unlink(fin)
    lines = out.split("\n")
    pos, res, n = 0, [], 0
    while pos < len(lines) and lines[pos].strip():
        Vm, pos = parse_matrix(lines, pos)                # vertices (either orientation)
        E, pos = parse_matrix(lines, pos)                 # facet equations <e,x> + 1 >= 0
        Pm, pos = parse_matrix(lines, pos)                # pairing matrix
        if Vm.shape[0] == 4 and Vm.shape[1] != 4:
            Vn = Vm.T
        elif Vm.shape[1] == 4 and Vm.shape[0] != 4:
            Vn = Vm
        else:
            raise ValueError("ambiguous vertex matrix shape")
        assert Pm.shape == (E.shape[0], Vn.shape[0])
        if n % 50 == 0:
            assert (E @ Vn.T + 1 == Pm).all()
        V = [tuple(int(x) for x in row) for row in Vn]
        U = [tuple(-int(x) for x in e) for e in E]
        inc = [frozenset(np.flatnonzero(z).tolist()) for z in (Pm == 0)]
        kd = classify(V, U, inc)
        if kd is not None:
            assert (E @ Vn.T + 1 == Pm).all()
            res.append({"row": start + n, "V": [list(v) for v in V], "two_faces": kd})
        n += 1
    assert n == len(polys), (n, len(polys))
    return len(polys), res


def init(poly):
    global POLY
    POLY = poly


def compare(found, nv):
    path = os.path.join(HERE, "certificates", f"heights_v{nv:02d}.jsonl.gz")
    with gzip.open(path, "rt") as f:
        cert = [json.loads(l) for l in f]
    a = {r["row"]: (sorted(map(tuple, r["V"])), r["two_faces"]) for r in found}
    b = {r["row"]: (sorted(map(tuple, r["V"])), r["two_faces"]) for r in cert}
    same = a == b
    print(f"compare with {os.path.relpath(path, REPO)}: {len(a)} rows found, {len(b)} certificates; "
          f"{'identical rows, vertex sets and two-face counts' if same else 'DIFFERENT'}", flush=True)
    return same


def main():
    global POLY
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("file", help="data/ks/polytopes-4d-NN-vertices.parquet")
    ap.add_argument("--out", default=None, help="write the admissible polytopes to this JSON-lines file")
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--chunk", type=int, default=20000)
    ap.add_argument("--limit", type=int, default=None, help="only the first N rows (for a quick test)")
    ap.add_argument("--palp", default=None, help="path of PALP's poly-4d.x")
    ap.add_argument("--compare", action="store_true", help="compare with the certificate file")
    A = ap.parse_args()
    POLY = find_palp(A.palp)
    import pyarrow.parquet as pq
    pf = pq.ParquetFile(A.file)

    def jobs():
        off = 0
        for batch in pf.iter_batches(batch_size=A.chunk, columns=["vertices"]):
            if A.limit is not None and off >= A.limit:
                return
            rows = batch.column(0).to_pylist()
            if A.limit is not None:
                rows = rows[:A.limit - off]
            yield (off, rows)
            off += batch.num_rows

    t0 = time.time()
    tot, found = 0, []
    with Pool(A.procs, initializer=init, initargs=(POLY,)) as pool:
        for k, (m, res) in enumerate(pool.imap(work, jobs(), chunksize=1)):
            tot += m
            found += res
            if (k + 1) % 50 == 0:
                print(f"[{time.time() - t0:7.0f}s] rows {tot}, admissible {len(found)}", flush=True)
    print(f"{A.file}: rows {tot}, admissible {len(found)}, with a hexagon "
          f"{sum(r['two_faces']['hexagon'] > 0 for r in found)} [{time.time() - t0:.0f}s]", flush=True)
    if A.out:
        with open(A.out, "w") as f:
            for r in found:
                f.write(json.dumps(r) + "\n")
    if A.compare:
        if A.limit is not None:
            print("--compare ignored with --limit")
        else:
            nv = len(found[0]["V"]) if found else int(os.path.basename(A.file).split("-")[2])
            sys.exit(0 if compare(found, nv) else 1)


if __name__ == "__main__":
    main()
