#!/usr/bin/env python3
"""
Fast BOTH-SIDES-ANY scan (paper 3, Sec. 5): find EVERY both-sides unit
reflexive 4-polytope -- Delta unit-edged AND max over 2-faces of l(F*) == 1
(<=> Delta* unit-edged, by the transverse-multiplicity identity) -- and test
the triangle-or-zonotope conjecture on each hit.

Unlike src/both_sides_ks.py (pure-python reference path, v5-v9), this reuses
the validated FAST engine of src/ks_sweep.py (numpy facet batches + cached
integer polygon classification), so files up to v12-v13 are feasible.  Hits
are rare (5 in v<=9), so each hit is re-verified with the exact reference
toolkit (batyrev_global.analyze) and its faces are tested for
triangle-or-centrally-symmetric directly.

Run:  ./venv/bin/python src/both_sides_fast.py data/ks/polytopes-4d-10-vertices.parquet \
          --procs 8 --json output/both_sides_v10.json
"""
from itertools import combinations
from multiprocessing import Pool
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ks_sweep import (CensusMatcher, classify_polytope_fast, facets_np,     # noqa: E402
                      file_meta, iter_vert_chunks)
from batyrev_global import analyze                                          # noqa: E402

_MATCHER = None


def _init_worker():
    global _MATCHER
    _MATCHER = CensusMatcher()


def _work_chunk(args):
    """Count both-sides-any polytopes in a chunk; return hits' vertex lists."""
    global _MATCHER
    verts, subsets = args
    n_seen, hits = 0, []
    B = verts.shape[0]
    C, n = subsets.shape[0], verts.shape[1]
    step = max(16, min(B, 40_000_000 // max(1, C * n * 8)))
    for s in range(0, B, step):
        vb = verts[s:s + step]
        all_facs = facets_np(vb, subsets)
        for bi in range(vb.shape[0]):
            V = [tuple(int(x) for x in row) for row in vb[bi]]
            r = classify_polytope_fast(V, all_facs[bi], _MATCHER)
            n_seen += 1
            if r["all_unit"] and r["max_dual"] == 1:
                hits.append(V)
    return n_seen, hits


def is_zonotope(evs):
    """Unit-edge polygon centrally symmetric <=> edge vectors in +/- pairs."""
    from collections import Counter
    c = Counter(tuple(e) for e in evs)
    return all(c[(-a, -b)] == m for (a, b), m in c.items())


def verify_hit(V):
    """Exact re-check of one both-sides hit + conjecture test per face.
    Returns dict with face kinds and any asymmetric (counterexample) faces."""
    rep = analyze("hit", [tuple(v) for v in V], verbose=False)
    assert rep is not None, "hit not reflexive?!"
    assert all(l == 1 for l in rep["edges"].values()), "hit not unit-edged?!"
    assert all(f["npoints"] == 1 for f in rep["faces"]), "hit has l(F*)>1?!"
    kinds, bad = [], []
    for f in rep["faces"]:
        if f["k"] == 3:
            kinds.append("triangle:" + f["status"])
        elif is_zonotope(f["edges2d"]):
            kinds.append(f"zonotope k={f['k']}:" + f["status"])
        else:
            kinds.append(f"ASYMMETRIC k={f['k']}:" + f["status"])
            bad.append(f["verts"])
    return dict(verts=V, facekinds=sorted(set(kinds)), asymmetric=bad)


def iter_range(path, chunk, start, stop):
    """iter_vert_chunks restricted to rows [start, stop): skip whole chunks
    before `start`, trim the boundary chunks."""
    done = 0
    for verts in iter_vert_chunks(path, chunk):
        lo, hi = done, done + verts.shape[0]
        done = hi
        if hi <= start:
            continue
        if lo >= stop:
            return
        yield verts[max(0, start - lo): min(verts.shape[0], stop - lo)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--procs", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument("--json", default=None)
    ap.add_argument("--start", type=int, default=0, help="first row (segmented runs)")
    ap.add_argument("--stop", type=int, default=None, help="one past last row")
    args = ap.parse_args()
    import numpy as np
    import threading
    out = {}
    for path in args.files:
        t0 = time.time()
        nrows, n = file_meta(path)
        stop = nrows if args.stop is None else min(args.stop, nrows)
        subsets = np.array(list(combinations(range(n), 4)), dtype=np.int64)
        C = subsets.shape[0]
        chunk = max(200, min(20000, 60_000_000 // max(1, C * n),
                             (nrows + 4 * args.procs - 1) // (4 * args.procs)))
        sem = threading.BoundedSemaphore(args.procs * 3)

        def jobs():
            for verts in iter_range(path, chunk, args.start, stop):
                sem.acquire()
                yield (verts, subsets)

        total, hits = 0, []
        with Pool(processes=args.procs, initializer=_init_worker) as pool:
            for n_seen, hs in pool.imap_unordered(_work_chunk, jobs()):
                sem.release()
                total += n_seen
                hits.extend(hs)
        verified = [verify_hit(V) for V in hits]
        n_asym = sum(1 for v in verified if v["asymmetric"])
        key = os.path.basename(path)
        if args.start or stop < nrows:
            key += f"[{args.start}:{stop}]"
        out[key] = dict(
            n=total, both_sides_any=len(hits), asymmetric_hits=n_asym,
            hits=verified, secs=round(time.time() - t0, 1))
        print(f"{key}: n={total}  both-sides-any={len(hits)}"
              f"  ASYMMETRIC(=counterexample)={n_asym}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        for v in verified:
            print("   hit:", v["facekinds"],
                  ("COUNTEREXAMPLE!" if v["asymmetric"] else "(conjecture ok)"),
                  flush=True)
        if args.json:
            with open(args.json, "w") as f:
                json.dump(out, f, indent=1)
    print("BOTH-SIDES-FAST-DONE", flush=True)


if __name__ == "__main__":
    main()
