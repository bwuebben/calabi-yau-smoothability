#!/usr/bin/env python3
"""
Option 2: scan the LOW-vertex Kreuzer-Skarke list for BOTH-SIDES-UNIT reflexive
4-polytopes with a non-smoothable 2-face -- over ALL reflexive polytopes (not
just the planted family of Option 1).

Condition (see both_sides_search.py):
  Delta unit-edged  AND  Delta* unit-edged  AND  a non-smoothable 2-face.
Since edges of Delta* <-> 2-faces F of Delta with length l(F*), "Delta* unit"
<=> max over 2-faces of l(F*) == 1.  We compute that max directly.

Uses the pure-python classification of ks_sweep (cross-validated there against
the dataset Hodge numbers and the fast engine).  Reads parquet with pyarrow;
streams in chunks; multiprocessed.

Run (venv with numpy+pyarrow):
  ./venv/bin/python src/both_sides_ks.py data/ks/polytopes-4d-05-vertices.parquet ...
"""
from collections import Counter
from itertools import combinations
from multiprocessing import Pool
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from toric_census import rigid, smoothing_components                        # noqa: E402
from ks_sweep import facets_reflexive, two_faces_of, polygon_of_face        # noqa: E402
from batyrev_global import vgcd, vsub                                       # noqa: E402


def classify_bs(V):
    """Both-sides classification of one reflexive polytope V (list of 4-tuples).
    Returns all_unit (Delta edges), max_dual = max l(F*) over 2-faces
    (= max edge length of Delta*), and the non-smoothable unit-edge faces."""
    facs = facets_reflexive(V)
    tf = two_faces_of(V, facs)
    all_unit = True
    max_dual = 1
    nbad_unit = 0
    bad_faces = []                                    # (status, k, i) per bad face
    for I, (a, b) in tf.items():
        u1, u2 = facs[a][0], facs[b][0]
        evs, lens = polygon_of_face(V, I, u1, u2)
        dl = vgcd(vsub(u1, u2))                       # l(F*) of this 2-face
        if dl > max_dual:
            max_dual = dl
        if any(l >= 2 for l in lens):
            all_unit = False
            continue
        k = len(evs)
        if k == 3:                                    # standard triangle => smooth
            aa, bb = evs[0], evs[1]
            if abs(aa[0] * bb[1] - aa[1] * bb[0]) == 1:
                continue
        if smoothing_components(list(evs)) == 0:      # rigid or def-only
            nbad_unit += 1
            A2, x, y = 0, 0, 0
            for e in evs:
                A2 += x * e[1] - y * e[0]; x += e[0]; y += e[1]
            A2 = abs(A2)
            ii = (A2 - k + 2) // 2
            bad_faces.append(("rigid" if rigid(list(evs)) else "def-only", k, ii))
    both = all_unit and (max_dual == 1)
    # a NON-TRIANGLE bad face on a both-sides polytope is the prize (non-quotient)
    nontri = both and any(k >= 4 for _, k, _ in bad_faces)
    return dict(all_unit=all_unit, max_dual=max_dual, nbad_unit=nbad_unit,
                bad_faces=bad_faces, both_sides=both,
                hit=(both and nbad_unit > 0), nontri_hit=nontri)


def selftest():
    """max l(F*) must match the mirror-run dual-edge maxima:
    Delta_A -> 9, Delta_B -> 9, Delta_N -> 7; and none is both-sides."""
    from toric_census import ccw_sort, verts_from_edges
    F = [(x, y, 1, 0) for (x, y) in [(0, 0), (-2, -1), (-2, -2), (-1, -2)]]
    A = F + [(0, 0, 0, 1), (1, 1, -1, -1)]
    P = [(x, y, 1, 0) for (x, y) in [(0, 0), (-2, -1), (-3, -2), (-2, -3), (-1, -2)]]
    B = P + [(0, 0, 0, 1), (1, 1, -1, -1)]
    sq = verts_from_edges(ccw_sort([(1, 0), (0, 1), (-1, 0), (0, -1)]))
    N = [(x, y, 1, 0) for (x, y) in sq] + [(0, 0, 0, 1), (0, 1, -1, 0), (-1, -2, 1, -1)]
    for nm, V, exp_max, exp_bad in [("A", A, 9, 1), ("B", B, 9, 1), ("N", N, 7, 0)]:
        r = classify_bs([tuple(v) for v in V])
        assert r["max_dual"] == exp_max, (nm, r["max_dual"], exp_max)
        assert r["nbad_unit"] == exp_bad, (nm, r["nbad_unit"], exp_bad)
        assert not r["both_sides"], nm
    print("selftest: max l(F*) = 9,9,7 on Delta_A,B,N; none both-sides — OK",
          flush=True)


def _empty():
    return dict(n=0, nonsmoothable=0, both_sides_any=0, hits=0, nontri_hits=0,
                maxdual_hist_ns=Counter(), hit_facekinds=Counter(),
                hit_polys=[], nontri_polys=[])


def _work(rows):
    agg = _empty()
    for V in rows:
        V = [tuple(int(x) for x in v) for v in V]
        r = classify_bs(V)
        agg["n"] += 1
        if r["both_sides"]:
            agg["both_sides_any"] += 1
        if r["nbad_unit"] > 0:
            agg["nonsmoothable"] += 1
            agg["maxdual_hist_ns"][r["max_dual"]] += 1
            if r["hit"]:
                agg["hits"] += 1
                for (st, k, ii) in r["bad_faces"]:
                    agg["hit_facekinds"][(st, k)] += 1
                if len(agg["hit_polys"]) < 30:
                    agg["hit_polys"].append((V, r["bad_faces"]))
                if r["nontri_hit"]:
                    agg["nontri_hits"] += 1
                    if len(agg["nontri_polys"]) < 30:
                        agg["nontri_polys"].append((V, r["bad_faces"]))
    return agg


def _merge(a, b):
    a["n"] += b["n"]; a["nonsmoothable"] += b["nonsmoothable"]
    a["both_sides_any"] += b["both_sides_any"]; a["hits"] += b["hits"]
    a["nontri_hits"] += b["nontri_hits"]
    a["maxdual_hist_ns"] += b["maxdual_hist_ns"]
    a["hit_facekinds"] += b["hit_facekinds"]
    a["hit_polys"].extend(b["hit_polys"][:max(0, 30 - len(a["hit_polys"]))])
    a["nontri_polys"].extend(b["nontri_polys"][:max(0, 30 - len(a["nontri_polys"]))])
    return a


def sweep(path, procs, chunk=3000, limit=None):
    import pyarrow.parquet as pq
    pf = pq.ParquetFile(path)
    agg = _empty()
    buf, done = [], 0
    with Pool(procs) as pool:
        pending = []
        for batch in pf.iter_batches(batch_size=chunk, columns=["vertices"]):
            rows = batch.column(0).to_pylist()
            if limit is not None:
                rows = rows[:max(0, limit - done)]
            done += len(rows)
            pending.append(pool.apply_async(_work, (rows,)))
            if len(pending) >= procs * 3:
                agg = _merge(agg, pending.pop(0).get())
            if limit is not None and done >= limit:
                break
        for p in pending:
            agg = _merge(agg, p.get())
    return agg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--procs", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    selftest()
    grand = _empty()
    for path in args.files:
        t0 = time.time()
        agg = sweep(path, args.procs, limit=args.limit)
        name = os.path.basename(path)
        hist = dict(sorted(agg["maxdual_hist_ns"].items()))
        min_md = min(agg["maxdual_hist_ns"]) if agg["maxdual_hist_ns"] else None
        print(f"{name}: n={agg['n']}  non-smoothable(unit face)={agg['nonsmoothable']}  "
              f"both-sides-unit(any)={agg['both_sides_any']}  "
              f"BOTH-SIDES+NON-SMOOTHABLE={agg['hits']}  "
              f"** NON-TRIANGLE (non-quotient) hits = {agg['nontri_hits']} **  "
              f"[{time.time()-t0:.0f}s]", flush=True)
        print(f"    min max l(F*) among non-smoothable = {min_md}; "
              f"hit face kinds (status,k)->count = {dict(agg['hit_facekinds'])}",
              flush=True)
        for V, bf in agg["nontri_polys"][:10]:
            print(f"    *** NON-TRIANGLE HIT: bad_faces={bf}  V={V}", flush=True)
        grand = _merge(grand, agg)

    print("\n" + "=" * 72)
    hist = dict(sorted(grand["maxdual_hist_ns"].items()))
    min_md = min(grand["maxdual_hist_ns"]) if grand["maxdual_hist_ns"] else None
    print(f"TOTAL: n={grand['n']}  non-smoothable={grand['nonsmoothable']}  "
          f"both-sides-unit(any)={grand['both_sides_any']}  "
          f"BOTH-SIDES+NON-SMOOTHABLE={grand['hits']}  "
          f"NON-TRIANGLE hits={grand['nontri_hits']}")
    print(f"  hit face kinds (status,k)->count = {dict(grand['hit_facekinds'])}")
    print(f"  among non-smoothable polytopes: min max l(F*) = {min_md}   "
          f"histogram of max l(F*) = {hist}")
    if grand["hits"] == 0:
        print("  => NO both-sides-unit reflexive polytope (<=9 vtx) has a "
              "non-smoothable face.")
    elif grand["nontri_hits"] == 0:
        print("  => both-sides-unit + non-smoothable EXISTS, but EVERY such face "
              "is a triangle\n     (quotient singularity): the novel non-quotient "
              "(F1) and def-only germs\n     are ALWAYS non-isolated on the Batyrev "
              "mirror (<=9 vtx).")
    else:
        print(f"  => {grand['nontri_hits']} both-sides-unit polytopes carry a "
              "NON-TRIANGLE non-smoothable face:\n     a clean two-sided mirror "
              "statement is possible beyond the quotient case.")
    if args.json:
        with open(args.json, "w") as f:
            json.dump({"n": grand["n"], "nonsmoothable": grand["nonsmoothable"],
                       "both_sides_any": grand["both_sides_any"],
                       "hits": grand["hits"], "nontri_hits": grand["nontri_hits"],
                       "hit_facekinds": {f"{s}_k{k}": v
                                         for (s, k), v in grand["hit_facekinds"].items()},
                       "hist": {str(k): v for k, v in hist.items()},
                       "nontri_polys": grand["nontri_polys"]}, f, indent=1)
        print(f"  results -> {args.json}")


if __name__ == "__main__":
    main()
