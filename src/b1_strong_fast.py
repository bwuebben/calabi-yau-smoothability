#!/usr/bin/env python3
"""
B1-strong test on the FAST integer engine (facets_np + polygon_of_face_int, the
same primitives ks_sweep validates row-by-row against its pure-python reference).

For every both-sides-unit reflexive 4-polytope (Delta unit-edge AND max l(F*)==1),
flag any non-triangle 2-face that is NOT centrally symmetric (edge multiset !=
its negation).  0 flags => B1-strong holds (=> B1, mirror non-isolation).

Uses sweep_file's exact backpressure (BoundedSemaphore) so memory stays flat
(the earlier apply_async scanner had none and thrashed on v08).

Run: ./venv/bin/python src/b1_strong_fast.py data/ks/polytopes-4d-09-vertices.parquet
"""
from collections import Counter
from itertools import combinations
from multiprocessing import Pool
import argparse, os, sys, threading, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from toric_census import smoothing_components
from ks_sweep import (facets_np, polygon_of_face_int, vgcd, vsub,
                      iter_vert_chunks, file_meta, selftest_engines)


def centrally_symmetric(evs):
    return Counter(tuple(e) for e in evs) == Counter((-e[0], -e[1]) for e in evs)


def _empty():
    return dict(n=0, both=0, viol=0, viol_nonsmooth=0, both_polys=[], viol_polys=[])


def _work(args):
    import numpy as np
    verts, subsets = args
    agg = _empty()
    B = verts.shape[0]
    C, n = subsets.shape[0], verts.shape[1]
    step = max(16, min(B, 40_000_000 // max(1, C * n * 8)))
    for s in range(0, B, step):
        vb = verts[s:s + step]
        all_facs = facets_np(vb, subsets)
        for bi in range(vb.shape[0]):
            V = [tuple(int(x) for x in row) for row in vb[bi]]
            facs = all_facs[bi]
            # 2-faces = facet pairs meeting in >=3 shared vertices
            tf = {}
            for (a, b) in combinations(range(len(facs)), 2):
                I = facs[a][1] & facs[b][1]
                if len(I) >= 3 and I not in tf:
                    tf[I] = (a, b)
            all_unit = True
            max_dual = 1
            faces = []
            for I, (a, b) in tf.items():
                u1, u2 = facs[a][0], facs[b][0]
                dl = vgcd(vsub(u1, u2))
                if dl > max_dual:
                    max_dual = dl
                evs, lens = polygon_of_face_int(V, sorted(I), u1, u2)
                if any(l >= 2 for l in lens):
                    all_unit = False
                faces.append((len(evs), tuple(evs)))
            agg["n"] += 1
            if not (all_unit and max_dual == 1):
                continue
            agg["both"] += 1
            if len(agg["both_polys"]) < 60:
                agg["both_polys"].append(V)
            for k, evs in faces:
                if k >= 4 and not centrally_symmetric(evs):
                    agg["viol"] += 1
                    if smoothing_components(list(evs)) == 0:
                        agg["viol_nonsmooth"] += 1
                    if len(agg["viol_polys"]) < 30:
                        agg["viol_polys"].append((k, list(evs), V))
    return agg


def _merge(a, b):
    for k in ("n", "both", "viol", "viol_nonsmooth"):
        a[k] += b[k]
    a["both_polys"].extend(b["both_polys"][:max(0, 60 - len(a["both_polys"]))])
    a["viol_polys"].extend(b["viol_polys"][:max(0, 30 - len(a["viol_polys"]))])
    return a


def sweep(path, procs, chunk_rows=None):
    import numpy as np
    nrows, n = file_meta(path)
    subsets = np.array(list(combinations(range(n), 4)), dtype=np.int64)
    C = subsets.shape[0]
    chunk = chunk_rows or max(200, min(20000, 60_000_000 // max(1, C * n),
                                       (nrows + 4 * procs - 1) // (4 * procs)))
    sem = threading.BoundedSemaphore(procs * 3)

    def jobs():
        for verts in iter_vert_chunks(path, chunk):
            sem.acquire()
            yield (verts, subsets)

    agg = _empty()
    with Pool(processes=procs) as pool:
        for part in pool.imap_unordered(_work, jobs()):
            sem.release()
            agg = _merge(agg, part)
    return agg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--procs", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument("--selftest", type=int, default=60)
    args = ap.parse_args()
    grand = _empty()
    for path in args.files:
        if args.selftest:
            nv = selftest_engines(path, args.selftest)
            print(f"{os.path.basename(path)}: engine selftest (ref==fast) OK on {nv} rows", flush=True)
        t0 = time.time()
        agg = sweep(path, args.procs)
        print(f"{os.path.basename(path)}: n={agg['n']}  both-sides-unit={agg['both']}  "
              f"non-tri-non-symmetric faces={agg['viol']} "
              f"(non-smoothable={agg['viol_nonsmooth']})  [{time.time()-t0:.0f}s]", flush=True)
        grand = _merge(grand, agg)
    print("\n" + "=" * 60)
    print(f"TOTAL: n={grand['n']}  both-sides-unit={grand['both']}  "
          f"non-triangle NON-symmetric both-sides faces={grand['viol']} "
          f"(non-smoothable={grand['viol_nonsmooth']})")
    if grand["viol"] == 0:
        print("==> B1-STRONG HOLDS: every both-sides 2-face is a triangle or a\n"
              "    centrally-symmetric zonotope.  (Implies B1: mirror non-isolation.)")
    else:
        print("==> B1-STRONG violated; examples:")
        for k, evs, V in grand["viol_polys"][:10]:
            sm = "smoothable" if smoothing_components(list(evs)) else "NON-SMOOTHABLE"
            print(f"    k={k} {sm} edges={evs}  V={V}")


if __name__ == "__main__":
    main()
