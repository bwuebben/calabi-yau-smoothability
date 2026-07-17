#!/usr/bin/env python3
"""
Test the LOCAL balancing lemma (Track B1 (a)): is there a reflexive 4-polytope
with a 2-face F that is non-triangle, NON-centrally-symmetric, has l(F*)=1, and
whose entire STAR (every 2-face sharing >=1 vertex with F) also has l=1?

0 such F  =>  the local lemma holds: an asymmetric l=1 face forces a long dual
edge in its own star (a first-neighbour obstruction), which implies B1 / B1-strong
by the transverse-multiplicity Lemma.  This local condition is WEAKER than
both-sides (only F's star must be unit-dual, not the whole polytope), so it is a
strictly harder test.

Fast integer engine (facets_np) + per-file ref selftest.
Run: ./venv/bin/python src/b1_local.py data/ks/polytopes-4d-0{5,6,7}-vertices.parquet
"""
from collections import Counter
from itertools import combinations
from multiprocessing import Pool
import argparse, os, sys, threading, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ks_sweep import (facets_np, polygon_of_face_int, vgcd, vsub,
                      iter_vert_chunks, file_meta, selftest_engines)


def centrally_symmetric(evs):
    return Counter(tuple(e) for e in evs) == Counter((-e[0], -e[1]) for e in evs)


def _empty():
    return dict(n=0, cand=0, viol=0, viol_polys=[])


def _work(args):
    import numpy as np
    verts, subsets = args
    agg = _empty()
    C, n = subsets.shape[0], verts.shape[1]
    step = max(16, min(verts.shape[0], 40_000_000 // max(1, C * n * 8)))
    for s in range(0, verts.shape[0], step):
        vb = verts[s:s + step]
        all_facs = facets_np(vb, subsets)
        for bi in range(vb.shape[0]):
            V = [tuple(int(x) for x in row) for row in vb[bi]]
            facs = all_facs[bi]
            faces = []                       # (idxset, l, k, evs, unit)
            for (a, b) in combinations(range(len(facs)), 2):
                I = facs[a][1] & facs[b][1]
                if len(I) >= 3:
                    u1, u2 = facs[a][0], facs[b][0]
                    l = vgcd(vsub(u1, u2))
                    evs, lens = polygon_of_face_int(V, sorted(I), u1, u2)
                    unit = not any(x >= 2 for x in lens)
                    faces.append((I, l, len(evs), tuple(evs), unit))
            agg["n"] += 1
            for (I, l, k, evs, unit) in faces:
                # TRULY-LOCAL condition: F unit-edge; ignore far edges of Delta
                if k >= 4 and l == 1 and unit and not centrally_symmetric(evs):
                    agg["cand"] += 1
                    # star: every OTHER 2-face sharing >=1 vertex must be unit-edge AND l==1
                    star_ok = True
                    for (J, lJ, kJ, eJ, uJ) in faces:
                        if J is I:
                            continue
                        if (I & J) and (lJ >= 2 or not uJ):
                            star_ok = False
                            break
                    if star_ok:
                        agg["viol"] += 1
                        if len(agg["viol_polys"]) < 20:
                            agg["viol_polys"].append((k, list(evs), V))
    return agg


def _merge(a, b):
    for k in ("n", "cand", "viol"):
        a[k] += b[k]
    a["viol_polys"].extend(b["viol_polys"][:max(0, 20 - len(a["viol_polys"]))])
    return a


def sweep(path, procs):
    import numpy as np
    nrows, n = file_meta(path)
    subsets = np.array(list(combinations(range(n), 4)), dtype=np.int64)
    C = subsets.shape[0]
    chunk = max(200, min(20000, 60_000_000 // max(1, C * n),
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
    args = ap.parse_args()
    grand = _empty()
    for path in args.files:
        selftest_engines(path, 40)
        t0 = time.time()
        agg = sweep(path, args.procs)
        print(f"{os.path.basename(path)}: n={agg['n']}  asymmetric-l1-nontriangle faces={agg['cand']}  "
              f"LOCAL-LEMMA VIOLATIONS (star all unit)={agg['viol']}  [{time.time()-t0:.0f}s]", flush=True)
        grand = _merge(grand, agg)
    print("\n" + "=" * 60)
    print(f"TOTAL: n={grand['n']}  candidate asymmetric l=1 non-triangle faces={grand['cand']}  "
          f"local-lemma violations={grand['viol']}")
    if grand["viol"] == 0:
        print("==> LOCAL BALANCING LEMMA HOLDS: an asymmetric l=1 face always has a\n"
              "    long dual edge in its own star.  (First-neighbour obstruction => B1.)")
    else:
        for k, evs, V in grand["viol_polys"][:10]:
            print(f"    VIOLATION k={k} edges={evs}  V={V}")


if __name__ == "__main__":
    main()
