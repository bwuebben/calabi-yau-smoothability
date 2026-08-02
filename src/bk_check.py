#!/usr/bin/env python3
"""
Batyrev-Kreuzer cross-check (paper 2): find the BK population — reflexive
4-polytopes ALL of whose 2-faces are standard triangles or unit
parallelograms (=> X has only ODPs) — and tabulate the dual-edge lengths
l(F°) of the square faces.  BK (arXiv:0802.3376) report 198,849 such
polytopes, of which 30,241 pass Namikawa's smoothability criterion.

Hypothesis under test (notes_globaldefect.md): per-face Friedman relations
among the l equal exceptional classes are satisfiable iff l >= 2, matching
the ∂-datum threshold; so #{population: every square face has l >= 2}
should approximate (lower-bound gap = cross-face rescues) BK's 30,241.

Usage: ./venv/bin/python src/bk_check.py data/ks/polytopes-4d-0[5-9]*.parquet --procs 4
"""
from itertools import combinations
from multiprocessing import Pool
import argparse, json, os, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ks_sweep import (facets_np, polygon_of_face_int, file_meta,
                      iter_vert_chunks)
from batyrev_global import vgcd, vsub


def classify_bk(V, facs):
    """None if some 2-face is not a standard triangle / unit parallelogram;
    else list of dual-edge lengths of the square faces."""
    tf = {}
    for (a, b) in combinations(range(len(facs)), 2):
        I = facs[a][1] & facs[b][1]
        if len(I) >= 3 and I not in tf:
            tf[I] = (a, b)
    ells = []
    for I, (a, b) in tf.items():
        u1, u2 = facs[a][0], facs[b][0]
        evs, lens = polygon_of_face_int(V, sorted(I), u1, u2)
        if any(l != 1 for l in lens):
            return None
        k = len(evs)
        if k == 3:
            e0, e1 = evs[0], evs[1]
            if abs(e0[0]*e1[1] - e0[1]*e1[0]) == 1:
                continue                     # standard triangle
            return None
        if k == 4:
            s = sorted(evs)
            if (s[0][0] == -s[3][0] and s[0][1] == -s[3][1]
                    and s[1][0] == -s[2][0] and s[1][1] == -s[2][1]
                    and abs(s[0][0]*s[1][1] - s[0][1]*s[1][0]) == 1):
                ells.append(vgcd(vsub(u1, u2)))   # unit parallelogram
                continue
            return None
        return None
    return ells


def _work(args):
    import numpy as np
    verts, subsets = args
    agg = dict(n=0, pop=0, all_ge2=0, sq_hist={}, ell1_polys=0)
    C, n = subsets.shape[0], verts.shape[1]
    step = max(16, min(verts.shape[0], 40_000_000 // max(1, C * n * 8)))
    for s in range(0, verts.shape[0], step):
        vb = verts[s:s + step]
        all_facs = facets_np(vb, subsets)
        for bi in range(vb.shape[0]):
            agg["n"] += 1
            V = [tuple(int(x) for x in row) for row in vb[bi]]
            ells = classify_bk(V, all_facs[bi])
            if ells is None:
                continue
            k = len(ells)
            agg["sq_hist"][k] = agg["sq_hist"].get(k, 0) + 1
            if k == 0:
                continue                     # all-triangle: X already smooth
            agg["pop"] += 1
            if all(e >= 2 for e in ells):
                agg["all_ge2"] += 1
            else:
                agg["ell1_polys"] += 1
    return agg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--json", default=None)
    args = ap.parse_args()
    import numpy as np
    tot = dict(n=0, pop=0, all_ge2=0, sq_hist={}, ell1_polys=0)
    for path in args.files:
        t0 = time.time()
        nrows, n = file_meta(path)
        subsets = np.array(list(combinations(range(n), 4)), dtype=np.int64)
        chunk = max(200, min(20000, 60_000_000 // max(1, subsets.shape[0] * n),
                             (nrows + 4 * args.procs - 1) // (4 * args.procs)))
        import threading
        sem = threading.BoundedSemaphore(args.procs * 3)
        def jobs():
            for verts in iter_vert_chunks(path, chunk):
                sem.acquire()
                yield (verts, subsets)
        agg = dict(n=0, pop=0, all_ge2=0, sq_hist={}, ell1_polys=0)
        with Pool(processes=args.procs) as pool:
            for part in pool.imap_unordered(_work, jobs()):
                sem.release()
                agg["n"] += part["n"]; agg["pop"] += part["pop"]
                agg["all_ge2"] += part["all_ge2"]
                agg["ell1_polys"] += part["ell1_polys"]
                for k, v in part["sq_hist"].items():
                    agg["sq_hist"][k] = agg["sq_hist"].get(k, 0) + v
        print(f"{os.path.basename(path)}: n={agg['n']}  BK-population={agg['pop']}"
              f"  all-squares-l>=2={agg['all_ge2']}  has-l=1-square={agg['ell1_polys']}"
              f"  #squares hist={dict(sorted(agg['sq_hist'].items()))}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        for k in tot:
            if k == "sq_hist":
                for kk, v in agg["sq_hist"].items():
                    tot["sq_hist"][kk] = tot["sq_hist"].get(kk, 0) + v
            else:
                tot[k] += agg[k]
    print(f"\nTOTAL: population {tot['pop']} (BK 2008: 198,849 over the full DB)"
          f"  all-squares-l>=2: {tot['all_ge2']} (BK smoothable: 30,241)"
          f"  with an l=1 square: {tot['ell1_polys']}")
    if args.json:
        json.dump(tot, open(args.json, "w"), indent=1)


if __name__ == "__main__":
    main()
