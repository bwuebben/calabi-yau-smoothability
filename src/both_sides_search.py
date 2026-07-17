#!/usr/bin/env python3
"""
Synthetic search (Option 1) for BOTH-SIDES-UNIT reflexive 4-polytopes with a
non-smoothable 2-face.

Target: Delta reflexive with
  (i)   every edge of Delta unit         (X has only ISOLATED point sings),
  (ii)  every edge of Delta* unit        (X* has only ISOLATED point sings),
  (iii) a non-smoothable 2-face (R/D)    (X non-smoothable).
Such a Delta would let the SAME R/D/S trichotomy run on X and on the Batyrev
mirror X*, so a clean global mirror statement becomes possible.

Key reduction: edges of Delta* <-> 2-faces F of Delta, with lattice length
l(F*).  So (ii) <=> max over 2-faces of l(F*) == 1.  In particular the PLANTED
non-smoothable face must itself have l(F*) = 1, which the two-vertex completion
family can never achieve (it forces l(F*) >= 2).  So we search only the
l(F*)=1 completion families of plant_search:
  * the single-point design  extras = {e4, w2, w3=(*,*,1,-r)}  (lo=0, hi=1);
  * general three-vertex completions {e4, w2, w3} over [-B,B]^4.
Each planted candidate is graded on the FULL both-sides condition, not just on
Delta's own edges.  We report any both-sides hit, and per class the best
(minimal) achievable max l(F*) so we can see how close the family gets.
"""
from itertools import product
from multiprocessing import Pool
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plant_search import nonsmoothable_classes, embed, planted_facets   # noqa: E402
from batyrev_global import (two_faces, face_lattice_polygon,            # noqa: E402
                            classify_polygon, vgcd, vsub)


def grade_both_sides(V, facs, want_idx, npoints):
    """Full both-sides grade of a planted candidate."""
    tf = two_faces(V, facs)
    all_unit = True
    others_smooth = True
    dual_lens = []
    for I, fpair in tf:
        u1, u2 = facs[fpair[0]][0], facs[fpair[1]][0]
        _, evs, lens = face_lattice_polygon(V, I, u1, u2)
        if any(l >= 2 for l in lens):
            all_unit = False
        dual_lens.append(vgcd(vsub(u1, u2)))         # l(F*) of this 2-face
        if I != want_idx and classify_polygon(evs, lens)["status"] != "smooth":
            others_smooth = False
    max_dual = max(dual_lens)
    return dict(npoints=npoints, all_unit=all_unit, others_smooth=others_smooth,
                max_dual=max_dual, dual_multiset=sorted(dual_lens),
                both_sides_unit=(all_unit and max_dual == 1),
                clean=(all_unit and others_smooth and max_dual == 1))


def _consider(best, g, w):
    """Keep the record with (both_sides, clean, then smallest max_dual)."""
    g = dict(g); g["w"] = [list(x) for x in w]
    key = (g["both_sides_unit"], g["clean"], -g["max_dual"])
    if best is None:
        return g
    bkey = (best["both_sides_unit"], best["clean"], -best["max_dual"])
    return g if key > bkey else best


def search_class(edges, B=2, RMAX=2, QMAX=2, B3=None):
    """All l(F*)=1 completions of the planted census class `edges`; return the
    best record (prefers a genuine both-sides hit)."""
    F4 = embed(edges)
    fidx = frozenset(range(len(F4)))
    e4 = (0, 0, 0, 1)
    best = None

    # (1) single-point design: w3 = (c1,c2,1,-r) pins lo=0, e4 pins hi=1.
    for c1 in range(-B, B + 1):
        for c2 in range(-B, B + 1):
            for rr in range(1, RMAX + 1):
                w3 = (c1, c2, 1, -rr)
                for b1 in range(-B, B + 1):
                    for b2 in range(-B, B + 1):
                        for b3 in range(-B, 0):
                            for e in range(-QMAX, QMAX + 1):
                                w2 = (b1, b2, b3, e)
                                if vgcd(w2) != 1:
                                    continue
                                r = planted_facets(F4, [e4, w2, w3])
                                if r is None:
                                    continue
                                facs, npts = r
                                if npts != 1:
                                    continue
                                g = grade_both_sides(F4 + [e4, w2, w3], facs,
                                                     fidx, npts)
                                best = _consider(best, g, [e4, w2, w3])
                                if g["both_sides_unit"]:
                                    return best

    # (2) general three-vertex completions over [-B,B]^4 (keep only l(F*)=1).
    B3 = B if B3 is None else B3
    W = [w for w in product(range(-B3, B3 + 1), repeat=4)
         if w != (0, 0, 0, 0) and vgcd(w) == 1
         and not (w[2] == 1 and w[3] == 0) and w != e4]
    for i, w2 in enumerate(W):
        for w3 in W[i + 1:]:
            if w2[3] >= 0 and w3[3] >= 0:
                continue
            if w2[2] >= 0 and w3[2] >= 0:
                continue
            if w2[2] + w2[3] == 1 and w3[2] + w3[3] == 1:
                continue
            r = planted_facets(F4, [e4, w2, w3])
            if r is None:
                continue
            facs, npts = r
            if npts != 1:                            # both-sides needs l(F*)=1
                continue
            g = grade_both_sides(F4 + [e4, w2, w3], facs, fidx, npts)
            best = _consider(best, g, [e4, w2, w3])
            if g["both_sides_unit"]:
                return best
    return best


def _worker(job):
    edges, meta = job
    t0 = time.time()
    best = search_class(edges)
    return dict(meta=meta, edges=[list(e) for e in edges], best=best,
                secs=round(time.time() - t0, 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()

    t0 = time.time()
    classes = nonsmoothable_classes()
    assert len(classes) == 87
    if args.quick:
        classes = classes[:12]
    print(f"both-sides search over {len(classes)} non-smoothable census "
          f"classes (l(F*)=1 completions only); procs={args.procs}", flush=True)

    with Pool(processes=args.procs) as pool:
        results = pool.map(_worker, [(e, m) for e, m in classes])

    hits = [r for r in results if r["best"] and r["best"]["both_sides_unit"]]
    reached1 = [r for r in results if r["best"] and r["best"]["npoints"] == 1]
    print(f"\n{'class':40s} {'l(F*)=1?':9s} {'min max l(F*)':13s} best dual multiset")
    for r in sorted(results, key=lambda r: (r["meta"]["A2"], r["meta"]["k"])):
        m, b = r["meta"], r["best"]
        name = f"[{m['status']:8s} k={m['k']:2d} i={m['i']:2d} 2A={m['A2']:3d}]"
        if b is None:
            print(f"{name:40s} {'no l=1':9s} {'-':13s} (no single-point completion)")
            continue
        got1 = "yes" if b["npoints"] == 1 else "no"
        tag = "  <<< BOTH-SIDES UNIT" if b["both_sides_unit"] else ""
        ms = b["dual_multiset"]
        msshort = (str(ms) if len(ms) <= 10 else
                   f"[min={min(ms)},max={max(ms)},n={len(ms)}]")
        print(f"{name:40s} {got1:9s} {b['max_dual']:<13d} {msshort}{tag}", flush=True)

    print("\n" + "=" * 72)
    print(f"BOTH-SIDES-UNIT hits: {len(hits)}/{len(classes)}   "
          f"(classes reaching l(F*)=1 at all: {len(reached1)})   "
          f"[{time.time()-t0:.0f}s]")
    for r in hits:
        b = r["best"]
        print(f"  ** {r['edges']}  w={b['w']}  dual multiset {b['dual_multiset']}"
              f"  clean={b['clean']}")
    if not hits:
        # how close did we get?
        cand = [r for r in reached1]
        if cand:
            bestr = min(cand, key=lambda r: r["best"]["max_dual"])
            b = bestr["best"]
            print(f"  closest: class {bestr['edges']} reached max l(F*)="
                  f"{b['max_dual']} (dual multiset {b['dual_multiset']})")
        print("  => NO both-sides-unit polytope in the 2/3-vertex planting family.")


if __name__ == "__main__":
    main()
