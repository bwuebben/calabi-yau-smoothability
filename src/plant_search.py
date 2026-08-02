#!/usr/bin/env python3
"""
Planting search at scale: which of the 87 non-smoothable census classes occur
as 2-faces of (unit-edge) reflexive 4-polytopes?  (Paper Sec. 6, Q1)  And can
the dual edge have lattice length 1, i.e. a compact CY with a SINGLE
non-smoothable point?  (Paper Sec. 6, Q2)

Gauge (RESEARCH_LOG day 2): embed the polygon Q in the plane {x3=1, x4=0}
(saturated, so the induced lattice is intrinsic); in-plane translations are
shears fixing 0, and a completing vertex with x4 = +-1 can be moved to e4 by a
unimodular map fixing the plane pointwise.  w1 = e4 is WLOG for two-vertex
completions (up to completions whose extra vertices all have |x4| >= 2), and
an ansatz for three-vertex completions.

STRUCTURE THEOREM FOR PLANTED CANDIDATES (used heavily; proof in comments):
Let V = F x {(1,0)} u extras.  A facet of conv(V) either
  (a) contains the whole plane aff(F): its normal is u = (0,0,1,t), t in Z
      (the unique functionals that restrict to 1 on that plane).  On an extra
      vertex w = (a,b,c,d) the value is c + t d, so u is supporting iff
      c + t d <= 1 for every extra; the set of such integer t is an interval
      [lo, hi], and (0,0,1,t) is a FACET iff some extra achieves equality
      (the face F u {equality extras} is then 3-dimensional).  Since no extra
      of ours has (c,d) = (1,0), equality at a given extra pins t, so the
      pencil facets are exactly t = lo and t = hi when those bounds are
      achieved.  Consequently:
        * F is a 2-face  <=>  lo < hi and both bounds achieved (then the
          2-face is exactly conv(F): the only V-points on both facets are F's);
        * the dual edge F* = conv{(0,0,1,lo),(0,0,1,hi)} has lattice length
          hi - lo  =  the number of C(F)-points on X.  No facet computation
          needed for any of this — O(#extras) per candidate.
  (b) does not contain aff(F): it meets that plane in at most a line, hence
      contains at most 2 vertices of F (a line meets the vertex set of a
      strictly convex polygon in at most 2 points).  With <= 3 extras, every
      such facet is spanned by a 4-subset with <= 2 F-vertices — so it
      suffices to enumerate those subsets: O(k^2) instead of O(k^4).
Reflexivity = every facet at level 1 with primitive normal.  Pencil facets are
at level 1 automatically; the (b)-facets are checked (early abort on the first
supporting functional at level != 1, which also catches 0-not-interior).
Both pencil facets at level 1 and a complete level-1 facet list imply 0 is
interior; full-dimensionality holds identically in our families (e4 and a
vertex with x3 < 0 are always present — see full_dim_note in code).

Dual-edge arithmetic in this gauge (paper Sec. 6 Q2): with extras {e4, w2},
hi = 1 (achieved by e4) and lo <= -1 (a t = 0 facet needs a second vertex at
x3 = 1 off the plane), so l(F*) >= 2 — the day-2 family can never give a
single point.  Adding w3 = (c1, c2, 1, -r), r >= 1, forces lo = 0 (achieved by
w3) and hi = 1: l(F*) = 1 automatically.  A CLEAN such hit (all edges unit,
every other 2-face smooth) is a compact CY threefold whose singular locus is
a SINGLE non-smoothable point.

Run:  python3 -u src/plant_search.py [--quick] [--json PATH] [--procs N]
"""
from itertools import combinations
from multiprocessing import Pool
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from toric_census import (ccw_sort, valid_polygon, verts_from_edges,        # noqa: E402
                          primitive, rigid, smoothable, interior_points,
                          twice_area, dedup, gl2_bounded)
from batyrev_global import (dot, vsub, vgcd, hyperplane_normal,             # noqa: E402
                            two_faces, face_lattice_polygon, classify_polygon)


# ---------------- census: the 87 non-smoothable classes ----------------
def nonsmoothable_classes(B=2):
    prims = [(x, y) for x in range(-B, B + 1) for y in range(-B, B + 1)
             if primitive((x, y))]
    polys = []
    for k in range(3, len(prims) + 1):
        for combo in combinations(prims, k):
            edges = ccw_sort(list(combo))
            if valid_polygon(edges):
                polys.append(edges)
    classes = dedup(polys, gl2_bounded(4)).values()
    out = []
    for edges in classes:
        k, A2 = len(edges), twice_area(edges)
        if k == 3 and A2 == 1:
            continue
        if smoothable(edges):
            continue
        out.append((edges, dict(k=k, i=interior_points(edges), A2=A2,
                                status="rigid" if rigid(edges) else "def-only")))
    out.sort(key=lambda t: (t[1]["A2"], t[1]["k"]))
    return out


def embed(edges):
    return [(x, y, 1, 0) for (x, y) in verts_from_edges(edges)]


# ---------------- the fast planted-candidate check ----------------
def pencil_interval(extras):
    """Integer interval [lo, hi] of supporting pencil functionals (0,0,1,t),
    with flags whether each bound is achieved by an extra vertex.
    Returns (lo, hi, lo_achieved, hi_achieved) or None (no valid t)."""
    lo, hi = None, None
    for (_, _, c, d) in extras:
        if d > 0:
            t = (1 - c) // d                     # floor
            hi = t if hi is None else min(hi, t)
        elif d < 0:
            t = -((1 - c) // (-d))               # ceil((1-c)/d), d<0
            lo = t if lo is None else max(lo, t)
        else:
            if c > 1:
                return None
    if lo is None or hi is None or lo > hi:
        return None                              # unbounded or empty
    lo_ach = any(d and c + lo * d == 1 for (_, _, c, d) in extras)
    hi_ach = any(d and c + hi * d == 1 for (_, _, c, d) in extras)
    return lo, hi, lo_ach, hi_ach


def planted_facets(F4, extras):
    """Complete facet list of conv(F4 + extras) IF it is reflexive with F a
    2-face; else None.  Uses the structure theorem in the module docstring.
    Returns (facets, npoints) with facets = [(u, 1, idx-frozenset)]."""
    pi = pencil_interval(extras)
    if pi is None:
        return None
    lo, hi, lo_ach, hi_ach = pi
    if lo >= hi or not (lo_ach and hi_ach):
        return None                              # F not a 2-face
    V = F4 + extras
    nf, ne = len(F4), len(extras)
    facs = []
    for t in (lo, hi):
        idx = frozenset(list(range(nf)) +
                        [nf + j for j, (_, _, c, d) in enumerate(extras)
                         if c + t * d == 1])
        facs.append(((0, 0, 1, t), 1, idx))
    # (b) non-pencil facets: 4-subsets with <= 2 F-vertices
    seen = {(0, 0, 1, lo), (0, 0, 1, hi)}
    fidx_range = range(nf)
    eidx_range = range(nf, nf + ne)
    subsets = []
    if ne >= 2:
        for fpair in combinations(fidx_range, 2):
            for epair in combinations(eidx_range, 2):
                subsets.append(fpair + epair)
    if ne >= 3:
        for f1 in fidx_range:
            for etri in combinations(eidx_range, 3):
                subsets.append((f1,) + etri)
        if ne >= 4:
            for equad in combinations(eidx_range, 4):
                subsets.append(equad)
    for S in subsets:
        u = hyperplane_normal(*[V[i] for i in S])
        if u == (0, 0, 0, 0):
            continue
        c = dot(u, V[S[0]])
        vals = [dot(u, v) for v in V]
        if all(x <= c for x in vals):
            pass
        elif all(x >= c for x in vals):
            u = tuple(-x for x in u); c = -c; vals = [-x for x in vals]
        else:
            continue
        g = vgcd(u)
        if c != g:
            return None                          # facet not at level 1
        u = tuple(x // g for x in u)
        if u in seen:
            continue
        seen.add(u)
        facs.append((u, 1, frozenset(i for i, x in enumerate(vals) if x == g)))
    return facs, hi - lo


def grade_hit(V, facs, want_face_idx, npoints):
    """all edges unit?  all other 2-faces smooth?"""
    tf = two_faces(V, facs)
    all_unit, others_smooth = True, True
    for I, fpair in tf:
        u1, u2 = facs[fpair[0]][0], facs[fpair[1]][0]
        _, evs, lens = face_lattice_polygon(V, I, u1, u2)
        if any(l >= 2 for l in lens):
            all_unit = False
        if I != want_face_idx and classify_polygon(evs, lens)["status"] != "smooth":
            others_smooth = False
    return dict(npoints=npoints, all_unit=all_unit, others_smooth=others_smooth,
                clean=(all_unit and others_smooth))


# ---------------- per-class searches ----------------
def two_vertex_search(edges, B, QMAX):
    """extras = {e4, w2 = (b1,b2,b3,-q)}, b3 < 0, 1 <= q <= QMAX.
    (Full-dim: F spans the plane, e4 adds x4, w2 has x3 - 1 = b3 - 1 < 0 and
    the 2x2 block [(-1,1),(b3-1,-q)] has det q+1-b3 > 0 — always rank 4.)"""
    F4 = embed(edges)
    fidx = frozenset(range(len(F4)))
    e4 = (0, 0, 0, 1)
    best = None
    for b1 in range(-B, B + 1):
        for b2 in range(-B, B + 1):
            for b3 in range(-B, 0):
                for q in range(1, QMAX + 1):
                    w2 = (b1, b2, b3, -q)
                    if vgcd(w2) != 1:
                        continue
                    r = planted_facets(F4, [e4, w2])
                    if r is None:
                        continue
                    facs, npts = r
                    g = grade_hit(F4 + [e4, w2], facs, fidx, npts)
                    g["w"] = [list(e4), list(w2)]
                    if best is None or (g["clean"], -g["npoints"]) > (best["clean"], -best["npoints"]):
                        best = g
                    if g["clean"] and g["npoints"] == 2:
                        return best              # 2 is the family minimum
    return best


def single_point_search(edges, B, QMAX, RMAX):
    """extras = {e4, w2 = (b1,b2,b3,e), b3 < 0, w3 = (c1,c2,1,-r), r >= 1}:
    the l(F*) = 1 design (lo = 0 via w3, hi = 1 via e4 — automatic)."""
    F4 = embed(edges)
    fidx = frozenset(range(len(F4)))
    e4 = (0, 0, 0, 1)
    best = None
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
                                    continue     # cannot happen by design
                                g = grade_hit(F4 + [e4, w2, w3], facs, fidx, npts)
                                g["w"] = [list(e4), list(w2), list(w3)]
                                if best is None or g["clean"] > best["clean"]:
                                    best = g
                                if g["clean"]:
                                    return best
    return best


def three_vertex_search(edges, B=2):
    """General completions {e4, w2, w3} with w2, w3 primitive in [-B,B]^4.
    Exclusions: in-plane vectors (x3=1, x4=0) would break the pencil argument
    (they satisfy (c,d)=(1,0)); need x4 < 0 and x3 < 0 somewhere; full-dim
    <=> some extra off the hyperplane x3 + x4 = 1 (= aff span of F u {e4})."""
    from itertools import product
    F4 = embed(edges)
    fidx = frozenset(range(len(F4)))
    e4 = (0, 0, 0, 1)
    W = [w for w in product(range(-B, B + 1), repeat=4)
         if w != (0, 0, 0, 0) and vgcd(w) == 1
         and not (w[2] == 1 and w[3] == 0) and w != e4]
    best = None
    for i, w2 in enumerate(W):
        for w3 in W[i + 1:]:
            if w2[3] >= 0 and w3[3] >= 0:
                continue
            if w2[2] >= 0 and w3[2] >= 0:
                continue
            if w2[2] + w2[3] == 1 and w3[2] + w3[3] == 1:
                continue                         # degenerate (not full-dim)
            r = planted_facets(F4, [e4, w2, w3])
            if r is None:
                continue
            facs, npts = r
            g = grade_hit(F4 + [e4, w2, w3], facs, fidx, npts)
            g["w"] = [list(e4), list(w2), list(w3)]
            if best is None or (g["clean"], -g["npoints"]) > (best["clean"], -best["npoints"]):
                best = g
            if g["clean"]:
                return best
    return best


def worker3(job):
    edges, meta = job
    t0 = time.time()
    hit = three_vertex_search(edges, B=2)
    return dict(meta=meta, edges=[list(e) for e in edges], hit=hit,
                secs=round(time.time() - t0, 1))


def worker(job):
    """Per-class pipeline: two-vertex boxes escalating, then the l=1 design."""
    edges, meta, quick = job
    res = dict(meta=meta, edges=[list(e) for e in edges])
    t0 = time.time()
    hit = two_vertex_search(edges, B=3, QMAX=3)
    res["A"] = hit
    stageB = None
    if hit is None and not quick:
        for (B, Q) in [(4, 4), (5, 5)]:
            stageB = two_vertex_search(edges, B=B, QMAX=Q)
            if stageB is not None:
                res["Bbox"] = [B, Q]
                break
    res["B"] = stageB
    res["L1"] = single_point_search(edges, B=2, QMAX=2, RMAX=2)
    if res["L1"] is None and meta["k"] <= 8 and not quick:
        res["L1"] = single_point_search(edges, B=3, QMAX=2, RMAX=2)
        if res["L1"] is not None:
            res["L1box"] = 3
    res["secs"] = round(time.time() - t0, 1)
    return res


def selftest():
    """Fast path == day-2 exact path on the two headline polytopes."""
    from batyrev_global import facets as slow_facets, is_reflexive
    for name, Q in [("Delta_A", [(-2, -1), (0, -1), (1, 0), (1, 2)]),
                    ("Delta_B", [(-2, -1), (-1, -1), (1, -1), (1, 1), (1, 2)])]:
        edges = ccw_sort(Q)
        F4 = embed(edges)
        extras = [(0, 0, 0, 1), (1, 1, -1, -1)]
        r = planted_facets(F4, extras)
        assert r is not None, f"{name}: fast path rejected a known hit"
        facs, npts = r
        assert npts == 3, f"{name}: expected 3 points, got {npts}"
        sf = slow_facets(F4 + extras)
        assert is_reflexive(sf)
        assert sorted(u for u, _, _ in facs) == sorted(u for u, _, _ in sf), \
            f"{name}: facet sets differ"
        g = grade_hit(F4 + extras, facs, frozenset(range(len(F4))), npts)
        assert g["clean"], f"{name}: headline hit should be clean"
    # negative control: quintic-style completion that is NOT reflexive with F planted
    assert planted_facets(embed(ccw_sort([(-2, -1), (0, -1), (1, 0), (1, 2)])),
                          [(0, 0, 0, 1), (3, 3, -1, -1)]) is None
    print("selftest: fast planted-facet path matches exact facets on "
          "Delta_A/Delta_B (clean, 3 pts each) + negative control — OK", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--stage3", action="store_true",
                    help="general 3-vertex completions on the classes the "
                         "2-vertex family misses")
    ap.add_argument("--json", default=None)
    ap.add_argument("--procs", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    args = ap.parse_args()

    t0 = time.time()
    selftest()

    if args.stage3:
        classes = nonsmoothable_classes()
        with Pool(processes=args.procs) as pool:
            planted2 = pool.starmap(two_vertex_search,
                                    [(e, 3, 3) for e, _ in classes])
        todo = [(e, m) for (e, m), p in zip(classes, planted2) if p is None]
        print(f"stage 3: general three-vertex completions on the "
              f"{len(todo)} classes the two-vertex family misses", flush=True)
        with Pool(processes=args.procs) as pool:
            res3 = pool.map(worker3, todo)
        newly = [r for r in res3 if r["hit"]]
        for r in res3:
            m, h = r["meta"], r["hit"]
            tag = (("CLEAN" if h["clean"] else "hit") + f" pts={h['npoints']}"
                   + f" w={h['w'][1:]}") if h else "-"
            print(f"  [{m['status']:8s} k={m['k']:2d} i={m['i']:2d} "
                  f"2A={m['A2']:3d}] {tag}   ({r['secs']}s)", flush=True)
        print(f"\nstage 3: {len(newly)}/{len(todo)} newly planted "
              f"({sum(1 for r in newly if r['hit']['clean'])} clean)   "
              f"[{time.time()-t0:.0f}s]")
        if args.json:
            with open(args.json, "w") as f:
                json.dump(res3, f, indent=1, default=str)
            print(f"results -> {args.json}")
        return
    classes = nonsmoothable_classes()
    assert len(classes) == 87, f"expected 87 non-smoothable classes, got {len(classes)}"
    print(f"census: 87 non-smoothable classes "
          f"({sum(1 for _, m in classes if m['status'] == 'rigid')} rigid, "
          f"{sum(1 for _, m in classes if m['status'] == 'def-only')} def-only); "
          f"procs={args.procs}", flush=True)
    if args.quick:
        classes = classes[:12]

    jobs = [(e, m, args.quick) for e, m in classes]
    with Pool(processes=args.procs) as pool:
        results = pool.map(worker, jobs)

    planted = clean = single = single_clean = 0
    print(f"\n{'class':44s} {'verdict':9s}  2-vertex           l=1 design")
    for r in results:
        m = r["meta"]
        hit = r["A"] or r["B"]
        l1 = r["L1"]
        if hit:
            planted += 1
            clean += hit["clean"]
        if l1:
            single += 1
            single_clean += l1["clean"]
        htag = ("CLEAN" if hit["clean"] else "hit  ") + f" pts={hit['npoints']}" if hit else "-       "
        ltag = ("CLEAN SINGLE PT" if l1["clean"] else "dirty single pt") if l1 else "-"
        name = f"[{m['status']:8s} k={m['k']:2d} i={m['i']:2d} 2A={m['A2']:3d}]"
        print(f"{name:44s} {htag:18s} {ltag}   ({r['secs']}s)", flush=True)

    print("\n" + "=" * 72)
    print(f"PLANTED: {planted}/{len(classes)} (clean {clean})    "
          f"SINGLE-POINT: {single} (clean {single_clean})    "
          f"[{time.time()-t0:.0f}s]")
    for r in results:
        l1 = r["L1"]
        if l1 and l1["clean"]:
            print(f"  ** CLEAN SINGLE-POINT: {r['edges']}  w={l1['w']}")

    if args.json:
        with open(args.json, "w") as f:
            json.dump(results, f, indent=1, default=str)
        print(f"results -> {args.json}")


if __name__ == "__main__":
    main()
