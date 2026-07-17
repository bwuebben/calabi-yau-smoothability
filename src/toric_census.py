#!/usr/bin/env python3
"""
Census of isolated Gorenstein toric 3-fold singularities: rigid / deformable /
smoothable, with smoothing-component counts.

A 3-dim Gorenstein toric singularity is the cone over a lattice polygon Q (ray
generators at height 1).  It is:
  * ISOLATED   <=> every edge of Q has lattice length 1  (unit-edge polygons);
  * TERMINAL   <=> Q has no lattice points except its vertices  (i = 0 interior);
  * CANONICAL  <=> always (Gorenstein height-1);  non-terminal <=> i >= 1.

Deformation theory (verified against the literature 2026-07-10):
  * Altmann, "The versal deformation of an isolated toric Gorenstein
    singularity", Invent. Math. 128 (1997) [arXiv:alg-geom/9403004]:
    irreducible components of the reduced miniversal base Def(V) are in 1-1
    correspondence with maximal INTEGRAL Minkowski decompositions of Q.
  * Corti-Filip-Petracci, "Mirror symmetry and smoothing Gorenstein toric
    affine 3-folds" [arXiv:2006.16885], Sec. 5.1: this restricts to a 1-1
    correspondence between SMOOTHING components of Def(V) and Minkowski
    decompositions of Q whose summands are UNIT SEGMENTS and STANDARD
    TRIANGLES (= Z^2 x| GL2(Z)-equivalent to conv{(0,0),(1,0),(0,1)}).
  * Filip [arXiv:2504.04486] Thm 1.1 re-derives/extends this (0-mutable
    Laurent polynomials; general fibre smooth <=> f 0-mutable).

Hence THREE classes, not two:
    RIGID       <=> edge multiset has no proper zero-sum subset
                    (integrally Minkowski-indecomposable; no deformations)
    def-only    <=> decomposable, but NO decomposition into unit segments +
                    standard triangles  (deformable yet NON-smoothable)
    SMOOTHABLE  <=> edge multiset partitions into {v,-v} pairs and standard
                    triples {a,b,c}: a+b+c=0, |det(a,b)|=1;
                    #partitions = #smoothing components.
The day-1 census conflated the last two ("decomposable => smoothable" is FALSE:
e.g. the pentagon = (1/3)(1,1,1)-triangle + unit segment is deformable but not
smoothable — its only decomposition has a rigid non-standard triangle summand).

Equivalence of singularities = affine unimodular equivalence of Q = GL2(Z)
acting on the edge-vector multiset (translation does not move edge vectors).
Both the rigidity and smoothability tests are GL2(Z)-invariants of the edge
multiset, so those verdicts are independent of the dedup.

Validation anchors: conifold (unit square) terminal & smoothable, 1 smoothing
component; cone over P^2 = (1/3)(1,1,1) rigid; cone over dP6 (hexagon)
smoothable with exactly 2 smoothing components (3 segments / 2 triangles —
Altmann's example); the T+segment pentagon deformable-but-non-smoothable.
"""
from math import atan2, pi
from itertools import combinations

# ---------- lattice helpers ----------
def ext_gcd(a, b):
    if b == 0:
        return (abs(a), (1 if a >= 0 else -1), 0)
    g, x, y = ext_gcd(b, a % b)
    return (g, y, x - (a // b) * y)

def primitive(v):
    return v != (0, 0) and ext_gcd(v[0], v[1])[0] == 1

def ccw_sort(edges):
    return sorted(edges, key=lambda e: atan2(e[1], e[0]))

def valid_polygon(edges):
    """distinct primitive directions, sum zero, not contained in a half-plane."""
    if sum(e[0] for e in edges) != 0 or sum(e[1] for e in edges) != 0:
        return False
    if len(set(edges)) != len(edges):
        return False
    angs = sorted(atan2(e[1], e[0]) for e in edges)
    gaps = [(angs[(i + 1) % len(angs)] - angs[i]) % (2 * pi) for i in range(len(angs))]
    return max(gaps) < pi - 1e-9

def verts_from_edges(edges):
    V = [(0, 0)]
    for e in edges[:-1]:
        V.append((V[-1][0] + e[0], V[-1][1] + e[1]))
    return V

def twice_area(edges):
    V = verts_from_edges(edges)
    s = 0
    n = len(V)
    for i in range(n):
        x1, y1 = V[i]; x2, y2 = V[(i + 1) % n]
        s += x1 * y2 - x2 * y1
    return abs(s)  # = 2 * Area

def interior_points(edges):
    # Pick: A = i + b/2 - 1, b = #boundary lattice pts = #edges (unit edges).
    A2 = twice_area(edges)          # 2A
    k = len(edges)                  # b
    # i = A - k/2 + 1 = (A2 - k + 2)/2
    return (A2 - k + 2) // 2

# ---------- deformations & smoothings (Altmann / CFP 5.1) ----------
def rigid(edges):
    """No proper non-empty zero-sum subset of edge vectors
    <=> integrally Minkowski-indecomposable <=> no deformations."""
    k = len(edges)
    for r in range(1, k):
        for S in combinations(range(k), r):
            if sum(edges[i][0] for i in S) == 0 and sum(edges[i][1] for i in S) == 0:
                return False
    return True

def smoothing_components(edges):
    """#Minkowski decompositions of Q into unit segments and standard triangles
    = #smoothing components of Def(V)  (Altmann 1997 + CFP arXiv:2006.16885
    Sec. 5.1).  A decomposition = a partition of the edge multiset into parts
    {v,-v} (unit segment) or {a,b,c} with a+b+c=0, |det(a,b)|=1 (standard
    triangle).  Returns 0 <=> the singularity is NOT smoothable."""
    edges = list(edges)
    if not edges:
        return 1
    a, rest = edges[0], edges[1:]
    n = 0
    na = (-a[0], -a[1])
    if na in rest:                                   # part {a, -a}
        n += smoothing_components([e for e in rest if e != na])
    for i in range(len(rest)):                       # part {a, b, c}, pos(c)>pos(b)
        b = rest[i]
        if abs(a[0] * b[1] - a[1] * b[0]) != 1:      # symmetric in b,c since c=-a-b
            continue
        c = (-a[0] - b[0], -a[1] - b[1])
        if c in rest[i + 1:]:
            n += smoothing_components(rest[:i] + [e for e in rest[i + 1:] if e != c])
    return n

def smoothable(edges):
    return smoothing_components(edges) > 0

# ---------- dedup under GL2(Z) on edge vectors ----------
# BUGFIX (2026-07-10, day 1): a "map an edge to (1,0) then shear" canonical
# min-key had a shear runaway (x -> -inf), so its minimum was not GL2(Z)-
# invariant and it OVER-SPLIT equivalent polygons.  Since every polygon here
# has edge coords in a fixed small box, the connecting GL2(Z) map is bounded;
# union-find over "some bounded map sends P onto P'" with transitive closure
# gives genuine GL2(Z) classes (validated below against the hand-verified
# square~parallelogram equivalence, and by K-stability).  For a larger or
# fully rigorous census, defer to PALP / Sage normal_form.
def _apply(M, v):
    (a, b), (c, d) = M
    return (a * v[0] + b * v[1], c * v[0] + d * v[1])

def gl2_bounded(K):
    return [((a, b), (c, d))
            for a in range(-K, K + 1) for b in range(-K, K + 1)
            for c in range(-K, K + 1) for d in range(-K, K + 1)
            if a * d - b * c in (1, -1)]

_GL2 = gl2_bounded(3)

def equiv(P, Q, G=None):
    """True if some bounded GL2(Z) map sends P's edge multiset onto Q's."""
    G = _GL2 if G is None else G
    tq = tuple(sorted(Q))
    return any(tuple(sorted(_apply(M, e) for e in P)) == tq for M in G)

def dedup(polys, G):
    """Union-find over bounded-GL2(Z) equivalence; components = classes."""
    sig = [tuple(sorted(e)) for e in polys]
    sigset = {s: i for i, s in enumerate(sig)}          # in-box sorted-multiset -> index
    parent = list(range(len(polys)))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        parent[find(a)] = find(b)
    for i, e in enumerate(polys):
        for M in G:
            t = tuple(sorted(_apply(M, v) for v in e))
            j = sigset.get(t)                            # image landed on another in-box poly?
            if j is not None:
                union(i, j)
    reps = {}
    for i in range(len(polys)):
        reps.setdefault(find(i), polys[i])
    return reps

# ================= VALIDATION =================
def _dump(name, edges):
    sc = smoothing_components(edges)
    print(f"  {name:26s} k={len(edges)} i={interior_points(edges)} "
          f"2A={twice_area(edges)} rigid={rigid(edges)} "
          f"smoothable={sc > 0} #smoothing-comps={sc} "
          f"terminal={interior_points(edges)==0}")

def main():
    print("VALIDATION — anchors")
    conifold = ccw_sort([(1, 0), (0, 1), (-1, 0), (0, -1)])
    P2tri    = ccw_sort([(-1, 1), (-1, -2), (2, 1)])
    dP6hex   = ccw_sort([(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)])
    Tpluss   = ccw_sort([(-2, -1), (1, -1), (1, 2), (1, 0), (-1, 0)])  # (1/3)(1,1,1) + segment
    _dump("conifold (square)", conifold)     # terminal, smoothable, 1 comp
    _dump("cone/P^2 (triangle)", P2tri)      # canonical, rigid
    _dump("cone/dP6 (hexagon)", dP6hex)      # smoothable, exactly 2 comps (Altmann)
    _dump("T+segment pentagon", Tpluss)      # deformable but NOT smoothable
    assert smoothing_components(conifold) == 1
    assert rigid(P2tri) and not smoothable(P2tri)
    assert not rigid(dP6hex) and smoothing_components(dP6hex) == 2
    assert not rigid(Tpluss) and not smoothable(Tpluss), "the def-only class exists"
    # rigid k-gons exist for ALL k (kills the "k>=7 => decomposable" seed):
    # edges (1,a_1),...,(1,a_{k-1}) + primitive closing vector — any proper zero-sum
    # subset needs the closing vector plus ALL the others (x-coordinate pigeonhole).
    kgon8 = ccw_sort([(1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 7), (-7, -22)])
    assert valid_polygon(kgon8) and rigid(kgon8) and smoothing_components(kgon8) == 0
    print("  asserts: conifold 1 comp; P^2 rigid; dP6 exactly 2 comps; "
          "T+seg deformable-not-smoothable; rigid 8-gon family — all pass")

    print("\nVALIDATION — dedup (bounded GL2(Z) equivalence)")
    sym8 = [((1,0),(0,1)),((0,-1),(1,0)),((-1,0),(0,-1)),((0,1),(-1,0)),
            ((-1,0),(0,1)),((1,0),(0,-1)),((0,1),(1,0)),((0,-1),(-1,0))]
    allsym = all(equiv(ccw_sort([_apply(M, e) for e in conifold]), conifold) for M in sym8)
    para = ccw_sort([(1,0),(1,1),(-1,0),(-1,-1)])   # = [[1,1],[0,1]] . square (hand-verified)
    print(f"  8 symmetries of the square all ~ square: {allsym}  (expect True)")
    print(f"  square ~ parallelogram: {equiv(conifold, para)}  (expect True) "
          f"<- the equivalence the old min-key wrongly split")
    print(f"  square NOT ~ P^2 triangle: {not equiv(conifold, P2tri)}  (expect True)")

    # ================= CENSUS =================
    B = 2   # edge-vector coordinate box
    prims = [(x, y) for x in range(-B, B + 1) for y in range(-B, B + 1) if primitive((x, y))]

    # enumerate ALL valid unit-edge polygons in the box (day 1 stopped at k<=6;
    # the box supports up to k=16 distinct primitive directions)
    polys = []
    for k in range(3, len(prims) + 1):
        for combo in combinations(prims, k):
            edges = ccw_sort(list(combo))
            if valid_polygon(edges):
                polys.append(edges)

    n3 = len(dedup(polys, gl2_bounded(3)))
    n4 = len(dedup(polys, gl2_bounded(4)))
    print(f"\nK-STABILITY  raw-polys={len(polys)}  classes@K=3: {n3}   classes@K=4: {n4}   "
          f"stable: {n3 == n4}")

    classes = dedup(polys, gl2_bounded(4))

    # tabulate, dropping the smooth unit triangle (k=3, i=0 => C^3)
    rows = []
    for edges in classes.values():
        k = len(edges); i = interior_points(edges)
        if k == 3 and i == 0:
            continue   # smooth C^3, not a singularity
        sc = smoothing_components(edges)
        rows.append((k, i, twice_area(edges), rigid(edges), sc, edges))

    # consistency: for a singular class, rigid => not smoothable; smoothable => decomposable
    for k, i, A2, rg, sc, edges in rows:
        assert not (rg and sc > 0), f"rigid+smoothable contradiction: {edges}"

    print(f"\nCENSUS  (ALL unit-edge lattice polygons, edge coords |.|<= {B}, up to GL2(Z))")
    print(f"  distinct singularity classes: {len(rows)}")
    term = [r for r in rows if r[1] == 0]
    cano = [r for r in rows if r[1] >= 1]
    nrig = sum(1 for r in cano if r[3])
    nsm  = sum(1 for r in cano if r[4] > 0)
    ndef = len(cano) - nrig - nsm
    print(f"  TERMINAL (i=0):        {len(term)}  -> "
          f"{[(r[0], 'smoothable' if r[4] > 0 else 'NOT-smoothable') for r in term]}")
    print(f"  CANONICAL non-term:    {len(cano)}")
    print(f"    RIGID (no deformations, non-smoothable):   {nrig}")
    print(f"    DEFORMABLE but NOT smoothable (def-only):  {ndef}")
    print(f"    SMOOTHABLE:                                {nsm}")

    print("\n  canonical non-terminal classes  [k=edges, i=interior pts, "
          "#sc=#smoothing components]:")
    print("   k   i   2A  #sc  status      rep. edge vectors")
    order = {True: 0}
    for k, i, A2, rg, sc, edges in sorted(
            cano, key=lambda r: (r[0], r[1], 0 if r[3] else (1 if r[4] == 0 else 2))):
        tag = "RIGID     " if rg else ("def-only  " if sc == 0 else "SMOOTHABLE")
        print(f"  {k:2d}  {i:2d}  {A2:3d}  {sc:3d}  {tag}  {edges}")

    nontri_rigid = sorted([r for r in cano if r[3] and r[0] > 3], key=lambda r: (r[2], r[0]))
    if nontri_rigid:
        r = nontri_rigid[0]
        print(f"\n  smallest RIGID non-triangle: k={r[0]}, i={r[1]}, edges={r[5]}")
    defonly = sorted([r for r in cano if not r[3] and r[4] == 0], key=lambda r: (r[2], r[0]))
    if defonly:
        r = defonly[0]
        print(f"  smallest DEFORMABLE-but-NON-SMOOTHABLE: k={r[0]}, i={r[1]}, "
              f"2A={r[2]}, edges={r[5]}")


if __name__ == "__main__":
    main()
