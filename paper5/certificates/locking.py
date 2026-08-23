#!/usr/bin/env python3
"""The rigidity criterion of Theorem A', and the proof that it is
degree-independent.

Theorem A on Delta_19 says: no single-degree Ilten-Vollmert family deforms
the dP7 germ, at any primitive degree.  Its case (C) rested on a numerical
fact -- a pentagon-bearing facet cell has one-dimensional summand cone --
re-established at every degree by a computation over Q(t).  That is a
statement about one polytope.

This file replaces it by a COMBINATORIAL criterion on the facet, and proves
the criterion transfers to every degree.  The point is the following.  For
a facet F of Delta and a degree R positive on every ray of sigma_F, the cell
        C_R = sigma_F cap [R = 1] = conv{ v / <R,v> : v in F }
is a cross-section of a pointed cone by a transversal hyperplane, hence a
3-polytope COMBINATORIALLY ISOMORPHIC to F.  Its edge vectors are not those
of F -- they are rescaled -- but the three forcing rules below use only two
facts about a 2-face, and both survive the rescaling:

    (a) the edge vectors around a 2-face sum to zero;
    (b) two consecutive edges of a convex polygon are linearly independent.

The rules, with their one-line proofs.  Let a 2-face have edge vectors
e_1..e_k in cyclic order, sum zero, and dilations t_1..t_k; the closing
condition is sum t_i e_i = 0.

  TRIANGLE   k = 3, nothing known.  Substituting e_3 = -e_1-e_2 gives
             (t_1-t_3)e_1 + (t_2-t_3)e_2 = 0, and e_1,e_2 are independent
             by (b), so t_1 = t_2 = t_3.

  SINGLE     t_i = L for every i except one, say k.  Then
             L(-e_k) + t_k e_k = 0 by (a), so t_k = L since e_k =/= 0.

  ADJACENT   t_i = L except for two CONSECUTIVE edges e_1,e_2.  Then
             (t_1-L)e_1 + (t_2-L)e_2 = 0 by (a), so t_1 = t_2 = L by (b).

Neither (a) nor (b) mentions the lattice or the degree, so a facet on which
the chain closes is rigid at EVERY degree in the interior of sigma_F^dual.
The 'adjacent' hypothesis is not decoration: two OPPOSITE edges of a
quadrilateral can be parallel, and then the rule is false -- the square
decomposes as segment + segment.  That is checked below.

Run:  python3 locking.py     (from certificates/)
"""
import itertools, json, os, sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper4", "certificates"))
from batyrev_global import facets, two_faces          # noqa: E402

CH = [0]
def ok(label, cond):
    CH[0] += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {label}")
    assert cond, label


# ----------------------------------------------------------------- combinatorics
def facet_complex(V, facs, tf, fi):
    """(vertices, 2-faces, edges, cyclic orders) of the 3-polytope conv(F_fi)."""
    twofaces = [sorted(I) for I, fp in tf if fi in fp]
    edges = set()
    for A, B in itertools.combinations(twofaces, 2):
        common = sorted(set(A) & set(B))
        if len(common) == 2:
            edges.add(tuple(common))
    cycles = []
    for A in twofaces:
        adj = {}
        for e in edges:
            if set(e) <= set(A):
                adj.setdefault(e[0], []).append(e[1])
                adj.setdefault(e[1], []).append(e[0])
        if not all(len(adj.get(v, [])) == 2 for v in A):
            return None
        cyc, prev, cur = [A[0]], None, A[0]
        while len(cyc) < len(A):
            nxt = [x for x in adj[cur] if x != prev][0]
            cyc.append(nxt); prev, cur = cur, nxt
        cycles.append(cyc)
    return sorted(facs[fi][2]), twofaces, sorted(edges), cycles


def lock(edges, cycles, seed_edges=None):
    """Run the three rules.  Returns (class of the seed, number of classes).

    The chain needs a seed, and the only rule that fires with nothing known
    is TRIANGLE, so the seed is the class the triangles generate.  With no
    triangle at all the chain cannot start and we report no locking.
    """
    par = {e: e for e in edges}
    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: par[ra] = rb

    seeded = False
    for cyc in cycles:
        if len(cyc) == 3:
            es = [tuple(sorted((cyc[j], cyc[(j + 1) % 3]))) for j in range(3)]
            union(es[0], es[1]); union(es[1], es[2])
            if not seeded: seed, seeded = es[0], True
    if seed_edges:
        for e in seed_edges[1:]: union(seed_edges[0], e)
        if not seeded: seed, seeded = seed_edges[0], True
    if not seeded:
        return set(), len({find(e) for e in edges})

    big = find(seed)
    progress = True
    while progress:
        progress = False
        for cyc in cycles:
            fe = [tuple(sorted((cyc[j], cyc[(j + 1) % len(cyc)])))
                  for j in range(len(cyc))]
            rest = [e for e in fe if find(e) != big]
            if not rest:
                continue
            if len(rest) == 1:
                union(rest[0], big); big = find(big); progress = True
            elif len(rest) == 2 and len(set(rest[0]) & set(rest[1])) == 1:
                union(rest[0], big); union(rest[1], big)
                big = find(big); progress = True
    big = find(seed)
    return {e for e in edges if find(e) == big}, len({find(e) for e in edges})


def is_locked(V, facs, tf, fi):
    fc = facet_complex(V, facs, tf, fi)
    if fc is None:
        return None
    _, _, edges, cycles = fc
    forced, ncl = lock(edges, cycles)
    return len(forced) == len(edges) and ncl == 1


# ----------------------------------------------------------------- linear algebra
def summand_dim(points, edges, cycles):
    """Dimension of the weak Minkowski summand cone, by the edge-cycle system.

    points: dict vertex-index -> tuple of Fractions.  One variable per edge,
    one closing condition per 2-face and per ambient coordinate.
    """
    eidx = {e: k for k, e in enumerate(edges)}
    n, rows = len(edges), []
    dim = len(next(iter(points.values())))
    for cyc in cycles:
        for c in range(dim):
            row = [Fraction(0)] * n
            for j in range(len(cyc)):
                a, b = cyc[j], cyc[(j + 1) % len(cyc)]
                e = tuple(sorted((a, b)))
                row[eidx[e]] += points[b][c] - points[a][c]
            rows.append(row)
    # rank by fraction-free elimination
    r, piv = 0, 0
    while r < len(rows) and piv < n:
        p = next((i for i in range(r, len(rows)) if rows[i][piv] != 0), None)
        if p is None:
            piv += 1; continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(r + 1, len(rows)):
            if rows[i][piv] != 0:
                f = rows[i][piv] / rows[r][piv]
                rows[i] = [x - f * y for x, y in zip(rows[i], rows[r])]
        r += 1; piv += 1
    return n - r


def rescale(V, idx, weights):
    return {i: tuple(Fraction(c, 1) / weights[i] for c in V[i]) for i in idx}


# ----------------------------------------------------------------- self-tests
print("== the rules are sound: validation on polytopes with known summand cones ==")
CASES = [
    ("tetrahedron, indecomposable", 1,
     {0: (0,0,0), 1: (1,0,0), 2: (0,1,0), 3: (0,0,1)},
     [[0,1,2],[0,1,3],[0,2,3],[1,2,3]]),
    ("triangular prism = triangle + segment", 2,
     {0:(0,0,0),1:(1,0,0),2:(0,1,0),3:(0,0,1),4:(1,0,1),5:(0,1,1)},
     [[0,1,2],[3,4,5],[0,1,4,3],[1,2,5,4],[0,2,5,3]]),
    ("cube = three segments", 3,
     {k:v for k,v in enumerate(
        [(x,y,z) for x in (0,1) for y in (0,1) for z in (0,1)])},
     [[0,1,3,2],[4,5,7,6],[0,1,5,4],[2,3,7,6],[0,2,6,4],[1,3,7,5]]),
    ("square pyramid", 1,
     {0:(0,0,0),1:(1,0,0),2:(1,1,0),3:(0,1,0),4:(0,0,1)},
     [[0,1,2,3],[0,1,4],[1,2,4],[2,3,4],[0,3,4]]),
    ("pentagonal pyramid, the configuration of Delta_9's facet 8", 1,
     {0:(1,0,0),1:(0,1,0),2:(-1,1,0),3:(-1,0,0),4:(0,-1,0),5:(0,0,1)},
     [[0,1,2,3,4],[0,1,5],[1,2,5],[2,3,5],[3,4,5],[4,0,5]]),
]
verdicts = []
for name, want, pts, faces in CASES:
    edges = sorted({tuple(sorted(set(A) & set(B)))
                    for A, B in itertools.combinations(faces, 2)
                    if len(set(A) & set(B)) == 2})
    pts = {k: tuple(Fraction(c) for c in v) for k, v in pts.items()}
    d = summand_dim(pts, edges, faces)
    forced, ncl = lock(edges, faces)
    locked = len(forced) == len(edges) and ncl == 1
    verdicts.append((name, d, locked))
    ok(f"{name}: summand dim {d}, want {want}", d == want)
ok("no locked polytope in the list is decomposable, which is the direction "
   "the criterion claims: " + "; ".join(f"{n} dim {d} locked {l}"
                                        for n, d, l in verdicts),
   all(d == 1 for _, d, l in verdicts if l))

print("\n== the ADJACENT hypothesis is not decoration ==")
pts = {k: tuple(Fraction(c) for c in v) for k, v in enumerate(
    [(x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1)])}
faces = [[0,1,3,2],[4,5,7,6],[0,1,5,4],[2,3,7,6],[0,2,6,4],[1,3,7,5]]
edges = sorted({tuple(sorted(set(A) & set(B)))
                for A, B in itertools.combinations(faces, 2)
                if len(set(A) & set(B)) == 2})
ok("a quadrilateral can have two OPPOSITE edges parallel, and then the two "
   "leftover dilations are not separated; the cube is the witness, with "
   f"summand dim {summand_dim(pts, edges, faces)}, so a rule firing on "
   "opposite pairs would be unsound", summand_dim(pts, edges, faces) == 3)

print("\n== DEGREE-INDEPENDENCE: locking survives every rescaling ==")
sys.path.insert(0, HERE)
from examples import V_19                                            # noqa: E402
V = [tuple(v) for v in V_19]
facs = facets(V_19)
tf = two_faces(V_19, facs)
PENT = frozenset([14, 15, 16, 17, 18])
pentfacets = sorted({f for I, fp in tf if I == PENT for f in fp})
ok(f"the pentagon of Delta_19 lies in exactly two facets, {pentfacets}",
   len(pentfacets) == 2)

WEIGHTS = [
    ("lattice facet",   lambda i: Fraction(1)),
    ("uniform 3",       lambda i: Fraction(3)),
    ("index+1",         lambda i: Fraction(i + 1)),
    ("primes",          lambda i: Fraction([2,3,5,7,11,13,17,19,23,29,31,37,
                                            41,43,47,53,59,61,67][i % 19])),
    ("lopsided",        lambda i: Fraction(1 if i % 3 else 1000)),
    ("reciprocals",     lambda i: Fraction(1, i + 2)),
]
for fi in pentfacets:
    idx, tfs, edges, cycles = facet_complex(V, facs, tf, fi)
    forced, ncl = lock(edges, cycles)
    locked = len(forced) == len(edges) and ncl == 1
    ok(f"facet {fi}: {len(idx)} vertices, {len(edges)} edges, 2-face types "
       f"{sorted(len(c) for c in cycles)}; LOCKED", locked)
    dims = [(nm, summand_dim(rescale(V, idx, {i: w(i) for i in idx}),
                             edges, cycles)) for nm, w in WEIGHTS]
    ok(f"facet {fi}: summand dim stays 1 under all {len(WEIGHTS)} rescalings "
       f"({', '.join(f'{n} {d}' for n, d in dims)})",
       all(d == 1 for _, d in dims))

print("\n== the CONTROL: an UNLOCKED facet's dimension does move ==")
moved, unlocked = [], 0
for fi in range(len(facs)):
    fc = facet_complex(V, facs, tf, fi)
    if fc is None:
        continue
    idx, tfs, edges, cycles = fc
    forced, ncl = lock(edges, cycles)
    if len(forced) == len(edges) and ncl == 1:
        continue
    unlocked += 1
    base = summand_dim(rescale(V, idx, {i: Fraction(1) for i in idx}),
                       edges, cycles)
    for nm, w in WEIGHTS[1:]:
        d = summand_dim(rescale(V, idx, {i: w(i) for i in idx}), edges, cycles)
        if d != base:
            moved.append((fi, base, nm, d))
            break
for fi, base, nm, d in moved[:4]:
    print(f"    facet {fi:>2}: summand dim {base} on the lattice facet, "
          f"{d} under '{nm}'")
ok(f"of the {len(facs)} facets of Delta_19, {unlocked} are unlocked and "
   f"{len(moved)} of those change summand dimension under rescaling.  So the "
   "numerical dimension is NOT a degree-independent invariant, and the "
   "criterion has to be combinatorial", len(moved) > 0)

print(f"\n{CH[0]} checks passed.")
print("""
CONCLUSION.  Locking is a property of the FACE LATTICE of a facet together
with which 2-faces are triangles.  It implies Minkowski rigidity of the cell
at every degree in the interior of sigma_F^dual, because the three rules use
only that edge vectors close up around a 2-face and that consecutive edges of
a convex polygon are independent.  It is sufficient and not necessary.  And it
is the right notion rather than the numerical summand dimension, which the
control above shows is degree-dependent.""")
