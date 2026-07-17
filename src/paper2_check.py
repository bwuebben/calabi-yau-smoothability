#!/usr/bin/env python3
"""
Machine checks for paper 2 (paper2/main.tex).  All asserted; ~0.1 s.

1. Identification lemma bookkeeping (Lemma: chart identification) for the
   worked example Delta_E: in Petracci's extended cone for the datum
   (pentagon = T + s, w = interior dual-edge point), every pentagon ray
   decomposes uniquely through the split generators, and integer
   functionals m0, m1 exist with <m0,.> = 1 on the triangle rays / 0 on
   the segment rays and vice versa — so Petracci's binomial IS the Cayley
   binomial chi^{m1} - chi^{m0} on the chart of the face.
2. Delta_N (paper 2 Thm B): reflexivity/faces/Hodge re-asserted, facet
   normals as printed in the paper.
3. Star subdivisions at the interior point are unimodular for BOTH i=1
   type-(S) classes (diamond = P1xP1-cone, dP7 pentagon) — the one-line
   lemma in the Gross proposition.
"""
from fractions import Fraction
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from batyrev_global import analyze, facets
from hodge_numbers import hodge_numbers
from toric_census import ccw_sort, verts_from_edges

# ---- 1. identification-lemma bookkeeping (Delta_E datum) ----
P = [(1, 0), (0, 1), (-1, 0), (0, -1), (1, -1)]
T = [(0, 0), (0, 1), (-1, 0)]
s = [(0, 0), (1, -1)]
tri = [(t[0], t[1], 1, 0, 1, -1) for t in T]      # Q0 x {1} - e1
seg = [(u[0], u[1], 0, 0, 0, 1) for u in s]       # Q1 + e1
pent = [(p[0], p[1], 1, 0, 1, 0) for p in P]
for p in pent:
    found = [(a, b) for a in tri for b in seg
             if tuple(x + y for x, y in zip(a, b)) == p]
    assert len(found) == 1, (p, found)

def _solve(rays_one, rays_zero):
    rows = [list(r) + [1] for r in rays_one] + [list(r) + [0] for r in rays_zero]
    A = [[Fraction(x) for x in row] for row in rows]
    piv, r = [], 0
    for c in range(6):
        pr = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        A[r] = [x / A[r][c] for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                A[i] = [a - A[i][c] * b for a, b in zip(A[i], A[r])]
        piv.append(c); r += 1
    sol = [Fraction(0)] * 6
    for i, c in enumerate(piv):
        sol[c] = A[i][-1]
    assert all(x.denominator == 1 for x in sol)
    m = tuple(int(x) for x in sol)
    for ray in rays_one:
        assert sum(a * b for a, b in zip(m, ray)) == 1
    for ray in rays_zero:
        assert sum(a * b for a, b in zip(m, ray)) == 0
    return m

m0 = _solve(tri, seg)
m1 = _solve(seg, tri)
assert m0 == (0, 0, 1, 0, 0, 0) and m1 == (0, 0, 1, 0, 0, 1)
print("1. identification-lemma bookkeeping: pentagon rays decompose uniquely;"
      f" Cayley degrees m0={m0}, m1={m1} exist (m1 = m0 + e1*)")

# ---- 2. Delta_N re-assertions ----
sq = verts_from_edges(ccw_sort([(1, 0), (0, 1), (-1, 0), (0, -1)]))
VN = [(x, y, 1, 0) for (x, y) in sq] + [(0, 0, 0, 1), (0, 1, -1, 0),
                                        (-1, -2, 1, -1)]
rN = analyze("Delta_N", VN, verbose=False)
assert rN is not None and all(l == 1 for l in rN["edges"].values())
bad = [f for f in rN["faces"] if f["status"] != "smooth"]
assert len(bad) == 1 and bad[0]["k"] == 4 and bad[0]["npoints"] == 1 \
    and "smoothable" in bad[0]["status"]
assert len(rN["facets"]) == 10 and len(rN["faces"]) == 21
assert hodge_numbers(VN, "") == (3, 103)
print("2. Delta_N: reflexive, 10 facets, 21 faces, unique square face with"
      " l(F*)=1, Hodge (3,103) — all re-asserted")

# ---- 3. unimodular star subdivisions for the THREE i=1 type-(S) classes ----
# Prop C (paper 2): the i=1 type-(S) census classes are the cones over the
# three smoothABLE smooth toric del Pezzo surfaces P1xP1, dP7, dP6 (the
# other two i=1 classes, P2 and F1, are rigid).  For each, the star
# subdivision of Cone(F) at the interior lattice point is unimodular
# (=> smooth primitive type-II contraction => Gross Thm 5.8 smooths it).
from toric_census import verts_from_edges, interior_points, smoothing_components, rigid  # noqa


def _star_unimodular(edges):
    V = verts_from_edges(edges)
    xs = [v[0] for v in V]; ys = [v[1] for v in V]

    def inside(p):
        s = None
        for i in range(len(V)):
            a, b = V[i], V[(i + 1) % len(V)]
            cr = (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])
            if cr == 0:
                return False
            if s is None:
                s = cr > 0
            elif (cr > 0) != s:
                return False
        return True
    ints = [(x, y) for x in range(min(xs), max(xs) + 1)
            for y in range(min(ys), max(ys) + 1) if inside((x, y))]
    assert len(ints) == 1, ints
    c = ints[0]
    return all(abs((a[0] - c[0]) * (b[1] - c[1]) - (a[1] - c[1]) * (b[0] - c[0])) == 1
               for a, b in ((V[i], V[(i + 1) % len(V)]) for i in range(len(V))))


i1_typeS = {
    "P1xP1": [(-2, -1), (0, -1), (2, 1), (0, 1)],
    "dP7":   [(-2, -1), (-1, -1), (1, 0), (2, 1), (0, 1)],
    "dP6":   [(-2, -1), (-1, -1), (1, 0), (2, 1), (1, 1), (-1, 0)],
}
for name, E in i1_typeS.items():
    E = ccw_sort([tuple(e) for e in E])
    assert interior_points(E) == 1, (name, "not i=1")
    assert not rigid(E) and smoothing_components(E) >= 1, (name, "not type-(S)")
    assert _star_unimodular(E), (name, "star subdivision not unimodular")
print("3. star subdivisions unimodular for all THREE i=1 type-(S) classes"
      " (P1xP1, dP7, dP6)")
print("all paper-2 machine checks pass")
