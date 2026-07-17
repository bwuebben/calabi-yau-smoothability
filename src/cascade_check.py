#!/usr/bin/env python3
"""
Machine checks for the GLOBAL CASCADE theorem (Track A2, notes/PROGRAM_GROSS.md).

Claim: paper-2's boundary deformation datum + chart-identification lemma apply
verbatim to a DEF-ONLY face F = T + s (T = the non-standard 1/3(1,1,1) triangle,
s = unit segment), the ONLY new point being that the two-dimensional summand T
need not be a *standard* triangle.  The general fibre then keeps exactly one
C(T) = 1/3(1,1,1) point per point over F (paper-1 cascade lemma), so X_B deforms
to a compact CY with l(F_B*) = 3 quotient points.

This verifies the two things that could have failed for a NON-standard summand:
  (a) the 5 pentagon vertices still decompose UNIQUELY as t + s   (datum cond.);
  (b) the Cayley separating functionals m0, m1 still exist over ZZ.
plus the surrounding classification facts (T rigid non-standard, Q_B def-only,
Delta_B unit-edge reflexive with l(F_B*) = 3 >= 2).  ~0.1 s; all asserted.
"""
from fractions import Fraction
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from batyrev_global import analyze, vgcd, vsub                       # noqa: E402
from toric_census import rigid, smoothing_components                 # noqa: E402

# ---- the def-only decomposition Q_B = T + s (T = 1/3(1,1,1), s a unit segment) ----
T = [(0, 0), (1, 2), (-1, 1)]          # E(T) = {(1,2),(-2,-1),(1,-1)} : 1/3(1,1,1)
s = [(0, 0), (1, 1)]                   # unit segment, direction (1,1) (not // any T edge)

def edges_of(poly):
    k = len(poly)
    return [tuple(poly[(i + 1) % k][j] - poly[i][j] for j in (0, 1)) for i in range(k)]

ET, Es = edges_of(T), edges_of(s)
# T is a NON-standard triangle: rigid, and |det| = 3 (not unimodular)
assert rigid([tuple(e) for e in ET]) and smoothing_components([tuple(e) for e in ET]) == 0
d = abs(ET[0][0] * ET[1][1] - ET[0][1] * ET[1][0])
assert d == 3, d                       # non-standard: |det(a,b)| = 3, the 1/3(1,1,1) core

# Minkowski sum Q_B = T + s, keep the extreme (vertex) points
sumpts = sorted({(a[0] + b[0], a[1] + b[1]) for a in T for b in s})
# convex hull (monotone chain) to get the pentagon vertices
def hull(pts):
    pts = sorted(set(pts))
    def half(pp):
        h = []
        for p in pp:
            while len(h) >= 2 and (h[-1][0]-h[-2][0])*(p[1]-h[-2][1]) - (h[-1][1]-h[-2][1])*(p[0]-h[-2][0]) <= 0:
                h.pop()
            h.append(p)
        return h[:-1]
    return half(pts) + half(pts[::-1])
QB = hull(sumpts)
assert len(QB) == 5, QB                 # T + s is a pentagon
EQB = [tuple(e) for e in edges_of(QB)]
assert smoothing_components(EQB) == 0 and not rigid(EQB)   # DEF-ONLY (decomposable, 0 smoothings)
print(f"Q_B = T + s is a def-only pentagon (verts {QB}); "
      f"T is the non-standard 1/3(1,1,1) triangle (|det|={d}), s a unit segment")

# ---- (a) unique vertex decomposition in Petracci's extended cone (Cayley encoding) ----
# tri rays: T x {1}, e0-grading (1,-1);  seg rays: s, e0-grading (0,1);
# pentagon rays: (T+s) x {1}, e0-grading (1,0) = (1,-1)+(0,1).
tri = [(t[0], t[1], 1, 0, 1, -1) for t in T]
seg = [(u[0], u[1], 0, 0, 0, 1) for u in s]
pent = [(p[0], p[1], 1, 0, 1, 0) for p in QB]
for p in pent:
    found = [(a, b) for a in tri for b in seg
             if tuple(x + y for x, y in zip(a, b)) == p]
    assert len(found) == 1, (p, found)
print("(a) all 5 pentagon rays decompose UNIQUELY as (triangle ray) + (segment ray)")

# ---- (b) integer Cayley separating functionals m0, m1 exist ----
def solve_functional(rays_one, rays_zero):
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
    assert all(x.denominator == 1 for x in sol), sol   # <-- integral: the point
    m = tuple(int(x) for x in sol)
    assert all(sum(a * b for a, b in zip(m, ray)) == 1 for ray in rays_one)
    assert all(sum(a * b for a, b in zip(m, ray)) == 0 for ray in rays_zero)
    return m

m0 = solve_functional(tri, seg)
m1 = solve_functional(seg, tri)
print(f"(b) integral Cayley degrees exist: m0={m0}, m1={m1}  (m1 - m0 = {tuple(b-a for a,b in zip(m0,m1))})")

# ---- surrounding facts: Delta_B unit-edge reflexive, l(F_B*) = 3 >= 2 ----
VB = [(x, y, 1, 0) for (x, y) in [(0, 0), (-2, -1), (-3, -2), (-2, -3), (-1, -2)]] \
     + [(0, 0, 0, 1), (1, 1, -1, -1)]
rB = analyze("Delta_B", VB, verbose=False)
assert rB is not None and all(l == 1 for l in rB["edges"].values())
badB = [f for f in rB["faces"] if f["status"] != "smooth"]
assert len(badB) == 1 and badB[0]["k"] == 5 and badB[0]["status"].startswith("def-only")
lFB = badB[0]["npoints"]
assert lFB == 3, lFB
print(f"Delta_B: unit-edge reflexive; unique bad face = def-only pentagon; l(F_B*) = {lFB} >= 2")

print()
print("==> A2 verified: the paper-2 boundary datum for Q_B = T + s satisfies its")
print("    hypotheses with a NON-standard summand, so the explicit binomial family")
print("    deforms X_B (3 x C(Q_B)) to a compact CY with 3 rigid 1/3(1,1,1) points.")
print("    A global cascade: a deformation between two non-smoothable CYs.")
