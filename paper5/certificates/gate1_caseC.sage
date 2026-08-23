#!/usr/bin/env sage
"""The case-(C) lines of gate 1, settled for EVERY t rather than over a range.

gate1_final.sage reduces gate 1 to thirty lines in M.  Twenty-two die
wholesale: the crosscut depends only on the pairing vector, which is constant
along a line, and there it admits no decomposition into unimodular simplices.
Eight lines are left -- four pairing vectors and their negatives.  On those the
crosscut could be smoothed, and what stops it is that a pentagon-bearing facet
cell is Minkowski-rigid, which is NOT a function of the pairing vector alone
and so has to be established at every t.

Method, the one slice_rigidity.sage uses on the (1,1,1,1,1) line: check that
the cell's combinatorial type does not move over a wide sample of t, then
compute the summand dimension once, symbolically, over the function field Q(t)
with that type fixed.

What makes it close is that the two pentagon-bearing facets have COMPLEMENTARY
ranges of boundedness -- facet 21 for t > 0, facet 24 for t < k with k >= 1 --
so between them they cover the whole line, and neither unbounded regime is ever
needed.

The negative pairing vectors need no separate work: R and -R give literally the
same cell, since the level is -1 for the second where it is +1 for the first,
and v/(-<R,v>) = v/<-R,v>.  That identity is checked below.

Run:  sage gate1_caseC.sage     (from paper5/)
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper4", "certificates"))
from batyrev_global import facets
from examples import V_19
exec(open(os.path.join(HERE, "slice_rigidity.sage")).read().split("print(\"== routine validated")[0]
     .split("from examples import V_19")[1])

LINES = [(1,1,1,1,1), (2,2,2,2,2), (3,3,3,3,3), (1,2,1,2,3)]
Asys = matrix(ZZ, [list(V[i]) for i in PENT])
dv = vector(ZZ, [ZZ(x) for x in Asys.right_kernel().basis()[0]])

def bounded_range(R0, F):
    """the open interval of t on which every ray of F pairs positively"""
    lo, hi = -Infinity, Infinity
    for i in F:
        a = R0.dot_product(vector(QQ, V[i]))
        b = vector(QQ, dv).dot_product(vector(QQ, V[i]))
        if b == 0:
            if a <= 0: return None
        elif b > 0: lo = max(lo, -a / b)
        else:       hi = min(hi, -a / b)
    return None if lo >= hi else (lo, hi)

def sample(lo, hi, n=10):
    """n rational points spread across (lo, hi), including toward the ends"""
    if lo == -Infinity and hi == Infinity: base = [QQ(k) for k in range(-n, n)]
    elif lo == -Infinity: base = [hi - QQ(2)**k for k in range(-4, n)]
    elif hi == Infinity:  base = [lo + QQ(2)**k for k in range(-4, n)]
    else: base = [lo + (hi - lo) * QQ((j, n + 1)) for j in range(1, n + 1)]
    return [x for x in base if lo < x < hi]

print("== R and -R give the same cell, so only the positive vectors need work ==")
same = True
for pv in LINES:
    Rp = vector(QQ, Asys.solve_right(vector(QQ, pv)))
    Rm = -Rp
    for i in PENT:
        same &= (vector(QQ, V[i]) / Rp.dot_product(vector(QQ, V[i]))
                 == vector(QQ, V[i]) / (-Rm.dot_product(vector(QQ, V[i]))))
ok("for each of the four vectors, the level-(+1) cell of R equals the "
   "level-(-1) cell of -R vertex for vertex", same)

print("\n== the four case-(C) lines, each rigid at every t ==")
for pv in LINES:
    R0 = vector(QQ, Asys.solve_right(vector(QQ, pv)))
    covered = []
    for fi in (21, 24):
        F = sorted(facs[fi][2])
        rng = bounded_range(R0, F)
        if rng is None:
            continue
        lo, hi = rng
        pts = sample(lo, hi)
        shape = None; dims = set()
        for tv in pts:
            R = R0 + tv * vector(QQ, dv)
            Q = Polyhedron(vertices=[list(vector(QQ, V[i])
                                          / R.dot_product(vector(QQ, V[i]))) for i in F],
                           base_ring=QQ)
            vv, ee, ff = cell_combinatorics(Q)
            sig = (Q.n_vertices(), len(ee), tuple(sorted(len(c) for c in ff)))
            shape = shape or sig
            if sig != shape: shape = None; break
            dims.add(summand_dim_from(vv, ee, ff))
        if shape is None or dims != {1}:
            print(f"    facet {fi} on ({lo}, {hi}): type not constant or not rigid, skipped")
            continue
        Ft = FractionField(PolynomialRing(QQ, 't')); tt = Ft.gen()
        t0 = pts[len(pts) // 2]
        R = R0 + t0 * vector(QQ, dv)
        vv, ee, ff = cell_combinatorics(Polyhedron(
            vertices=[list(vector(QQ, V[i]) / R.dot_product(vector(QQ, V[i]))) for i in F],
            base_ring=QQ))
        look = {tuple(vector(QQ, V[i]) / R.dot_product(vector(QQ, V[i]))): i for i in F}
        sym = []
        for v in vv:
            i = look[tuple(v)]
            den = (Ft(R0.dot_product(vector(QQ, V[i])))
                   + Ft(vector(QQ, dv).dot_product(vector(QQ, V[i]))) * tt)
            sym.append(vector(Ft, [Ft(c) / den for c in V[i]]))
        d = summand_dim_from(sym, ee, ff, field=Ft)
        ok(f"pairings {pv}, facet {fi}: cell bounded on ({lo}, {hi}), type {shape} "
           f"constant over {len(pts)} values, summand dim over Q(t) is {d}", d == 1)
        if d == 1:
            covered.append((lo, hi))
    ok(f"pairings {pv}: the rigid ranges cover the whole line -- "
       f"{covered[0]} and {covered[1]}",
       len(covered) == 2 and max(c[1] for c in covered) == Infinity
       and min(c[0] for c in covered) == -Infinity
       and max(c[0] for c in covered) < min(c[1] for c in covered))

print(f"\n{CH[0]} checks passed.")
print("""
Gate 1 now has no range left anywhere.  Twenty-two of the thirty lines die at
the crosscut, by a test that is constant along the line.  The other eight -- the
four pairing vectors above and their negatives, which give the same cells --
die because a pentagon-bearing facet cell is Minkowski-rigid over the whole
function field Q(t), on two overlapping ranges that between them cover every t.""")
