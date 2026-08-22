#!/usr/bin/env python3
"""Exact certificate for the mirror partner X of Paper 4 (the section on
partial smoothings): the relation matrix of X with the rigid F_1 cone kept
singular, reconstructed from the printed data alone.

Input: the 22 vertices of Delta_{F_1} (Paper 3, Theorem 6.1), the twenty unit
squares with the diagonals of the fixed crepant fan
(mirror_partner_fan.sage), and the two exceptional surfaces given by their
cyclic star rays in the induced plane lattice of the face.  Everything else
-- Picard groups, canonical classes, intersection forms, K^perp lattices,
root decompositions, restriction dictionaries -- is derived here.

Checks: the star fans are unimodular and complete; K^perp(dP6) = A1 + A2 and
K^perp(F_1) is the line of a primitive (-8)-class with no roots; every row
annihilates the lattice relations; the nodal subsystem has rank 14 with
coloops exactly n_9, n_10; the 24-row matrix has rank 16 with 8-dimensional
relation space; the coordinates of n_9 and n_10 vanish identically on that
space and on both dP6 branch profiles, while the F_1 coordinate is a coloop.

Pure Python, exact rational arithmetic.  Run:  python3 mirror_partner.py
"""

import sys
from fractions import Fraction as F
from itertools import product

PASS = 0


def ok(label, cond):
    global PASS
    if not cond:
        print(f"FAIL: {label}")
        sys.exit(1)
    PASS += 1
    print(f"  pass: {label}")


# ---------------------------------------------------------------- linear algebra
def rref(A):
    M = [[F(x) for x in row] for row in A]
    rows, cols = len(M), len(M[0]) if M else 0
    piv, r = [], 0
    for c in range(cols):
        p = next((i for i in range(r, rows) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        M[r] = [x / M[r][c] for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c] != 0:
                M[i] = [a - M[i][c] * b for a, b in zip(M[i], M[r])]
        piv.append(c)
        r += 1
    return M, piv


def rank(A):
    return len(rref(A)[1]) if A else 0


def kernel(A):
    """Basis of {x : A x = 0}."""
    M, piv = rref(A)
    cols = len(A[0])
    basis = []
    for f in (c for c in range(cols) if c not in piv):
        v = [F(0)] * cols
        v[f] = F(1)
        for r, c in enumerate(piv):
            v[c] = -M[r][f]
        basis.append(v)
    return basis


def inverse(A):
    n = len(A)
    aug = [[F(x) for x in A[i]] + [F(1) if j == i else F(0) for j in range(n)]
           for i in range(n)]
    M, piv = rref(aug)
    assert piv == list(range(n)), "singular matrix"
    return [row[n:] for row in M]


def dot(a, b):
    return sum(F(x) * F(y) for x, y in zip(a, b))


# ---------------------------------------------------------------- printed data
# Paper 3, Theorem 6.1, in the paper's order; 23 and 24 are the interior
# points of the hexagon face and of the F_1 face.
V = {
    1: (1, 0, 0, 0), 2: (0, 1, 0, 0), 3: (1, -1, 0, 0), 4: (0, 0, 1, 0),
    5: (0, 0, 0, 1), 6: (0, 0, 1, -1), 7: (0, 0, -1, 1), 8: (0, 0, 0, -1),
    9: (0, 0, -1, 0), 10: (-1, 1, 0, 0), 11: (0, -1, 0, 0),
    12: (-1, 1, -1, 1), 13: (0, -1, 0, -1), 14: (0, -1, -1, 0),
    15: (-1, 1, -1, 0), 16: (-1, 0, 0, -1), 17: (-1, 0, -1, 1),
    18: (-1, -1, 0, -1), 19: (-2, 1, -1, 1), 20: (-2, 1, -1, 0),
    21: (-1, -1, -1, 0), 22: (-2, 0, -1, 0),
    23: (-1, 0, -1, 0), 24: (-1, 0, 0, 0),
}

# The twenty unit squares, each with the diagonal of the fixed crepant fan.
SQUARES = (
    ((1, 2, 7, 12), (2, 7)), ((1, 2, 9, 15), (2, 9)), ((1, 3, 8, 13), (3, 8)),
    ((1, 3, 9, 14), (3, 9)), ((2, 10, 12, 19), (10, 12)),
    ((2, 10, 15, 20), (10, 15)), ((3, 11, 13, 18), (11, 13)),
    ((3, 11, 14, 21), (11, 14)), ((5, 11, 17, 21), (11, 17)),
    ((6, 10, 16, 20), (10, 16)), ((7, 9, 12, 15), (9, 12)),
    ((7, 12, 17, 19), (12, 17)), ((7, 14, 17, 21), (14, 17)),
    ((8, 9, 13, 14), (9, 13)), ((8, 13, 16, 18), (13, 16)),
    ((8, 15, 16, 20), (15, 16)), ((12, 15, 19, 20), (15, 19)),
    ((13, 14, 18, 21), (14, 18)), ((16, 18, 20, 22), (18, 20)),
    ((17, 19, 21, 22), (19, 21)),
)

# The two divisorial faces: cyclic boundary rays, the interior point, and the
# star rays in the induced plane lattice of the face, based at that point.
SURFACES = {
    "H6": dict(
        interior=23,
        rays=(21, 14, 9, 15, 20, 22),
        rays2d=((0, -1), (1, -1), (1, 0), (0, 1), (-1, 1), (-1, 0)),
    ),
    "F1": dict(
        interior=24,
        rays=(18, 4, 19, 22),
        rays2d=((0, -1), (1, 0), (-1, 1), (-1, 0)),
    ),
}

print("== the two exceptional surfaces ==")
for name, S in SURFACES.items():
    u = S["rays2d"]
    n = len(u)
    ok(f"{name}: the star fan is unimodular and complete "
       f"(consecutive rays span the lattice)",
       all(abs(u[i][0] * u[(i + 1) % n][1] - u[i][1] * u[(i + 1) % n][0]) == 1
           for i in range(n)))
    # The printed 2-D rays must be the face vertices measured from the
    # interior point, read in a basis of the plane lattice of the face.
    # Solve for that basis from the first two rays (their determinant is 1,
    # so the solution is integral) and verify it on all the others.
    w = [tuple(V[r][m] - V[S["interior"]][m] for m in range(4))
         for r in S["rays"]]
    det = u[0][0] * u[1][1] - u[0][1] * u[1][0]
    assert abs(det) == 1
    basis = [[F(u[1][1] * w[0][m] - u[0][1] * w[1][m], det) for m in range(4)],
             [F(-u[1][0] * w[0][m] + u[0][0] * w[1][m], det) for m in range(4)]]
    ok(f"{name}: the printed 2-D rays are the face vertices seen from the "
       f"interior point, in an integral basis of the plane lattice of the "
       f"face",
       all(x.denominator == 1 for b in basis for x in b)
       and all(tuple(u[j][0] * basis[0][m] + u[j][1] * basis[1][m]
                     for m in range(4)) == tuple(F(x) for x in w[j])
               for j in range(n)))
    # Pic(E) = Z^n / im(M -> Z^n): relations are the two coordinate rows
    rel = [[F(u[j][i]) for j in range(n)] for i in range(2)]
    M, piv = rref(rel)
    free = [c for c in range(n) if c not in piv]
    assert len(free) == n - 2

    def cls(vec, M=M, piv=piv, free=free):
        v = [F(x) for x in vec]
        for r, pc in enumerate(piv):
            if v[pc] != 0:
                coef = v[pc]
                v = [a - coef * b for a, b in zip(v, M[r])]
        return [v[c] for c in free]

    Dcls = [cls([1 if t == j else 0 for t in range(n)]) for j in range(n)]
    K = [-sum(Dcls[j][c] for j in range(n)) for c in range(n - 2)]
    # toric intersection numbers on the smooth complete surface
    I = [[F(0)] * n for _ in range(n)]
    for i in range(n):
        a, b, c = u[(i - 1) % n], u[i], u[(i + 1) % n]
        k = (F(a[0] + c[0]) / b[0]) if b[0] != 0 else (F(a[1] + c[1]) / b[1])
        assert a[0] + c[0] == k * b[0] and a[1] + c[1] == k * b[1]
        I[i][i] = -k
        I[i][(i + 1) % n] = I[(i + 1) % n][i] = F(1)
    G = [[I[free[a]][free[b]] for b in range(n - 2)] for a in range(n - 2)]
    ok(f"{name}: the Gram matrix of the Picard basis reproduces every "
       f"D_i . D_j", all(
        sum(Dcls[i][a] * G[a][b] * Dcls[j][b]
            for a in range(n - 2) for b in range(n - 2)) == I[i][j]
        for i in range(n) for j in range(n)))
    S.update(Dcls=Dcls, K=K, G=G, Ginv=inverse(G))
    # K^perp and its roots
    kp = kernel([[F(x) for x in K]])          # rational basis; saturate below
    den = [max(x.denominator for x in v) for v in kp]
    kp = [[x * d for x in v] for v, d in zip(kp, den)]

    def pair(a, b, S=S):
        m = len(S["K"])
        return sum(F(a[i]) * S["Ginv"][i][j] * F(b[j])
                   for i in range(m) for j in range(m))

    roots = []
    B = 4
    for xs in product(range(-B, B + 1), repeat=len(kp)):
        if not any(xs):
            continue
        v = [sum(F(xs[t]) * kp[t][a] for t in range(len(kp)))
             for a in range(len(K))]
        if pair(v, v) == F(-2):
            roots.append(v)
    S.update(kperp=kp, pair=pair, roots=roots)
    ok(f"{name}: K^perp has rank {n - 3}", len(kp) == n - 3)

# the A1 + A2 decomposition for the hexagon
H = SURFACES["H6"]
ok("H6: K^perp contains exactly 8 roots (the root system A1 + A2)",
   len(H["roots"]) == 8)
a1 = [v for v in H["roots"]
      if all(H["pair"](v, w) == 0 or rank([v, w]) == 1 for w in H["roots"])]
ok("H6: exactly one +/- pair of roots is orthogonal to all others (the A1 "
   "summand)", len(a1) == 2)
s3 = a1[0]
a2 = [v for v in H["roots"] if rank([s3, v]) == 2]
s1 = a2[0]
s2 = next(v for v in a2 if rank([s1, v]) == 2)
ok("H6: s1, s2 span the A2 summand (s1^2 = s2^2 = -2, s1.s2 = +/-1) and s3 "
   "spans A1 (s3^2 = -2)",
   H["pair"](s1, s1) == -2 and H["pair"](s2, s2) == -2
   and abs(H["pair"](s1, s2)) == 1 and H["pair"](s3, s3) == -2
   and H["pair"](s3, s1) == 0 and H["pair"](s3, s2) == 0)
H["rows"] = {"s1": s1, "s2": s2, "s3": s3}

Fq = SURFACES["F1"]
ok("F1: K^perp is a line and contains NO roots: the F_1 cone has no "
   "(-2)-class, unlike every del Pezzo cone of Paper 4",
   len(Fq["kperp"]) == 1 and Fq["roots"] == [])
t = Fq["kperp"][0]
ok("F1: its primitive generator has self-intersection -8",
   Fq["pair"](t, t) == -8
   and __import__("math").gcd(*[abs(int(x)) for x in t]) == 1)
Fq["rows"] = {"t": t}

# ---------------------------------------------------------------- the matrix
# Reduced-rigidity of the F_1 germ, over the LATTICE.  By Altmann the
# reduced miniversal base is a point exactly when the edge multiset has no
# proper nonempty zero-sum sub-multiset (no nontrivial Minkowski
# decomposition into lattice polygons).  The cone of RATIONAL summands is
# two-dimensional here, which is why dim T^1 = 1; the two are different
# statements and only the lattice one gives rigidity.
from itertools import combinations as _comb  # noqa: E402


def _edges(rays2d):
    m = len(rays2d)
    return [(rays2d[(i + 1) % m][0] - rays2d[i][0],
             rays2d[(i + 1) % m][1] - rays2d[i][1]) for i in range(m)]


_ev = _edges(SURFACES["F1"]["rays2d"])
ok("F1: the four edge vectors of the quadrilateral sum to zero",
   sum(e[0] for e in _ev) == 0 and sum(e[1] for e in _ev) == 0)
_proper = [t for r in (1, 2, 3) for t in _comb(range(4), r)
           if sum(_ev[i][0] for i in t) == 0
           and sum(_ev[i][1] for i in t) == 0]
ok("F1: NO proper nonempty subset of the edge vectors sums to zero, so the "
   "quadrilateral has no nontrivial Minkowski decomposition into lattice "
   "polygons and the germ is reduced-rigid", _proper == [])
_hv = _edges(SURFACES["H6"]["rays2d"])
_hp = [t for r in (1, 2, 3, 4, 5) for t in _comb(range(6), r)
       if sum(_hv[i][0] for i in t) == 0 and sum(_hv[i][1] for i in t) == 0]
ok("H6: the hexagon by contrast HAS proper zero-sum subsets, so it is "
   "Minkowski decomposable over the lattice and its germ is smoothable",
   len(_hp) > 0)

print("\n== the 24-row relation matrix of X ==")
rows, labels = [], []
for k, (sq, (a, c)) in enumerate(SQUARES, start=1):
    b, d = sorted(set(sq) - {a, c})
    assert tuple(x + y for x, y in zip(V[a], V[c])) == \
        tuple(x + y for x, y in zip(V[b], V[d])), (sq, "not a unit square")
    row = [F(0)] * 24
    row[a - 1] = row[c - 1] = F(-1)
    row[b - 1] = row[d - 1] = F(1)
    rows.append(row)
    labels.append(f"n{k}")
ok("all twenty squares satisfy their diagonal relation v_a + v_c = v_b + v_d "
   "(so the circuit row is independent of which diagonal the fan chose, up "
   "to sign)", len(rows) == 20)
for name, S in SURFACES.items():
    for lb, phi in S["rows"].items():
        row = [F(0)] * 24
        for j, r in enumerate(S["rays"]):
            row[r - 1] = dot(phi, S["Dcls"][j])
        row[S["interior"] - 1] = dot(phi, S["K"])
        rows.append(row)
        labels.append(f"{name}:{lb}")
ok("24 rows on the 24 boundary lattice points of Delta", len(rows) == 24)
ok("every row annihilates the lattice relations: sum_j row_j v_j = 0",
   all(sum(rows[r][j] * V[j + 1][m] for j in range(24)) == 0
       for r in range(24) for m in range(4)))
ok("every row has coefficient sum zero (it is a functional on Pic)",
   all(sum(r) == 0 for r in rows))

node_rows = rows[:20]
ok("the nodal subsystem has rank 14", rank(node_rows) == 14)
colo = [labels[i] for i in range(20)
        if rank(node_rows[:i] + node_rows[i + 1:]) == 13]
ok("its coloops are exactly n_9 = conv{v5,v11,v17,v21} and "
   "n_10 = conv{v6,v10,v16,v20}", colo == ["n9", "n10"])
ok("the full 24-row matrix has rank 16", rank(rows) == 16)


def forced(prohibit):
    allowed = [i for i, lb in enumerate(labels) if lb not in prohibit]
    kb = kernel([[rows[i][c] for i in allowed] for c in range(24)])
    f = {labels[allowed[k]] for k in range(len(allowed))
         if all(v[k] == 0 for v in kb)}
    return len(kb), f


d0, f0 = forced(set())
ok("the unrestricted relation space has dimension 8", d0 == 8)
ok("on it, the coordinates of n_9, n_10, F1:t and H6:s3 vanish identically",
   f0 == {"n9", "n10", "F1:t", "H6:s3"})
dL, fL = forced({"H6:s1", "H6:s2"})
ok("dP6 line profile (A1 branch): dimension 6; n_9, n_10 still forced",
   dL == 6 and {"n9", "n10"} <= fL)
dP, fP = forced({"H6:s3"})
ok("dP6 plane profile (A2 branch): dimension 8; n_9, n_10 still forced",
   dP == 8 and {"n9", "n10"} <= fP)
# The F_1 point enters B_Delta only through one coordinate and a ZERO
# column, so the forcing is computed from the nodes and the hexagon alone.
# Deleting the row by itself would prove nothing: its coordinate already
# vanishes on the whole relation space, so deleting it intersects with a
# hyperplane that contains that space.
iF1 = SURFACES["F1"]["interior"] - 1
iH6 = SURFACES["H6"]["interior"] - 1
ok("the columns of B_Delta at the two interior points v_23, v_24 are "
   "identically zero (node rows miss them; divisorial rows have entry "
   "<phi, K_E> = 0 because phi lies in K^perp)",
   all(r[iF1] == 0 for r in rows) and all(r[iH6] == 0 for r in rows))
_keep_r = [i for i, lb in enumerate(labels) if lb != "F1:t"]
_keep_c = [c for c in range(24) if c != iF1]
_kb = kernel([[rows[i][c] for i in _keep_r] for c in _keep_c])
_f = {labels[_keep_r[k]] for k in range(len(_keep_r))
      if all(v[k] == 0 for v in _kb)}
ok("discarding the F_1 row AND its column altogether still forces n_9 and "
   "n_10: the obstruction does not use the F_1 cone", {"n9", "n10"} <= _f)
# Strongest form: replace the F_1 block by every functional supported on the
# rays of the F_1 face, so that no identification of H_2 of its link is
# assumed at all.
_rows2 = [r for i, r in enumerate(rows) if labels[i] != "F1:t"]
_lab2 = [lb for lb in labels if lb != "F1:t"]
for _j, _r in enumerate(SURFACES["F1"]["rays"]):
    _row = [F(0)] * 24
    _row[_r - 1] = F(1)
    _rows2.append(_row)
    _lab2.append(f"F1free{_j}")
_kb2 = kernel([[_rows2[i][c] for i in range(len(_rows2))] for c in range(24)])
_f2 = {_lab2[k] for k in range(len(_lab2)) if all(v[k] == 0 for v in _kb2)}
ok("n_9 and n_10 remain forced even when the F_1 block is replaced by the "
   "full space of functionals on its face rays, so Theorem 7.2 does not "
   "depend on identifying H_2 of the F_1 link", {"n9", "n10"} <= _f2)
ok("the F_1 coordinate is itself a coloop: no relation class is nonzero on "
   "the F_1 link", "F1:t" in f0)

print(f"\nALL CHECKS PASSED ({PASS} assertions)")
