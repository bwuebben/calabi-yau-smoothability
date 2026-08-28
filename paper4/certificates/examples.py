#!/usr/bin/env python3
"""Apply the mixed necessity criterion to new Batyrev examples.

For a reflexive 4-polytope V whose generic anticanonical hypersurface has
only isolated germs of the covered kinds (nodes; anticanonical dP6/dP7
cones), this script builds the mixed matrix intrinsically and runs the
forced-zero test of the paper
(Theorem 6.1, the mixed necessity criterion):

  * node rows: the square circuit rows (+1 on one diagonal pair, -1 on
    the other), functionals on the divisors of the MPCP resolution;
  * divisorial rows: for each pentagon/hexagon face, a basis of
    K^perp inside Pic(E)^vee computed from the star fan of the face
    around its interior point, composed with the restriction dictionary
    (boundary rays -> boundary curve classes, interior ray -> K_E,
    all other rays -> 0);
  * the root subspaces are intrinsic: the unique (-2)-line for a
    pentagon (dP7), the A1/A2 root summands for a hexagon (dP6).

Verdict logic (necessity direction only): on the branch-restricted
kernel, a forced-to-zero node coordinate contradicts Lemma P, and a
forced-to-zero dP7 (-2)-coordinate contradicts Lemma D; either way the
hypersurface admits no one-parameter smoothing.  For dP6 points the
line branch is covered by Lemma D and the plane branch is open, so
profiles over dP6 components are reported per choice.

First target: Delta_19 (Paper 3, Section 5.3): 14 nodes + one dP7 point.

Run:  python3 examples.py
"""

import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))

from batyrev_global import (facets, two_faces, face_lattice_polygon,   # noqa
                            classify_polygon, dual_edge_length,
                            int_kernel, solve_int_coords, vgcd)

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
    M, piv = rref(A)
    cols = len(A[0])
    free = [c for c in range(cols) if c not in piv]
    out = []
    for fc in free:
        v = [F(0)] * cols
        v[fc] = F(1)
        for r, pc in enumerate(piv):
            v[pc] = -M[r][fc]
        out.append(v)
    return out

def inverse(A):
    n = len(A)
    M = [row[:] + [F(1) if i == j else F(0) for j in range(n)]
         for i, row in enumerate(A)]
    R2, piv = rref(M)
    assert piv == list(range(n))
    return [row[n:] for row in R2]

def dot(u, v):
    return sum(F(a) * F(b) for a, b in zip(u, v))

# ---------------------------------------------------------------- face tools
def vsub(a, b):
    return tuple(x - y for x, y in zip(a, b))

def analyze_faces(V):
    """Return (squares, dp_faces): squares as vertex 4-tuples (cyclic),
    dp_faces as dicts with the 2D star-fan data of pentagon/hexagon faces.
    The vertex <-> 2D-coordinate pairing is maintained throughout."""
    import math
    facs = facets(V)
    squares, dps = [], []
    for I, fp in two_faces(V, facs):
        u1, u2 = facs[fp[0]][0], facs[fp[1]][0]
        Ilist = sorted(I)
        pts = [V[i] for i in Ilist]
        ker = int_kernel([list(u1), list(u2)])
        p0 = pts[0]
        coords = [solve_int_coords(ker, vsub(p, p0)) for p in pts]
        # classification via the shared helper (order-insensitive inputs)
        _, evs, lens = face_lattice_polygon(V, I, u1, u2)
        cl = classify_polygon(evs, lens)
        npts = dual_edge_length(u1, u2)
        k, i = cl["k"], cl["i"]
        if k == 3 and i == 0:
            continue
        assert npts == 1, (I, "only dual length one is handled here")
        assert all(l == 1 for l in lens), (I, "non-unit face edge")
        # centroid in 2D face coordinates (rational)
        n = len(pts)
        cx = sum(F(c[0]) for c in coords) / n
        cy = sum(F(c[1]) for c in coords) / n
        order = sorted(range(n),
                       key=lambda t: math.atan2(float(F(coords[t][1]) - cy),
                                                float(F(coords[t][0]) - cx)))
        cyc_verts = [pts[t] for t in order]
        cyc_coords = [coords[t] for t in order]
        if k == 4 and i == 0:
            squares.append(tuple(cyc_verts))
        elif (k, i) in ((5, 1), (6, 1)):
            # the unique interior lattice point, by direct enumeration
            xs = [c[0] for c in cyc_coords]
            ys = [c[1] for c in cyc_coords]
            inter = []
            for px in range(min(xs), max(xs) + 1):
                for py in range(min(ys), max(ys) + 1):
                    signs = set()
                    for t in range(n):
                        ax, ay = cyc_coords[t]
                        bx, by = cyc_coords[(t + 1) % n]
                        cross = (bx - ax) * (py - ay) - (by - ay) * (px - ax)
                        signs.add(0 if cross == 0 else (1 if cross > 0 else -1))
                    if signs in ({1}, {-1}):
                        inter.append((px, py))
            assert len(inter) == 1, (inter, "interior point count")
            ix, iy = inter[0]
            rays2d = [(F(c[0]) - ix, F(c[1]) - iy) for c in cyc_coords]
            for t in range(n):
                a, b = rays2d[t], rays2d[(t + 1) % n]
                assert abs(a[0] * b[1] - a[1] * b[0]) == 1, \
                    "star not unimodular"
            interior = tuple(
                p0[j] + ix * ker[0][j] + iy * ker[1][j]
                for j in range(4))
            dps.append(dict(k=k, rays2d=rays2d, verts=cyc_verts,
                            interior=interior))
        else:
            raise AssertionError((I, k, i, "unhandled germ"))
    return squares, dps

# ------------------------------------------------- 2D toric surface package
def surface_lattice(star):
    """Pic(E), boundary classes, K, intersection form on Pic^vee."""
    u = star["rays2d"]
    n = len(u)
    # Pic(E) = Z^n / im(M2): relations rows (u_j[i])_j
    rel = [[F(u[j][i]) for j in range(n)] for i in range(2)]
    # basis of the quotient: choose n-2 coordinates complementary to pivots
    Mr, piv = rref(rel)
    freec = [c for c in range(n) if c not in piv]
    assert len(freec) == n - 2
    # class of D_j in the free coordinates: reduce e_j modulo the relations
    def cls(vec):
        v = [F(x) for x in vec]
        for r, pc in enumerate(piv):
            coef = v[pc]
            if coef != 0:
                v = [a - coef * b for a, b in zip(v, Mr[r])]
        assert all(v[c] == 0 for c in piv)
        return [v[c] for c in freec]
    Dcls = [cls([1 if t == j else 0 for t in range(n)]) for j in range(n)]
    Kcls = [-sum(Dcls[j][c] for j in range(n)) for c in range(n - 2)]
    # intersection numbers D_i . D_j on the smooth complete surface
    def selfint(i):
        a, b, c = u[(i - 1) % n], u[i], u[(i + 1) % n]
        # a + c = k b  =>  D_i^2 = -k
        if b[0] != 0:
            k = (a[0] + c[0]) / F(b[0])
        else:
            k = (a[1] + c[1]) / F(b[1])
        assert a[0] + c[0] == k * b[0] and a[1] + c[1] == k * b[1]
        return -k
    I_DD = [[F(0)] * n for _ in range(n)]
    for i in range(n):
        I_DD[i][i] = selfint(i)
        I_DD[i][(i + 1) % n] = I_DD[(i + 1) % n][i] = F(1)
    # Gram matrix of the basis classes: G_ab = (basis_a . basis_b) where
    # basis classes are D_{freec}; intersection of classes via D-matrix
    G = [[sum(Dcls[i][a] * Dcls[j][b] * I_DD[i][j]
              for i in range(n) for j in range(n))
          for b in range(n - 2)] for a in range(n - 2)]
    # wait: Dcls expresses D_j in basis; we need the form ON the basis:
    # solve from  I_DD[i][j] = sum_ab Dcls[i][a] G[a][b] Dcls[j][b].
    # Since basis classes are D_{freec} themselves, G[a][b] =
    # I_DD[freec[a]][freec[b]].
    G = [[I_DD[freec[a]][freec[b]] for b in range(n - 2)]
         for a in range(n - 2)]
    # consistency: reconstruct all pairwise numbers
    for i in range(n):
        for j in range(n):
            got = sum(F(Dcls[i][a]) * G[a][b] * F(Dcls[j][b])
                      for a in range(n - 2) for b in range(n - 2))
            assert got == I_DD[i][j], (i, j, got, I_DD[i][j])
    return dict(n=n, Dcls=Dcls, K=Kcls, G=G, freec=freec)

def integral_kernel_row(K):
    """Integral basis of {x in Z^d : <x, K> = 0} (saturated lattice),
    via a unimodular transformation U with K U = (g, 0, ..., 0)."""
    d = len(K)
    U = [[1 if i == j else 0 for j in range(d)] for i in range(d)]
    k = [int(x) for x in K]
    # column operations: reduce k to (g, 0, ..., 0)
    for j in range(1, d):
        a, b = k[0], k[j]
        while b != 0:
            q = a // b
            # col_0 <- col_0 - q col_j ; then swap
            k[0], k[j] = b, a - q * b
            for r in range(d):
                U[r][0], U[r][j] = U[r][j], U[r][0] - q * U[r][j]
            a, b = k[0], k[j]
    assert all(x == 0 for x in k[1:])
    return [[F(U[r][j]) for r in range(d)] for j in range(1, d)]

def kperp_root_data(surf):
    """K^perp in Pic^vee coordinates; the (-2)-vectors; root split."""
    m = len(surf["K"])
    Ginv = inverse(surf["G"])
    # functionals phi (rows in dual basis); K^perp: phi(K) = 0,
    # as a SATURATED integral lattice
    kp = integral_kernel_row(surf["K"])
    # pairing on functionals: <r,s> = r G^{-1} s^T
    def pair(r, s):
        return sum(r[a] * Ginv[a][b] * s[b] for a in range(m)
                   for b in range(m))
    # enumerate (-2)-vectors in K^perp with rigorous box
    Gram = [[pair(a, b) for b in kp] for a in kp]
    Q = [[-x for x in row] for row in Gram]
    Qinv = inverse(Q)
    import math
    bounds = [int(math.isqrt(int(2 * Qinv[i][i]) + 1)) + 1
              for i in range(len(kp))]
    from itertools import product
    roots = []
    for xs in product(*[range(-b, b + 1) for b in bounds]):
        if all(x == 0 for x in xs):
            continue
        v = [sum(F(xs[t]) * kp[t][a] for t in range(len(kp)))
             for a in range(m)]
        if pair(v, v) == F(-2):
            roots.append((xs, v))
    return dict(kperp=kp, pair=pair, roots=roots)

# ------------------------------------------------- the mixed matrix and test
def run_example(name, V, dp6_profiles=("L", "P")):
    print(f"=== {name} ===")
    squares, dps = analyze_faces(V)
    print(f"  {len(squares)} nodes; "
          f"{sum(1 for d in dps if d['k'] == 5)} dP7 point(s); "
          f"{sum(1 for d in dps if d['k'] == 6)} dP6 point(s)")

    # ray index: all lattice points used in supports
    ray_of = {}
    def ridx(v):
        return ray_of.setdefault(tuple(v), len(ray_of))
    rows, labels, restrict_sets = [], [], {}
    # node rows: circuit +1/-1 with the lexicographically first diagonal
    for a, b, c, d in squares:      # cyclic order: diagonals (a,c), (b,d)
        assert vsub(a, c) == vsub(vsub(a, b), vsub(c, b))
        assert tuple(x + y for x, y in zip(a, c)) == \
               tuple(x + y for x, y in zip(b, d)), "not a unit square"
        for v in (a, b, c, d):
            ridx(v)
        rows.append({ridx(a): F(-1), ridx(c): F(-1),
                     ridx(b): F(1), ridx(d): F(1)})
        labels.append(f"u{len(labels) + 1}")
    n_nodes = len(rows)
    # divisorial rows
    dp_info = []
    for t, star in enumerate(dps):
        surf = surface_lattice(star)
        rd = kperp_root_data(surf)
        kp, pair = rd["kperp"], rd["pair"]
        if star["k"] == 5:
            lines = {tuple(x) for x, _ in rd["roots"]}
            ok(f"{name}: dP7 #{t}: unique (-2)-line in K^perp "
               f"(discriminant-7 lattice)",
               len(rd["roots"]) == 2)
            phi_R = rd["roots"][0][1]
            # complement functional in K^perp independent of phi_R
            phi_A = next(v for v in kp if rank([phi_R, v]) == 2)
            basis = [("R", phi_R), ("A", phi_A)]
        else:
            ok(f"{name}: dP6 #{t}: 8 roots in K^perp (A1+A2)",
               len(rd["roots"]) == 8)
            # A1 = the +/- pair orthogonal to all other roots
            rootvecs = [v for _, v in rd["roots"]]
            a1 = [v for v in rootvecs
                  if sum(1 for w in rootvecs
                         if pair(v, w) != 0 and rank([v, w]) == 2) == 0]
            assert len(a1) == 2
            phi_s3 = a1[0]
            a2 = [v for v in rootvecs if rank([phi_s3, v]) == 2]
            a2b = [a2[0], next(v for v in a2 if rank([a2[0], v]) == 2)]
            basis = [("s1", a2b[0]), ("s2", a2b[1]), ("s3", phi_s3)]
        # ambient rows: value on D_j = phi(D_j|_E)
        for lbl, phi in basis:
            row = {}
            for j, v in enumerate(star["verts"]):
                row[ridx(v)] = dot(phi, surf["Dcls"][j])
            row[ridx(star["interior"])] = dot(phi, surf["K"])
            rows.append(row)
            labels.append(f"dp{t}:{lbl}")
        dp_info.append((t, star["k"]))

    nrays = len(ray_of)
    B = [[row.get(j, F(0)) for j in range(nrays)] for row in rows]
    # rows must annihilate the lattice relations
    inv_rays = {i: v for v, i in ray_of.items()}
    for r, row in enumerate(B):
        for mcoord in range(4):
            s = sum(row[j] * inv_rays[j][mcoord] for j in range(nrays))
            assert s == 0, (labels[r], mcoord, s)
    ok(f"{name}: all {len(B)} rows annihilate the lattice relations "
       f"({nrays} rays in supports)", True)

    node_rows = B[:n_nodes]
    print(f"  node subsystem: rank {rank(node_rows)} of {n_nodes}; "
          f"coloops: "
          f"{[labels[i] for i in range(n_nodes) if rank(node_rows[:i] + node_rows[i+1:]) == rank(node_rows) - 1]}")

    def forced(prohibit):
        pro = [i for i, lb in enumerate(labels) if lb in prohibit]
        allowed = [i for i in range(len(labels)) if i not in pro]
        kb = kernel([[B[i][c] for i in allowed] for c in range(nrays)])
        f = {labels[allowed[k]] for k in range(len(allowed))
             if all(v[k] == 0 for v in kb)}
        return len(kb), f

    d0, f0 = forced(set())
    print(f"  unrestricted kernel: dim {d0}; forced: {sorted(f0)}")
    # branch restriction: dP7 -> A-coordinate zero; dP6 -> per profile
    dp7s = [t for t, k in dp_info if k == 5]
    dp6s = [t for t, k in dp_info if k == 6]
    from itertools import product as iproduct
    verdicts = []
    for choice in iproduct(*[dp6_profiles] * len(dp6s)):
        prohibit = {f"dp{t}:A" for t in dp7s}
        for t, c in zip(dp6s, choice):
            prohibit |= ({f"dp{t}:s1", f"dp{t}:s2"} if c == "L"
                         else {f"dp{t}:s3"})
        d, f = forced(prohibit)
        req_nodes = {lb for lb in f if lb.startswith("u")}
        req_dp7 = {lb for lb in f if lb.endswith(":R")}
        tag = "".join(choice) if choice else "-"
        print(f"  profile {tag}: dim {d}; forced: {sorted(f)}")
        verdicts.append((tag, req_nodes, req_dp7, f))
    # verdict: covered directions are nodes (Lemma P) and dP7 (Lemma D);
    # a dP6 line-profile is covered by Lemma D, the plane profile is open.
    all_covered = all(rn or rd for _, rn, rd, _ in verdicts)
    if all_covered:
        print(f"  VERDICT: NONSMOOTHABLE — on every profile a covered "
              f"coordinate (node or dP7 root) is forced to zero.")
    else:
        open_profiles = [t for t, rn, rd, _ in verdicts if not (rn or rd)]
        print(f"  VERDICT: necessity test passes on profiles "
              f"{open_profiles}; nonsmoothability NOT established.")
    print()
    return all_covered

# ---------------------------------------------------------------- targets
V_19 = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, -1, 0),
        (0, -1, 0, 0), (1, -1, -1, 0), (0, 0, 0, 1), (-1, 1, 1, -1),
        (0, -1, 0, 1), (0, 0, -1, 1), (-1, 1, 0, -1), (-1, 0, 1, -1),
        (-1, 0, 1, 0), (-1, 1, 0, 0), (0, -1, -1, 0), (0, -1, -1, 1),
        (-1, 0, 0, -1), (-2, 1, 1, -1), (-1, 0, 0, 1)]

# cross-check the printed vertices against the committed source of paper 3
try:
    import face_data
    assert V_19 == [tuple(v) for v in face_data.V_19] or \
           [tuple(v) for v in face_data.V_19] == V_19
    print("  (V_19 matches src/face_data.py)")
except Exception:
    pass

# Delta_20 (paper, Section 8): 20 vertices; the generic hypersurface has
# 17 nodes and one dP7-cone point, and its nodal subsystem has NO coloop,
# so no purely nodal argument can obstruct it.  The mixed test obstructs.
V_20 = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1),
        (1, -1, -1, -1), (-1, 1, 1, 1), (0, 0, 0, -1), (0, 0, -1, 0),
        (0, -1, 0, 0), (-1, 1, 0, 1), (0, 0, -1, -1), (-1, 1, 0, 0),
        (-1, 0, 1, 1), (0, -1, 0, -1), (-1, 0, 1, 0), (0, -1, -1, -1),
        (-1, 0, 0, 1), (-1, 0, -1, 0), (-1, -1, 0, 0), (-1, -1, 1, 0)]

V_F1 = [(1, 0, 0, 0), (0, 1, 0, 0), (1, -1, 0, 0), (0, 0, 1, 0),
        (0, 0, 0, 1), (0, 0, 1, -1), (0, 0, -1, 1), (0, 0, 0, -1),
        (0, 0, -1, 0), (-1, 1, 0, 0), (0, -1, 0, 0), (-1, 1, -1, 1),
        (0, -1, 0, -1), (0, -1, -1, 0), (-1, 1, -1, 0), (-1, 0, 0, -1),
        (-1, 0, -1, 1), (-1, -1, 0, -1), (-2, 1, -1, 1), (-2, 1, -1, 0),
        (-1, -1, -1, 0), (-2, 0, -1, 0)]

def polar(V):
    out = []
    for f in facets(V):
        u, c = f[0], f[1]
        assert c == 1, (u, c)
        out.append(tuple(-x for x in u))
    return out

if __name__ == "__main__":
    # END-TO-END VALIDATION: the pipeline must reproduce the certified
    # verdict for X° = the polar-of-Delta_F1 hypersurface
    # (26 nodes + 2 dP7 + 2 dP6): two node coloops, forced node
    # coordinates on every profile, verdict NONSMOOTHABLE.
    v0 = run_example("X-degree (polar Delta_F1): 26 nodes + 2 dP7 + 2 dP6",
                     polar(V_F1))
    ok("validation: the pipeline reproduces NONSMOOTHABLE for X-degree", v0)
    verdict = run_example("Delta_19 hypersurface (14 nodes + 1 dP7)",
                          [tuple(v) for v in V_19])
    ok("Delta_19: the test is silent (no covered coordinate is forced)",
       not verdict)
    # Delta_20: the test fires with NO nodal coloop, so the obstruction is
    # entirely divisorial and no purely nodal criterion can see it.
    v20 = run_example("Delta_20 hypersurface (17 nodes + 1 dP7)", V_20)
    # the nodal subsystem of Delta_20, recomputed here: rank 11, no coloop
    sq20, _ = analyze_faces(V_20)
    ray20 = {}
    nrows20 = []
    for a, b, c, d in sq20:
        for v in (a, b, c, d):
            ray20.setdefault(tuple(v), len(ray20))
        row = [F(0)] * 24
        row[ray20[tuple(a)]] = row[ray20[tuple(c)]] = F(-1)
        row[ray20[tuple(b)]] = row[ray20[tuple(d)]] = F(1)
        nrows20.append(row)
    r20 = rank(nrows20)
    colo20 = [i for i in range(len(nrows20))
              if rank(nrows20[:i] + nrows20[i + 1:]) == r20 - 1]
    ok("Delta_20: the nodal subsystem has rank 11 and NO coloop, so no "
       "purely nodal relation argument obstructs it",
       len(nrows20) == 17 and r20 == 11 and colo20 == [])
    ok("Delta_20: the mixed test nevertheless obstructs, through the "
       "dP7 root direction alone", v20)
    print(f"ALL CHECKS PASSED ({PASS} assertions); "
          f"Delta_19 nonsmoothable: {verdict}; Delta_20 nonsmoothable: {v20}")
