#!/usr/bin/env python3
"""The examples X_9, X_19, X_20, X_88 and X_154 of Section 11 and Appendix F, checked from the data files alone.

Run from the root of the repository (standard library only, Python >= 3.10):

    python3 paper6/examples/examples_check.py                 # all five examples
    python3 paper6/examples/examples_check.py --only X9 X19   # some of them

Each data file paper6/examples/x*.json describes one example (the fields are documented in paper6/README.md).
For each example the script recomputes, in exact integer or rational arithmetic:

 1. The facets and two-faces of Delta from its vertices; that Delta is admissible (Definition 2.1: primitive
    edges; every two-face a unimodular triangle, a unit square, a dP_7 pentagon or a dP_6 hexagon; dual edge of
    lattice length one at every non-triangular two-face); the numbers of squares, pentagons and hexagons.
 2. The node parallelograms q = (a, b, c, d): a + c = b + d; each square is one node; each pentagon carries the two
    parallelograms of its unique subdivision; each hexagon carries the parallelograms of the profile (A_1 or A_2
    subdivision); the chosen diagonal {a, c} passes through the centre at a pentagon or hexagon.
 3. The circuits circ_q = e_b + e_d - e_a - e_c, a Z-basis of K = {lambda : sum lambda_q circ_q = 0} (Hermite
    normal form), rank K, and that the printed basis spans K; the forced nodes (coordinates vanishing on K); the
    printed relations with all coefficients nonzero lie in K; |Q| - rank K.
 4. The unimodular triangulation T of the boundary of Delta° and its height h: every recorded point is a nonzero
    lattice point of Delta° (and all of them are recorded), every tetrahedron lies in a facet and has determinant
    +-1, every triangle lies in two tetrahedra on opposite sides, each of three fixed rays (three of the five rays of
    condition (iv) of Appendix A) lies on no wall, that is, in the cone over no triangle of a tetrahedron, and lies in
    the cone over exactly one tetrahedron, and the bends of h and of h - phi across every wall are positive
    (Theorem F; Appendix A, conditions (i)-(v), with (iv) required of all three rays).
 5. The fan Sigma' of (Proj) and its height: the rays are lattice points of the boundary of Delta; the cells are the
    cells of the subdivision induced by the height (strict convexity, with a linear function on each cell below the
    height at every other ray); the cells of each facet of Delta fill it (volumes); on every non-triangular two-face
    they induce the profile.
 6. The cycles: the discriminant Gamma of the base from the stored polyhedral decomposition (a leg for each edge e
    with two-point lower support in a two-cell f between maximal cells with different facet labels); the boundary
    legs of each stored cycle lie on Gamma; the boundary graph (closed paths, or trivalent vertices at negative
    vertices of Gamma) and its segments; the canonical lift of the boundary network and the node coefficients
    rho_q (Definition 3.7: the lift after p_q minus the lift before is rho_q circ_q), with balance at every other
    point; positive vertices met; the node coefficients generate K (Smith normal form); for belt cycles, a disc
    with one field up to sign. The oriented triangle-chain boundary is checked against the recorded
    network in the tangent charts, and triangle fields are checked primitive and tangent.
 7. For X_9 in addition: the six vectors u_0, ..., u_5 of Section 11.1 in T; the parallelogram sides parallel to
    d_0; the disc S of 26 triangles, its three maximal cells with facet label u_0, its 18 legs, the order in which
    its boundary meets the two-faces F_{u_0 u_i} and the nodes, and rho(S) = q_1 - q_2 - q_3 + q_4.
 8. For X_88 and X_154: the tie at the hexagon, the projection of K to the hexagon nodes, and the other tilings of
    the hexagon.

Euler characteristics, four-valent vertices of the interior branching graph, and interior intersection
points with Gamma are also counted directly from the triangles. The link at each counted intersection
is checked to be a circle. This counts the asserted crossing points; it does not by itself prove affine
transversality or every clause of [CBM, Def. 7.2].
What is not recomputed here: all the clauses of [CBM, Def. 7.2] for the cycles whose domain has interior vertices or
crosses Gamma, and the construction of the base from T, h and Sigma'. The
stored decomposition is used as given; step 6 checks the statements about node coefficients, which depend only on
the boundary legs, their lower supports and the facet labels.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction
from math import gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent
T0 = time.time()
FAILURES = []


def check(cond, msg):
    if not cond:
        FAILURES.append(msg)
        print(f"  FAIL  {msg}", flush=True)
    return cond


def ok(msg):
    print(f"  ok    {msg}", flush=True)


# ----------------------------------------------------------------------------------------------- linear algebra
def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def det(M):
    """Exact determinant (Fraction elimination; entries may be int or Fraction)."""
    A = [[Fraction(x) for x in r] for r in M]
    n = len(A)
    d = Fraction(1)
    for k in range(n):
        p = next((i for i in range(k, n) if A[i][k] != 0), None)
        if p is None:
            return Fraction(0)
        if p != k:
            A[k], A[p] = A[p], A[k]
            d = -d
        d *= A[k][k]
        for i in range(k + 1, n):
            f = A[i][k] / A[k][k]
            if f:
                for j in range(k, n):
                    A[i][j] -= f * A[k][j]
    return d


def solve(M, b):
    """x with x . M = b (M square, rows = basis), exact."""
    n = len(M)
    A = [[Fraction(M[j][i]) for j in range(n)] + [Fraction(b[i])] for i in range(n)]
    for k in range(n):
        p = next(i for i in range(k, n) if A[i][k] != 0)
        A[k], A[p] = A[p], A[k]
        for i in range(n):
            if i != k and A[i][k]:
                f = A[i][k] / A[k][k]
                A[i] = [x - f * y for x, y in zip(A[i], A[k])]
    return [A[i][n] / A[i][i] for i in range(n)]


def rank(rows):
    A = [[Fraction(x) for x in r] for r in rows]
    r = 0
    ncol = len(A[0]) if A else 0
    for c in range(ncol):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        for i in range(len(A)):
            if i != r and A[i][c]:
                f = A[i][c] / A[r][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        r += 1
    return r


def left_kernel(M):
    """Z-basis of {x in Z^m : x M = 0} for an integer m x n matrix M (rows of a unimodular transform whose
    image rows vanish; the basis is saturated)."""
    m = len(M)
    n = len(M[0]) if m else 0
    A = [list(M[i]) + [1 if j == i else 0 for j in range(m)] for i in range(m)]
    r = 0
    for c in range(n):
        while True:
            nz = [i for i in range(r, m) if A[i][c] != 0]
            if not nz:
                break
            p = min(nz, key=lambda i: abs(A[i][c]))
            A[r], A[p] = A[p], A[r]
            done = True
            for i in range(r + 1, m):
                if A[i][c]:
                    q = A[i][c] // A[r][c]
                    A[i] = [x - q * y for x, y in zip(A[i], A[r])]
                    if A[i][c]:
                        done = False
            if done:
                break
        if any(A[i][c] for i in range(r, m)):
            r += 1
    return [row[n:] for row in A[r:]]


def elementary_divisors(M):
    """Nonzero elementary divisors of an integer matrix (Smith normal form)."""
    A = [list(r) for r in M]
    out = []
    while A and A[0]:
        entries = [(abs(A[i][j]), i, j) for i in range(len(A)) for j in range(len(A[0])) if A[i][j]]
        if not entries:
            break
        _, i, j = min(entries)
        A[0], A[i] = A[i], A[0]
        for row in A:
            row[0], row[j] = row[j], row[0]
        while True:
            changed = False
            for i in range(1, len(A)):
                if A[i][0]:
                    q = A[i][0] // A[0][0]
                    A[i] = [x - q * y for x, y in zip(A[i], A[0])]
                    if A[i][0]:
                        A[0], A[i] = A[i], A[0]
                        changed = True
            for j in range(1, len(A[0])):
                if A[0][j]:
                    q = A[0][j] // A[0][0]
                    for row in A:
                        row[j] -= q * row[0]
                    if A[0][j]:
                        for row in A:
                            row[0], row[j] = row[j], row[0]
                        changed = True
            if not changed:
                bad = [(i, j) for i in range(1, len(A)) for j in range(1, len(A[0])) if A[i][j] % A[0][0]]
                if not bad:
                    break
                i, _ = bad[0]
                A[0] = [x + y for x, y in zip(A[0], A[i])]
        out.append(abs(A[0][0]))
        A = [row[1:] for row in A[1:]]
    return out


def same_lattice(B1, B2):
    """The rows of B1 and of B2 span the same sublattice of Z^n."""
    if rank(B1) != rank(B2) or rank(B1 + B2) != rank(B1):
        return False
    prod = lambda ed: eval("*".join(map(str, ed)) or "1")  # noqa: E731
    p1, p2, p12 = prod(elementary_divisors(B1)), prod(elementary_divisors(B2)), prod(elementary_divisors(B1 + B2))
    return p1 == p12 and p2 == p12


def in_lattice(v, B):
    return same_lattice(B, B + [list(v)]) if B else not any(v)


def vgcd(v):
    g = 0
    for x in v:
        g = gcd(g, x)
    return g


# ----------------------------------------------------------------------------------------------- polytopes
def cross4(a, b, c):
    """A normal vector of the span of three vectors of R^4 (3x3 minors)."""
    M = [a, b, c]
    out = []
    for k in range(4):
        cols = [j for j in range(4) if j != k]
        out.append((-1) ** k * det([[M[i][j] for j in cols] for i in range(3)]))
    return tuple(int(x) for x in out)


def facets_of(V):
    """Facet normals u (integer) of the reflexive polytope conv(V): <u, x> >= -1 on V, = -1 on the facet."""
    out = {}
    for a, b, c, d in itertools.combinations(range(len(V)), 4):
        n = cross4(sub(V[b], V[a]), sub(V[c], V[a]), sub(V[d], V[a]))
        if not any(n):
            continue
        vals = [dot(n, v) for v in V]
        lo, hi = min(vals), max(vals)
        base = dot(n, V[a])
        if base == lo:
            u = n
        elif base == hi:
            u = tuple(-x for x in n)
            base = -base
        else:
            continue
        assert base != 0 and (-1) % 1 == 0
        if base % 1 or (-1 * 1) % 1:
            pass
        g = vgcd(u)
        u = tuple(x // g for x in u)
        base = dot(u, V[a])
        assert base < 0 and base == -1, "polytope not reflexive"
        out[u] = frozenset(i for i, v in enumerate(V) if dot(u, v) == -1)
    return out


def affine_rank(P):
    if len(P) <= 1:
        return len(P) - 1
    return rank([sub(p, P[0]) for p in P[1:]])


def two_faces_of(V, F):
    faces = set()
    items = list(F.items())
    for (u1, s1), (u2, s2) in itertools.combinations(items, 2):
        s = s1 & s2
        if len(s) >= 3 and affine_rank([V[i] for i in s]) == 2:
            faces.add(s)
    faces = {s for s in faces if not any(s < t for t in faces)}
    out = []
    for s in faces:
        us = [u for u, t in F.items() if s <= t]
        assert len(us) == 2
        out.append((sorted(s), us))
    return out


def tri_area2(p, q, r):
    """Twice the lattice area of the triangle pqr (gcd of the 2x2 minors of q-p, r-p)."""
    a, b = sub(q, p), sub(r, p)
    return vgcd([a[i] * b[j] - a[j] * b[i] for i, j in itertools.combinations(range(4), 2)])


def cyclic_order(P):
    """The points P (vertices of a convex lattice polygon in R^4) in cyclic order."""
    c = [Fraction(sum(p[i] for p in P), len(P)) for i in range(4)]
    for i, j in itertools.combinations(range(4), 2):
        pr = [(p[i] - c[i], p[j] - c[j]) for p in P]
        if affine_rank([(x, y, 0, 0) for x, y in pr]) == 2:
            break

    def half(v):
        return 0 if (v[1] > 0 or (v[1] == 0 and v[0] > 0)) else 1

    import functools

    def cmp(a, b):
        ha, hb = half(a[1]), half(b[1])
        if ha != hb:
            return ha - hb
        cr = a[1][0] * b[1][1] - a[1][1] * b[1][0]
        return -1 if cr > 0 else (1 if cr < 0 else 0)
    return [p for p, _ in sorted(zip(P, pr), key=functools.cmp_to_key(cmp))]


def polygon_data(P):
    cyc = cyclic_order(P)
    k = len(cyc)
    edges_primitive = all(vgcd(sub(cyc[(i + 1) % k], cyc[i])) == 1 for i in range(k))
    a2 = sum(tri_area2(cyc[0], cyc[i], cyc[i + 1]) for i in range(1, k - 1))
    return cyc, edges_primitive, a2


# ----------------------------------------------------------------------------------------------- the checks
def check_faces(ex, V, rays):
    F = facets_of(V)
    faces = two_faces_of(V, F)
    kinds = Counter()
    nontri = []
    adm = True
    for idx, us in faces:
        P = [V[i] for i in idx]
        cyc, prim, a2 = polygon_data(P)
        adm &= prim
        kind = {(3, 1): "triangle", (4, 2): "square", (5, 5): "pentagon", (6, 6): "hexagon"}.get((len(P), a2))
        adm &= kind is not None
        if kind and kind != "triangle":
            adm &= vgcd(sub(us[0], us[1])) == 1
            centre = None
            if kind in ("pentagon", "hexagon"):
                inside = [j for j, r in enumerate(rays) if j >= len(V) and in_polygon_interior(r, cyc)]
                check(len(inside) == 1, f"{ex}: one ray in the interior of the {kind}")
                centre = inside[0]
            nontri.append(dict(kind=kind, cyclic=[rays.index(p) for p in cyc], centre=centre))
        kinds[kind] += 1
    check(adm, f"{ex}: Delta is admissible")
    ok(f"{ex}: {len(F)} facets, two-faces {dict(kinds)}; admissible")
    return F, nontri, kinds


def in_polygon_interior(r, cyc):
    """r lies in the relative interior of the convex polygon with vertices cyc (in cyclic order)."""
    if affine_rank(cyc + [r]) != 2:
        return False
    k = len(cyc)
    c = [Fraction(sum(p[i] for p in cyc), k) for i in range(4)]
    signs = set()
    for i in range(k):
        p, q = cyc[i], cyc[(i + 1) % k]
        # orientation of (p, q, r) against (p, q, centroid) in the plane: compare 2x2 minors
        s = orient2(p, q, r, c)
        if s == 0:
            return False
        signs.add(s)
    return len(signs) == 1


def orient2(p, q, r, c):
    """+1 if r and c lie on the same side of the line pq in the plane of the polygon, -1 otherwise, 0 on it."""
    a, b, d = sub(q, p), [x - y for x, y in zip(r, p)], [x - y for x, y in zip(c, p)]
    for i, j in itertools.combinations(range(4), 2):
        m1 = a[i] * b[j] - a[j] * b[i]
        m2 = a[i] * d[j] - a[j] * d[i]
        if m2 != 0:
            return 0 if m1 == 0 else (1 if (m1 > 0) == (m2 > 0) else -1)
    raise AssertionError


def tiles_of(face, nodes_in_face):
    """The tiles (as sets of ray indices) of the subdivision of a pentagon or hexagon: node parallelograms and
    the triangles {o, w_j, w_{j+1}} not covered by them."""
    w, o = face["cyclic"], face["centre"]
    k = len(w)
    tiles = [frozenset(q) for q in nodes_in_face]
    covered = set()
    for q in nodes_in_face:
        s = [x for x in q if x != o]
        for j in range(k):
            if {w[j], w[(j + 1) % k]} <= set(s):
                covered.add(j)
    for j in range(k):
        if j not in covered:
            tiles.append(frozenset({o, w[j], w[(j + 1) % k]}))
    return tiles


def check_nodes(ex, D, rays, nontri):
    nodes = [tuple(q) for q in D["nodes"]]
    R = sorted({x for q in nodes for x in q})
    ok_par = all(add(rays[a], rays[c]) == add(rays[b], rays[d]) for a, b, c, d in nodes)
    check(ok_par, f"{ex}: every node parallelogram satisfies n_a + n_c = n_b + n_d")
    squares = [f for f in nontri if f["kind"] == "square"]
    stars = [f for f in nontri if f["kind"] != "square"]
    sq_sets = {frozenset(f["cyclic"]) for f in squares}
    node_sets = [frozenset(q) for q in nodes]
    check(sq_sets <= set(node_sets), f"{ex}: every unit square is a node parallelogram")
    face_nodes = [i + 1 for i, q in enumerate(node_sets) if q not in sq_sets]
    check(face_nodes == D["nodes_in_pentagon_or_hexagon"], f"{ex}: the nodes in pentagons or hexagons are "
          f"{['q%d' % i for i in face_nodes]}")
    branch = None
    for f in stars:
        inq = [q for q in nodes if set(q) <= set(f["cyclic"]) | {f["centre"]} and f["centre"] in q]
        check(all(q[0] == f["centre"] or q[2] == f["centre"] for q in inq),
              f"{ex}: the chosen diagonal passes through the centre of the {f['kind']}")
        w, o, k = f["cyclic"], f["centre"], len(f["cyclic"])
        allowed = [frozenset({o, w[j], w[(j + 1) % k], w[(j + 2) % k]}) for j in range(k)
                   if add(rays[w[j]], rays[w[(j + 2) % k]]) == add(rays[o], rays[w[(j + 1) % k]])]
        check(all(frozenset(q) in allowed for q in inq), f"{ex}: the {f['kind']} nodes are parallelograms "
              "{o, w_j, w_j+1, w_j+2} of the face")
        starts = sorted(j for j in range(k) for q in inq
                        if frozenset(q) == frozenset({o, w[j], w[(j + 1) % k], w[(j + 2) % k]}))
        if k == 5:
            check(len(inq) == 2 and (starts[1] - starts[0]) % 5 in (2, 3), f"{ex}: the pentagon carries the two "
                  "parallelograms of its subdivision")
        else:
            if len(inq) == 3 and all((starts[i] - starts[0]) % 2 == 0 for i in range(3)):
                branch = "A2"
            elif len(inq) == 2 and (starts[1] - starts[0]) % 6 == 3:
                branch = "A1"
            check(branch is not None, f"{ex}: the hexagon carries an A_1 or A_2 subdivision")
            ok(f"{ex}: the hexagon profile is on the {branch} branch")
        f["nodes"] = [node_sets.index(frozenset(q)) for q in inq]
    check(len(nodes) == len(squares) + sum(len(f["nodes"]) for f in stars),
          f"{ex}: every node lies in a square, pentagon or hexagon")
    ok(f"{ex}: |Q| = {len(nodes)} ({len(squares)} squares, "
       f"{sum(len(f['nodes']) for f in stars)} in pentagons and hexagons); R has {len(R)} points")
    return nodes, branch


def circuits(nodes, nrays):
    C = []
    for a, b, c, d in nodes:
        v = [0] * nrays
        v[b] += 1
        v[d] += 1
        v[a] -= 1
        v[c] -= 1
        C.append(v)
    return C


def forced_of(Kb, n):
    return [q for q in range(n) if all(r[q] == 0 for r in Kb)]


def check_K(ex, D, nodes, nrays):
    C = circuits(nodes, nrays)
    Kb = left_kernel(C)
    n = len(nodes)
    check(all(not any(dot(r, [C[q][j] for q in range(n)]) for j in range(nrays)) for r in Kb),
          f"{ex}: the computed basis consists of relations")
    rk = len(Kb)
    check(same_lattice(Kb, D["K_basis"]), f"{ex}: the printed Z-basis of K spans K")
    check(len(D["K_basis"]) == rk, f"{ex}: rank K = {rk}")
    forced = forced_of(Kb, n)
    lam = D.get("relation_with_all_coefficients_nonzero")
    if lam:
        check(in_lattice(lam, Kb) and all(lam), f"{ex}: the printed relation with all coefficients nonzero lies "
              "in K")
    ok(f"{ex}: rank K = {rk}, |Q| - rank K = {n - rk}, forced nodes "
       f"{['q%d' % (q + 1) for q in forced] or 'none'}")
    return C, Kb, forced


def check_triangulation(ex, V, Fdelta):
    """T and h: conditions (i)-(v) of Appendix A for the stored tetrahedra."""
    T = D_cur["triangulation"]
    pts = [tuple(p) for p in T["points"]]
    h = T["height"]
    nz = [i for i, p in enumerate(pts) if any(p)]
    # (i) the recorded nonzero points are exactly the nonzero lattice points of Delta°
    U = list(Fdelta)                     # vertices of Delta° = facet normals of Delta
    box = [range(min(u[i] for u in U), max(u[i] for u in U) + 1) for i in range(4)]
    lat = {m for m in itertools.product(*box) if any(m) and all(dot(m, v) >= -1 for v in V)}
    check(lat == {pts[i] for i in nz}, f"{ex}: the points of T are the {len(lat)} nonzero lattice points of Delta°")
    check(all(min(dot(pts[i], v) for v in V) == -1 for i in nz), f"{ex}: they lie on the boundary")
    check(all(isinstance(x, int) for x in h), f"{ex}: integral heights")
    tets = [tuple(t) for t in T["tetrahedra"]]
    # (ii)
    for t in tets:
        P = [pts[i] for i in t]
        if not check(any(all(dot(p, v) == -1 for p in P) for v in V), f"{ex}: tetrahedron {t} in a facet"):
            return
        if not check(abs(det(P)) == 1, f"{ex}: tetrahedron {t} unimodular"):
            return
    # (iii)
    walls = defaultdict(list)
    for t in tets:
        for tau in itertools.combinations(sorted(t), 3):
            walls[tau].append(next(i for i in t if i not in tau))
    good = all(len(a) == 2 for a in walls.values())
    good &= all(det([pts[i] for i in tau] + [pts[a[0]]]) * det([pts[i] for i in tau] + [pts[a[1]]]) < 0
                for tau, a in walls.items())
    check(good, f"{ex}: every triangle lies in two tetrahedra, on opposite sides")
    # (iv)
    rays = [(1000003, 999983, -1000037, 1000039), (-7919, 104729, 1299709, -15485863),
            (104723, -104717, 104711, 104707)]
    for r in rays:
        coords = [solve([pts[i] for i in t], r) for t in tets]
        check(not any(min(c) == 0 for c in coords), f"{ex}: a fixed ray lies on no wall of T")
        hits = sum(1 for c in coords if min(c) > 0)
        check(hits == 1, f"{ex}: a fixed ray lies in the cone over exactly one tetrahedron of T")
    # (v)
    minb = None
    for tau, (a, b) in walls.items():
        c = solve([pts[i] for i in tau] + [pts[a]], pts[b])
        bend = h[b] - sum(ci * h[i] for ci, i in zip(c, list(tau) + [a]))
        e = 1 - sum(c)                    # the bend of phi, which is 1 on the boundary
        check(bend > 0 and bend - e > 0, f"{ex}: h and h - phi bend positively across the wall {tau}")
        minb = bend if minb is None else min(minb, bend)
    ok(f"{ex}: T has {len(tets)} unimodular tetrahedra on {len(nz)} points; h and h - phi strictly convex "
       f"(max |h| = {max(abs(x) for x in h)})")
    return pts, tets


def volume3(P):
    """Normalised volume of conv(0, P) for points P spanning a 3-dimensional polytope in an affine hyperplane
    not through 0 (= |det| summed over a triangulation)."""
    P = [tuple(p) for p in set(map(tuple, P))]
    g = [Fraction(sum(p[i] for p in P), len(P)) for i in range(4)]
    faces = set()
    for a, b, c in itertools.combinations(P, 3):
        s = [det([a, b, c, x]) for x in P]
        if any(s) and (all(x >= 0 for x in s) or all(x <= 0 for x in s)):
            faces.add(frozenset(x for x, y in zip(P, s) if y == 0))
    vol = Fraction(0)
    for Fc in faces:
        Fl = list(Fc)
        # polygon edges: pairs with every other point of the face on one side (within the face's plane)
        verts = [x for x in Fl if not any(between(x, y, z) for y, z in itertools.combinations(Fl, 2)
                                          if x not in (y, z))]
        for x, y in itertools.combinations(verts, 2):
            s = [det([g, x, y, z]) for z in verts if z not in (x, y)]
            if all(v >= 0 for v in s) or all(v <= 0 for v in s):
                if x != verts[0] and y != verts[0]:
                    vol += abs(det([g, verts[0], x, y]))
    return vol


def between(x, y, z):
    """x lies strictly between y and z on a line."""
    d1, d2 = sub(x, y), sub(z, y)
    if rank([d1, d2]) != 1:
        return False
    k = next(i for i in range(4) if d2[i])
    t = Fraction(d1[k], d2[k])
    return 0 < t < 1


def check_fan(ex, V, Fdelta, rays, nontri):
    fan = D_cur["fan"]
    cells = [tuple(c) for c in fan["cells"]]
    h = fan["height"]
    on_bd = all(min(dot(u, r) for u in Fdelta) == -1 for r in rays)
    check(on_bd, f"{ex}: every ray of Sigma' is a lattice point of the boundary of Delta")
    extra = [j for j in range(len(V), len(rays)) if all(j != f["centre"] for f in nontri)]
    for j in extra:
        on = [u for u in Fdelta if dot(u, rays[j]) == -1]
        check(len(on) == 1, f"{ex}: the further ray n_{j} lies in the interior of a facet of Delta")
    # each cell lies in a facet; a linear function equal to h on the cell and below h at every other ray
    for c in cells:
        fac = [u for u in Fdelta if all(dot(u, rays[i]) == -1 for i in c)]
        if not check(len(fac) >= 1 and rank([rays[i] for i in c]) == 4, f"{ex}: cell {c} is 3-dimensional "
                     "and lies in a facet"):
            return
        base = []
        for i in c:
            if rank([rays[j] for j in base] + [rays[i]]) > len(base):
                base.append(i)
        coef = solve([[rays[j][k] for k in range(4)] for j in base], [0, 0, 0, 0]) if False else None
        # the linear function l with l(n_j) = h_j on the basis: solve B^T l = h
        B = [rays[j] for j in base]
        lvec = solve([[B[r][k] for r in range(4)] for k in range(4)], [h[j] for j in base])
        lval = lambda r: sum(a * b for a, b in zip(lvec, r))  # noqa: E731
        check(all(lval(rays[i]) == h[i] for i in c), f"{ex}: the height is linear on the cell {c}")
        check(all(h[i] > lval(rays[i]) for i in range(len(rays)) if i not in c),
              f"{ex}: the height is strictly convex at the cell {c}")
    # the cells in each facet fill it
    for u, vs in Fdelta.items():
        incell = [c for c in cells if all(dot(u, rays[i]) == -1 for i in c)]
        vf = volume3([V[i] for i in vs])
        vc = sum(volume3([rays[i] for i in c]) for c in incell)
        check(vf == vc, f"{ex}: the cells in the facet {u} fill it (volume {vf})")
    # induced subdivision of the non-triangular two-faces
    for f in nontri:
        pts = set(f["cyclic"]) | ({f["centre"]} if f["centre"] is not None else set())
        got = set()
        for c in cells:
            s = frozenset(set(c) & pts)
            if len(s) >= 3 and affine_rank([rays[i] for i in s]) == 2:
                got.add(s)
        got = {s for s in got if not any(s < t for t in got)}
        if f["kind"] == "square":
            want = {frozenset(f["cyclic"])}
        else:
            want = set(tiles_of(f, [D_cur["nodes"][q] for q in f["nodes"]]))
        check(got == want, f"{ex}: Sigma' induces the profile on the {f['kind']} {f['cyclic']}")
    ok(f"{ex}: Sigma' has {len(cells)} cells on {len(rays)} rays, induced by the stored height, filling every "
       "facet and inducing the profile")


# ----------------------------------------------------------------------------------------------- the base
class Base:
    def __init__(self, cx):
        self.vx = [tuple(v[:4]) for v in cx["vertices"]]
        self.low = [v[4] for v in cx["vertices"]]
        self.E = [frozenset(e) for e in cx["edges"]]
        self.F = [frozenset(f) for f in cx["two_cells"]]
        self.C = [frozenset(c) for c in cx["three_cells"]]
        self.lab = [tuple(u) for u in cx["facet_labels"]]
        v2C = defaultdict(set)
        for i, c in enumerate(self.C):
            for v in c:
                v2C[v].add(i)
        self.v2C = v2C
        self.cof = [sorted(set.intersection(*(v2C[v] for v in f))) for f in self.F]
        assert all(len(x) == 2 for x in self.cof), "a two-cell lies in two maximal cells"
        v2E = defaultdict(list)
        for i, e in enumerate(self.E):
            for v in e:
                v2E[v].append(i)
        self.legs = {}
        for fi, f in enumerate(self.F):
            ca, cb = self.cof[fi]
            if self.lab[ca] == self.lab[cb]:
                continue
            for ei in {ei for v in f for ei in v2E[v] if self.E[ei] <= f}:
                lo = sorted({self.low[v] for v in self.E[ei]})
                if len(lo) == 2:
                    self.legs[(ei, fi)] = tuple(lo)
        self.legs_at = defaultdict(list)
        for (ei, fi) in self.legs:
            self.legs_at[("e", ei)].append((ei, fi))
            self.legs_at[("f", fi)].append((ei, fi))

    def lower(self, cell):
        return frozenset(self.low[v] for v in cell)

    def labels(self, cell):
        return frozenset(self.lab[i] for i in set.intersection(*(self.v2C[v] for v in cell)))

    def kind(self, p, node_by_support):
        k, i = p
        n = len(self.legs_at[p])
        if k == "e":
            return "positive" if n == 3 else "edge point"
        lo = self.lower(self.F[i])
        if lo in node_by_support:
            return "node"
        if len(lo) == 3 and n == 3:
            return "negative"
        return "two-cell point"


def check_cycles(ex, D, nodes, C, Kb, rays):
    B = Base(D["complex"])
    node_by_support = {frozenset(q): k for k, q in enumerate(nodes)}
    nodepts = [("f", fi) for fi, f in enumerate(B.F) if B.lower(f) in node_by_support
               and B.lab[B.cof[fi][0]] != B.lab[B.cof[fi][1]]]
    check(len(nodepts) == len(nodes), f"{ex}: Gamma has one node point on the two-cell over each node "
          "parallelogram")
    four = [p for p in nodepts if len(B.legs_at[p]) == 4]
    check(len(four) == len(nodes), f"{ex}: each node point is a four-valent vertex of Gamma")
    ok(f"{ex}: the decomposition has f-vector {[len(B.vx), len(B.E), len(B.F), len(B.C)]}; Gamma has "
       f"{len(B.legs)} legs")
    rhos = []
    for cyc in D["cycles"]:
        rhos.append(check_one_cycle(ex, D, B, cyc, nodes, C, node_by_support, rays))
    if rhos:
        ed = elementary_divisors(rhos)
        rk = rank(rhos)
        check(rk == len(Kb) and all(x == 1 for x in ed), f"{ex}: the node coefficient vectors of the "
              f"{len(rhos)} cycles generate K (rank {rk}, elementary divisors all 1)")
        if "cycle_coordinates_printed" in D:
            coords = [[int(x) for x in solve_in_basis(r, D["K_basis"])] for r in rhos]
            check(coords == D["cycle_coordinates_printed"], f"{ex}: coordinates of the rho(S_i) in the printed "
                  f"basis: {coords}")
        ok(f"{ex}: the node coefficients of {len(rhos)} cycles generate K with index one")
    return rhos, B


def solve_in_basis(v, Bs):
    """Coordinates of v in the basis rows Bs (exact; the rows are independent)."""
    n = len(v)
    cols = []
    for j in range(n):
        if rank([[r[k] for k in cols + [j]] for r in Bs]) > len(cols):
            cols.append(j)
        if len(cols) == len(Bs):
            break
    x = solve([[r[k] for k in cols] for r in Bs], [v[k] for k in cols])
    assert all(sum(xi * r[k] for xi, r in zip(x, Bs)) == v[k] for k in range(n))
    return x


def check_one_cycle(ex, D, B, cyc, nodes, C, node_by_support, rays):
    name = cyc["name"]
    check_chain(ex, cyc, B, rays)
    legs = [(e, f, c) for e, f, c in cyc["boundary_legs"]]
    check(all((e, f) in B.legs for e, f, _ in legs), f"{ex} {name}: every boundary leg is a leg of Gamma")
    adj = defaultdict(list)
    for e, f, c in legs:
        adj[("e", e)].append((("f", f), (e, f), c))
        adj[("f", f)].append((("e", e), (e, f), c))
    kinds = {p: B.kind(p, node_by_support) for p in adj}
    # canonical lift and its boundary: sum over legs at a point of c (e_n' - e_n), with sign + at the two-cell end
    nr = len(rays)
    bd = defaultdict(lambda: [0] * nr)
    for e, f, c in legs:
        n, n2 = B.legs[(e, f)]
        for p, s in ((("f", f), 1), (("e", e), -1)):
            bd[p][n2] += s * c
            bd[p][n] -= s * c
    rho = [0] * len(nodes)
    bal = True
    for p, v in bd.items():
        if kinds[p] == "node":
            q = node_by_support[B.lower(B.F[p[1]])]
            m = [-x for x in v]
            if m == C[q]:
                rho[q] = 1
            elif m == [-x for x in C[q]]:
                rho[q] = -1
            elif any(m):
                check(False, f"{ex} {name}: the lift jumps by a multiple of circ at q{q + 1}")
        else:
            bal &= not any(v)
    check(bal, f"{ex} {name}: the canonical lift is balanced at every point of Gamma other than the nodes")
    check(not any(dot(rho, [C[q][j] for q in range(len(nodes))]) for j in range(nr)),
          f"{ex} {name}: its node coefficients lie in K")
    # boundary graph: strands between negative vertices, or closed paths
    negs = sorted(p for p in adj if kinds[p] == "negative" and len(adj[p]) == 3)
    val_ok = all(len(adj[p]) == (3 if p in negs else 2) for p in adj)
    check(val_ok, f"{ex} {name}: the boundary graph has valency 3 at its negative vertices and 2 elsewhere")
    strands = walk_strands(B, adj, kinds, negs, node_by_support, legs)
    npos = sum(s["positive"] for s in strands)
    if not negs:
        shape = "closed path" if len(strands) == 1 else f"{len(strands)} closed paths"
    elif len(negs) == 2 and len(strands) == 3:
        shape = "two vertices joined by three paths"
    elif len(negs) == 4 and len(strands) == 6 and len({tuple(sorted((s["start"], s["end"]))) for s in strands}) == 6:
        shape = "complete graph on four vertices"
    else:
        shape = f"{len(negs)} trivalent vertices, {len(strands)} paths"
    # belt cycles: one field up to sign, a disc
    fields = {tuple(t["field"]) for t in cyc["triangles"]}
    signs = {tuple(abs(x) for x in f) for f in fields}
    extra = ""
    cyc["_belt"] = False
    if cyc.get("belt_expected", False):
        check(not negs and len(strands) == 1, f"{ex} {name}: a belt boundary is a single closed path")
        v = next(iter(fields))
        const = all(f == v or f == tuple(-x for x in v) for f in fields)
        dirs = {tuple(sorted((sub(rays[B.legs[(e, f)][1]], rays[B.legs[(e, f)][0]]),
                              sub(rays[B.legs[(e, f)][0]], rays[B.legs[(e, f)][1]])))) for e, f, _ in legs}
        vv = tuple(sorted((v, tuple(-x for x in v))))
        check(const and dirs == {vv}, f"{ex} {name}: one field +-{v}, and every boundary leg has vanishing vector "
              "+-v")
        chi, nb = surface_invariants(cyc)
        check(chi == 1 and nb == 1, f"{ex} {name}: the triangles form a disc (Euler characteristic {chi}, "
              f"{nb} boundary circle)")
        boundary_only = gamma_on_boundary(cyc, B)
        check(boundary_only, f"{ex} {name}: the disc meets Gamma only on its boundary")
        cyc["_belt"] = const and dirs == {vv} and chi == 1 and nb == 1 and boundary_only
        extra = f", field +-{v}"
    pr = cyc.get("printed")
    if pr:
        check(rho == pr["node_coefficients"], f"{ex} {name}: node coefficients as printed")
        check(npos == pr["positive_vertices_on_boundary"], f"{ex} {name}: {npos} positive vertices on the boundary, "
              "as printed")
        if "field" in pr:
            v = tuple(pr["field"])
            check(signs == {tuple(abs(x) for x in v)}, f"{ex} {name}: the printed field")
        want = canonical_segments(pr["boundary_segments"], D)
        got = canonical_segments(strands_as_tokens(strands, negs, D, B), D)
        check(want == got, f"{ex} {name}: boundary segments as printed")
    profile = surface_profile(cyc, B)
    cyc['_surface_profile'] = profile
    for field in ('euler_characteristic', 'interior_four_valent', 'interior_gamma_points'):
        if pr and field in pr:
            check(profile[field] == pr[field], f"{ex} {name}: {field} agrees with the printed value")
    check(profile['gamma_links_circular'], f"{ex} {name}: interior Gamma intersections have circular links")
    nz = {f"q{q + 1}": r for q, r in enumerate(rho) if r}
    ok(f"{ex} {name}: {len(cyc['triangles'])} triangles (subdivision {cyc['subdivision']}), boundary graph: "
       f"{shape}, {len(legs)} legs, {npos} positive vertices{extra}; rho = {nz}")
    cyc["_rho"], cyc["_strands"], cyc["_shape"], cyc["_negs"] = rho, strands, shape, negs
    return rho


def check_chain(ex, cyc, B, rays):
    """Check the integral chain boundary in the base's tangent charts.

    Within a maximal cell the chart is its tangent hyperplane. On a lower
    face containing a vertex w it is N/Z n(w); equality there means that
    the difference is an integral multiple of the primitive ray n(w).
    This is a chain check, not all local embedding conditions for CBM cycles.
    """
    residual = defaultdict(lambda: [0]*4)
    cells = [{i} for i in range(len(B.vx))], B.E, B.F, B.C

    def add_edge(a, b, field, sign=1):
        if b < a:
            a, b, sign = b, a, -sign
        for j in range(4):residual[(a,b)][j] += sign*field[j]

    valid_flags = primitive = tangent = True
    for t in cyc['triangles']:
        flags = [tuple(tuple(x) for x in f) for f in t['flags']]
        support = sorted(set(x for f in flags for x in f))
        valid_flags &= all(cells[d][i] < cells[e][j] for (d,i),(e,j) in zip(support,support[1:]))
        if cyc['subdivision'] == 2:
            valid_flags &= all(set(a) < set(b) for a,b in zip(flags,flags[1:]))
        v = t['field']
        primitive &= gcd(*v) == 1
        dim, idx = support[-1]
        representatives = [idx] if dim == 3 else B.cof[idx]
        tangent &= any(dot(B.lab[c],v)==0 for c in representatives)
        for i in range(3):
            a,b = [f for j,f in enumerate(flags) if j != i]
            add_edge(a,b,v,(-1)**i)
    check(valid_flags, f"{ex} {cyc['name']}: triangles are flags of incident base cells")
    check(primitive and tangent, f"{ex} {cyc['name']}: primitive fields in their tangent hyperplanes")
    for e,f,c in cyc['boundary_legs']:
        n,m=B.legs[(e,f)]
        v=[c*(b-a) for a,b in zip(rays[n],rays[m])]
        a,b=((1,e),),((2,f),)
        if cyc['subdivision']==1:add_edge(a,b,v,-1)
        else:
            mid=((1,e),(2,f))
            add_edge(a,mid,v,-1);add_edge(mid,b,v,-1)
    failures=0
    for (a,b),v in residual.items():
        support=sorted(set(a+b))
        if support[-1][0]==3:
            good=not any(v)
        else:
            dim,idx=support[0]
            w=min(cells[dim][idx]);n=rays[B.low[w]]
            good=all(v[i]*n[j]==v[j]*n[i] for i in range(4) for j in range(i+1,4))
        failures += not good
    check(failures==0, f"{ex} {cyc['name']}: triangle-chain boundary equals the recorded discriminant network ({failures} inconsistent edges)")


def surface_profile(cyc, B):
    """Invariants of the stored simplicial domain, including its branching graph."""
    vertices, edges, links = set(), Counter(), defaultdict(list)
    triangles = []
    for t in cyc['triangles']:
        vs = [tuple(sorted(tuple(x) for x in flag)) for flag in t['flags']]
        triangles.append(frozenset(vs))
        vertices.update(vs)
        edges.update(frozenset(e) for e in itertools.combinations(vs, 2))
        for v in vs:
            links[v].append(tuple(w for w in vs if w != v))
    check(len(set(triangles)) == len(triangles), 'surface triangles are distinct')
    check(all(n in (1, 2, 3) for n in edges.values()), 'surface edges have one, two or three incident triangles')
    boundary = set().union(*(set(e) for e, n in edges.items() if n == 1))
    branching = Counter(v for e, n in edges.items() if n == 3 for v in e)
    def on_gamma(flag):
        if len(flag) == 1:
            dim, i = flag[0]
            return ('e' if dim == 1 else 'f' if dim == 2 else '', i) in B.legs_at
        return len(flag) == 2 and [x[0] for x in flag] == [1, 2] and (flag[0][1], flag[1][1]) in B.legs
    intersections = [v for v in vertices - boundary if on_gamma(v)]
    def circle(es):
        adj = defaultdict(set)
        for a, b in es:
            adj[a].add(b); adj[b].add(a)
        if not adj or any(len(x) != 2 for x in adj.values()):
            return False
        seen, todo = set(), [next(iter(adj))]
        while todo:
            x = todo.pop()
            if x not in seen:
                seen.add(x); todo.extend(adj[x] - seen)
        return len(seen) == len(adj)
    return dict(euler_characteristic=len(vertices)-len(edges)+len(triangles),
                interior_four_valent=sum(n == 4 and v not in boundary for v, n in branching.items()),
                interior_gamma_points=len(intersections),
                gamma_links_circular=all(circle(links[v]) for v in intersections))


def walk_strands(B, adj, kinds, negs, node_by_support, legs):
    used = set()
    strands = []

    def pair_of(frm, leg, c):
        n, n2 = B.legs[leg]
        s = c if frm[0] == "e" else -c
        return (n, n2) if s == 1 else (n2, n)

    def walk(start, first):
        p = start
        nxt, leg, c = first
        seq = []
        npos = 0
        while True:
            used.add(leg)
            pr = pair_of(p, leg, c)
            if not seq or seq[-1]["pair"] != pr:
                seq.append(dict(pair=pr, node_after=None, faces=[]))
            fl = B.labels(B.F[leg[1]])
            if not seq[-1]["faces"] or seq[-1]["faces"][-1] != fl:
                seq[-1]["faces"].append(fl)
            p = nxt
            if kinds[p] == "positive":
                npos += 1
            if kinds[p] == "node":
                seq[-1]["node_after"] = node_by_support[B.lower(B.F[p[1]])]
            if p in negs or p == start:
                return seq, p, npos
            cand = [x for x in adj[p] if x[1] != leg]
            nxt, leg, c = cand[0]

    for t in negs:
        for x in adj[t]:
            if x[1] in used:
                continue
            seq, end, npos = walk(t, x)
            strands.append(dict(start=t, end=end, seq=seq, positive=npos, closed=False))
    while True:
        rest = [(e, f) for e, f, _ in legs if (e, f) not in used]
        if not rest:
            break
        e, f = rest[0]
        start = ("f", f) if kinds[("f", f)] == "node" else ("e", e)
        # start at a node of the component if there is one
        comp, stack = set(), [start]
        while stack:
            p = stack.pop()
            if p in comp:
                continue
            comp.add(p)
            stack += [x[0] for x in adj[p]]
        nd = sorted((node_by_support[B.lower(B.F[p[1]])], p) for p in comp if kinds[p] == "node")
        if nd:
            start = nd[0][1]
        seq, end, npos = walk(start, adj[start][0])
        strands.append(dict(start=start, end=end, seq=seq, positive=npos, closed=True))
    return strands


def strands_as_tokens(strands, negs, D, B):
    names = {}
    pr = D.get("negative_vertices_printed", [])
    for t in negs:
        f = B.F[t[1]]
        key = (frozenset(B.lower(f)), frozenset(B.labels(f)))
        for rec in pr:
            if (frozenset(rec["lower_support"]), frozenset(tuple(u) for u in rec["face_label"])) == key:
                names[t] = rec["name"]
    out = []
    for s in strands:
        toks = []
        if not s["closed"]:
            toks.append(names.get(s["start"], "?"))
        for i, x in enumerate(s["seq"]):
            toks.append(list(x["pair"]))
            if x["node_after"] is not None:
                toks.append(f"q{x['node_after'] + 1}")
        if not s["closed"]:
            toks.append(names.get(s["end"], "?"))
        out.append(toks)
    return out


def canonical_segments(segs, D):
    """A normal form of a list of boundary segments, independent of the direction of listing and, for closed
    segments, of the starting point."""
    out = []
    for toks in segs:
        closed = not isinstance(toks[0], str) or not toks[0].startswith("t")
        forms = []
        for rev in (False, True):
            t = list(toks)
            if rev:
                t = [x[::-1] if isinstance(x, list) else x for x in t[::-1]]
            if closed:
                # rotate: pair, node, pair, node ...; after reversal the node precedes its pair
                if rev:
                    t = t[1:] + t[:1] if isinstance(t[0], str) else t
                units = [tuple(t[i:i + 2]) for i in range(0, len(t), 2)]
                units = [(tuple(a), b) for a, b in units]
                forms += [tuple(units[k:] + units[:k]) for k in range(len(units))]
            else:
                forms.append(tuple(tuple(x) if isinstance(x, list) else x for x in t))
        out.append(min(forms, key=repr))
    return sorted(out, key=repr)


def gamma_on_boundary(cyc, B):
    """Gamma is a subcomplex of the first barycentric subdivision.
    In subdivision two its additional vertices are midpoints of its legs.
    Check that no vertex of that subcomplex lies in the interior of the disc.
    """
    verts, edges = set(), Counter()
    for t in cyc["triangles"]:
        vs = [tuple(sorted(tuple(x) for x in flag)) for flag in t["flags"]]
        verts.update(vs)
        for x, y in itertools.combinations(vs, 2):
            edges[frozenset((x,y))] += 1
    boundary = set().union(*(set(e) for e,n in edges.items() if n == 1))
    def on_gamma(flag):
        if len(flag) == 1:
            dim,i = flag[0]
            return ("e" if dim == 1 else "f" if dim == 2 else "", i) in B.legs_at
        if len(flag) == 2 and [x[0] for x in flag] == [1,2]:
            return (flag[0][1],flag[1][1]) in B.legs
        return False
    return all(v in boundary for v in verts if on_gamma(v))


def surface_invariants(cyc):
    """Euler characteristic and number of boundary circles of the simplicial surface formed by the triangles."""
    V, E, F = set(), Counter(), 0
    for t in cyc["triangles"]:
        vs = [tuple(tuple(x) for x in fl) for fl in t["flags"]]
        V.update(vs)
        for a, b in itertools.combinations(vs, 2):
            E[frozenset((a, b))] += 1
        F += 1
    chi = len(V) - len(E) + F
    bd = [e for e, n in E.items() if n == 1]
    adj = defaultdict(list)
    for e in bd:
        a, b = tuple(e)
        adj[a].append(b)
        adj[b].append(a)
    if any(len(x) != 2 for x in adj.values()) or any(n > 2 for n in E.values()):
        return chi, -1
    seen, comps = set(), 0
    for v in adj:
        if v in seen:
            continue
        comps += 1
        stack = [v]
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            stack += adj[x]
    return chi, comps


# ----------------------------------------------------------------------------------------------- per example
def check_x9(D, rays, nodes, C, Kb, T, B):
    ex = "X9"
    d0 = (4, 2, 1, 1)
    n = rays
    check(sub(n[3], n[6]) == sub(n[8], n[10]) == sub(n[4], n[5]) == sub(n[7], n[2]) == d0,
          f"{ex}: n3 - n6 = n8 - n10 = n4 - n5 = n7 - n2 = d0 = {d0}")
    check([list(q) for q in nodes] == [[3, 6, 10, 8], [2, 6, 3, 7], [4, 5, 10, 8], [2, 5, 4, 7]],
          f"{ex}: q1 = (n3,n6,n10,n8), q2 = (n2,n6,n3,n7), q3 = (n4,n5,n10,n8), q4 = (n2,n5,n4,n7)")
    check(Kb == [[1, -1, -1, 1]] or Kb == [[-1, 1, 1, -1]], f"{ex}: K = Z (q1 - q2 - q3 + q4)")
    u0 = (0, 1, -1, -1)
    us = [(1, -1, -1, -1), (0, 0, 1, -1), (-1, 2, 1, -1), (-1, 2, 0, 0), (-1, 2, -1, 1)]
    check(all(dot(u, d0) == 0 for u in [u0] + us), f"{ex}: u0, ..., u5 annihilate d0")
    pts, tets = T
    idx = {p: i for i, p in enumerate(pts)}
    i0 = idx[u0]
    nbrs = {pts[j] for t in tets if i0 in t for j in t if j != i0}
    check(set(us) <= nbrs, f"{ex}: in T the vertex u0 is joined by edges to u1, ..., u5")
    tris = {frozenset(tau) for t in tets for tau in itertools.combinations(t, 3)}
    check(all(frozenset({i0, idx[us[i]], idx[us[(i + 1) % 5]]}) in tris for i in range(5)),
          f"{ex}: [u0, u_i, u_i+1] is a triangle of T for i = 1, ..., 5")
    # the disc
    S = D["cycles"][0]
    check(S["subdivision"] == 1 and len(S["triangles"]) == 26 and len(S["boundary_legs"]) == 18,
          f"{ex}: S consists of 26 triangles (b_e, b_f, b_C) of the barycentric subdivision; 18 legs")
    Cs = {t["flags"][2][0][1] for t in S["triangles"]}
    check(all(t["flags"][0][0][0] == 1 and t["flags"][1][0][0] == 2 and t["flags"][2][0][0] == 3
              for t in S["triangles"]), f"{ex}: every triangle is of type (edge, two-cell, maximal cell)")
    check(len(Cs) == 3 and {B.lab[c] for c in Cs} == {u0}, f"{ex}: the triangles lie in three maximal cells, "
          "all with facet label u0")
    legset = {(e, f) for e, f, _ in S["boundary_legs"]}
    meets = all((t["flags"][0][0][1], t["flags"][1][0][1]) in legset for t in S["triangles"]
                if len(B.labels(B.F[t["flags"][1][0][1]])) > 1)
    check(meets, f"{ex}: S meets the two-faces of nabla only along its boundary")
    check(all({tuple(t["field"]), tuple(-x for x in t["field"])} == {d0, tuple(-x for x in d0)}
              for t in S["triangles"]), f"{ex}: the field of S is +-d0")
    # naive boundary in the single chart of the interior of F_u0 (all triangles lie in the closed facet)
    bd = defaultdict(lambda: [0, 0, 0, 0])
    for t in S["triangles"]:
        e, f, c = (x[0][1] for x in t["flags"])
        for edge, s in (((2, f, 3, c), 1), ((1, e, 3, c), -1), ((1, e, 2, f), 1)):
            bd[edge] = [a + s * b for a, b in zip(bd[edge], t["field"])]
    nonzero = {k: v for k, v in bd.items() if any(v)}
    legd = {}
    for e, f, cc in S["boundary_legs"]:
        n1, n2 = B.legs[(e, f)]
        legd[(1, e, 2, f)] = [cc * x for x in sub(rays[n2], rays[n1])]
    check(nonzero == legd, f"{ex}: the chain sum N_t [t] with field +-d0 has boundary sum_E c_E d_E [E] on the 18 "
          "legs (a relative tropical 2-cycle in the chart of the facet)")
    st = S["_strands"]
    seq = st[0]["seq"]
    order = [x["node_after"] + 1 for x in seq if x["node_after"] is not None]
    check(order in ([2, 4, 3, 1], [3, 4, 2, 1]), f"{ex}: the boundary meets the nodes in the order q1, q2, q4, q3 "
          "(cyclically, in one of the two directions)")
    # faces F_{u0 u_i} met along the boundary, in order, with the lower supports
    walk = []
    for x in seq:
        for fl in x["faces"]:
            other = [u for u in fl if u != u0]
            walk.append((us.index(other[0]) + 1, frozenset(x["pair"])))
        if x["node_after"] is not None:
            walk.append(("q", x["node_after"] + 1))
    eps = {1: frozenset({6, 3}), 2: frozenset({10, 8}), 3: frozenset({5, 4}), 4: frozenset({2, 7})}
    expect = [(1, eps[1]), (5, eps[1]), ("q", 2), (5, eps[4]), (4, eps[4]), (3, eps[4]), ("q", 4),
              (3, eps[3]), (2, eps[3]), (1, eps[3]), ("q", 3), (1, eps[2]), ("q", 1)]
    rot = lambda L: [L[(i + L.index(("q", 1)) + 1) % len(L)] for i in range(len(L))]  # noqa: E731
    got = rot(walk)
    rev = rot(walk[::-1])
    check(got == expect or rev == expect, f"{ex}: from p_q1 along ([u0,u1],e1), ([u0,u5],e1) to p_q2; "
          "([u0,u5],e4), ([u0,u4],e4), ([u0,u3],e4) to p_q4; ([u0,u3],e3), ([u0,u2],e3), ([u0,u1],e3) to p_q3; "
          "([u0,u1],e2) back to p_q1")
    posf = []
    for (e, f) in legset:
        if B.kind(("e", e), {}) == "positive":
            posf.append(frozenset(B.labels(B.E[e])))
    posf = set(posf)
    want = {frozenset({u0, us[i], us[(i + 1) % 5]}) for i in range(5)}
    check(posf == want, f"{ex}: the five positive vertices of the boundary lie on the edges F_u0 cap F_ui cap F_ui+1")
    check(S["_rho"] == [1, -1, -1, 1], f"{ex}: rho(S) = q1 - q2 - q3 + q4")
    ok(f"{ex}: the facet F_u0, the closed path and the disc S of Section 11.1")


def check_hex_extras(ex, D, rays, nodes, C, Kb, nontri, branch):
    hexf = [f for f in nontri if f["kind"] == "hexagon"][0]
    check(hexf["cyclic"] in rotations_of(D["hexagon_cyclic"]), f"{ex}: the hexagon in cyclic order as printed")
    hq = hexf["nodes"]
    o = hexf["centre"]
    only = all(o not in q for i, q in enumerate(nodes) if i not in hq)
    check(only and all(C[q][o] == -1 for q in hq), f"{ex}: the centre occurs only in the circuits of the hexagon "
          "nodes, with coefficient -1")
    check(all(sum(r[q] for q in hq) == 0 for r in Kb), f"{ex}: the tie sum over the hexagon nodes of lambda_q = 0 "
          "on K")
    if len(hq) == 3:
        proj = [[r[q] for q in hq] for r in Kb]
        plane = [[1, -1, 0], [0, 1, -1]]
        check(same_lattice([p for p in proj if any(p)], plane), f"{ex}: K projects onto the whole tie lattice "
              "{x in Z^3 : x1 + x2 + x3 = 0}")
    # the other tilings of the hexagon
    w, k = hexf["cyclic"], 6
    others = []
    for j in range(3):
        others.append(("A1", [(o, w[(j + i) % k], w[(j + i + 1) % k], w[(j + i + 2) % k]) for i in (0, 3)]))
    for j0 in range(2):
        others.append(("A2", [(o, w[j], w[(j + 1) % k], w[(j + 2) % k]) for j in range(j0, 6, 2)]))
    res = []
    for br, par in others:
        nd = [q for i, q in enumerate(nodes) if i not in hq] + par
        Cn = circuits(nd, len(rays))
        Kn = left_kernel(Cn)
        fo = forced_of(Kn, len(nd))
        res.append((br, len(nd), len(Kn), len(fo)))
    ok(f"{ex}: the five tilings of the hexagon: (branch, |Q|, rank K, forced nodes) = {res}")
    return res


def rotations_of(c):
    out = []
    for s in (c, c[::-1]):
        out += [s[i:] + s[:i] for i in range(len(s))]
    return out


D_cur = None


def run_example(key, path):
    global D_cur
    D = json.load(open(path))
    D_cur = D
    rays = [tuple(r) for r in D["rays"]]
    V = rays[:D["number_of_vertices"]]
    print(f"\n== {D['name']} ({D['polytope']}: {len(V)} vertices; {D['source']})", flush=True)
    Fd, nontri, kinds = check_faces(key, V, rays)
    nodes, branch = check_nodes(key, D, rays, nontri)
    C, Kb, forced = check_K(key, D, nodes, len(rays))
    T = check_triangulation(key, V, Fd)
    check_fan(key, V, Fd, rays, nontri)
    rhos, B = check_cycles(key, D, nodes, C, Kb, rays)
    n, rk = len(nodes), len(Kb)
    # the statements of Section 11 on this example
    if key == "X9":
        check((kinds["square"], kinds["pentagon"], kinds["hexagon"]) == (2, 1, 0) and (n, rk, n - rk) == (4, 1, 3)
              and not forced, f"{key}: two squares, one pentagon; |Q| = 4, rank K = 1, relative Picard number 3; "
              "no node forced")
        pent = [f for f in nontri if f["kind"] == "pentagon"][0]
        check(sorted(pent["cyclic"]) == [3, 4, 5, 6, 8] and pent["centre"] == 10,
              f"{key}: the pentagon has vertices n3, n4, n5, n6, n8 and centre n10")
        check({frozenset(f["cyclic"]) for f in nontri if f["kind"] == "square"} == {frozenset({2, 6, 3, 7}),
              frozenset({2, 5, 4, 7})}, f"{key}: the squares are {{n2,n6,n3,n7}} and {{n2,n5,n4,n7}}")
        check_x9(D, rays, nodes, C, Kb, T, B)
    elif key == "X19":
        check((len(V), kinds["square"], kinds["pentagon"], n, rk, n - rk) == (19, 14, 1, 16, 3, 13) and not forced,
              f"{key}: 19 vertices, 14 squares, one pentagon; |Q| = 16, rank K = 3, |Q| - rank K = 13; no node "
              "forced")
        r = [c["_rho"] for c in D["cycles"]]
        lam = [3 * a + 2 * b + c for a, b, c in zip(*r)]
        check(lam == D["relation_with_all_coefficients_nonzero"] and all(lam),
              f"{key}: lambda = 3 rho(S1) + 2 rho(S2) + rho(S3) has all 16 coefficients nonzero")
    elif key == "X20":
        check((len(V), kinds["square"], kinds["pentagon"], n, rk, n - rk) == (20, 17, 1, 19, 6, 13),
              f"{key}: 20 vertices, 17 squares, one pentagon; |Q| = 19, rank K = 6, |Q| - rank K = 13")
        check(forced == [0, 5], f"{key}: the pentagon nodes q1 and q6 are forced (coordinate 0 on the basis of K)")
    elif key == "X88":
        check((len(V), kinds["square"], kinds["pentagon"], kinds["hexagon"], n, rk) == (21, 22, 0, 1, 24, 10)
              and not forced and branch == "A1", f"{key}: 21 vertices, 22 squares, one hexagon on the A1 branch, no "
              "pentagon; |Q| = 24, rank K = 10; no node forced; nodal locus of codimension 14")
        res = check_hex_extras(key, D, rays, nodes, C, Kb, nontri, branch)
        check(all(r[1:] == (24, 10, 0) for r in res if r[0] == "A1"), f"{key}: all three A1 tilings give |Q| = 24, "
              "rank K = 10 and no forced node")
        shapes = Counter(c["_shape"] for c in D["cycles"])
        check(len(D["cycles"]) == 10 and sum(c["_belt"] for c in D["cycles"]) == 6 and shapes["two vertices joined by three paths"]
              == 2 and shapes["complete graph on four vertices"] == 1, f"{key}: ten cycles: six belt cycles, two "
              f"with two trivalent vertices, one with a complete graph on four vertices, one more ({dict(shapes)})")
        hq = [f for f in nontri if f["kind"] == "hexagon"][0]["nodes"]
        through = [c for c in D["cycles"] if any(c["_rho"][q] for q in hq)]
        check(len(through) == 2 and all(sorted(c["_rho"][q] for q in hq) == [-1, 1] for c in through),
              f"{key}: two of the ten pass through the hexagon nodes, with node coefficients +-(1,-1) there")
        check(D['cycles'][9]['_surface_profile']['interior_gamma_points'] == 18,
              f"{key}: the tenth cycle meets Gamma at 18 interior points")
    elif key == "X154":
        check((len(V), kinds["square"], kinds["pentagon"], kinds["hexagon"], n, rk) == (22, 20, 0, 1, 23, 9)
              and not forced and branch == "A2", f"{key}: 22 vertices, 20 squares, one hexagon on the A2 branch, no "
              "pentagon; |Q| = 23, rank K = 9; no node forced; nodal locus of codimension 14")
        res = check_hex_extras(key, D, rays, nodes, C, Kb, nontri, branch)
        check(all(r[3] == 2 for r in res if r[0] == "A1"), f"{key}: on the A1 branch two nodes are forced")
        check(all(r[3] == 0 for r in res if r[0] == "A2"), f"{key}: on both A2 tilings no node is forced")
        shapes = Counter(c["_shape"] for c in D["cycles"])
        check(len(D["cycles"]) == 9 and sum(c["_belt"] for c in D["cycles"]) == 9, f"{key}: nine belt cycles")
        hq = [f for f in nontri if f["kind"] == "hexagon"][0]["nodes"]
        through = [[c["_rho"][q] for q in hq] for c in D["cycles"] if any(c["_rho"][q] for q in hq)]
        roots = all(sorted(x) == [-1, 0, 1] for x in through)
        check(len(through) == 4 and roots and same_lattice(through, [[1, -1, 0], [0, 1, -1]]),
              f"{key}: four belt cycles pass through two hexagon nodes each, with coefficients roots of the tie "
              f"plane that span it: {through}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--only", nargs="+", choices=["X9", "X19", "X20", "X88", "X154"], default=None)
    A = ap.parse_args()
    keys = ["X9", "X19", "X20", "X88", "X154"]
    if A.only:
        keys = [k for k in keys if k in A.only]
    for k in keys:
        t = time.time()
        run_example(k, HERE / f"{k.lower()}.json")
        print(f"   [{k}: {time.time() - t:.1f} s]", flush=True)
    print(f"\n{'ALL CHECKS PASSED' if not FAILURES else f'{len(FAILURES)} FAILURES'} "
          f"[{time.time() - T0:.1f} s]")
    sys.exit(1 if FAILURES else 0)


if __name__ == "__main__":
    main()
