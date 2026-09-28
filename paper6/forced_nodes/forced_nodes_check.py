#!/usr/bin/env python3
"""Hypothesis (Proj), the relation lattice K and forced nodes on the admissible polytopes of the list.

This program re-checks, from the data file list_polytopes.json.gz alone, the counts of Section 11.4 of the paper
(Table "The admissible polytopes of the list" and the counts for condition (b) of Corollary B).  It uses only the
Python standard library (Python >= 3.10) and exact integer arithmetic throughout.

    python3 paper6/forced_nodes/forced_nodes_check.py            # full check, about a minute
    python3 paper6/forced_nodes/forced_nodes_check.py --sources  # also rebuild the list from the repository files

Notation (Section 2 of the paper).  Delta is a reflexive 4-polytope in N = Z^4.  Delta is admissible if every edge is
primitive, every two-face is a unimodular triangle, a unit square, a dP7 pentagon or a dP6 hexagon, and the dual edge
of every non-triangular two-face has lattice length one.  R is the set of vertices of Delta together with the centres
of its pentagons and hexagons.  A profile sigma subdivides every singular two-face: a square is not subdivided, a
pentagon receives its unique subdivision into two unit parallelograms and a triangle, and a hexagon w_0..w_5 with
centre o receives one of five subdivisions,
    A1:j (j = 0, 1, 2): the parallelograms {o, w_j, w_j+1, w_j+2} and {o, w_j+3, w_j+4, w_j+5}, and two triangles;
    A2:j (j = 0, 1):    the three parallelograms {o, w_i, w_i+1, w_i+2}, i = j, j+2, j+4.
The branch pattern of a profile records which hexagons are on the A1 and which on the A2 branch.  A node
parallelogram q is written (a, b, c, d) with a + c = b + d, the chosen diagonal {a, c} passing through the centre at a
pentagon or hexagon (for a square in cyclic order v_0..v_3 we take {v_0, v_2}); its circuit is
circ_q = e_b + e_d - e_a - e_c in Z^R, and
    K = { lambda in Z^Q : sum_q lambda_q circ_q = 0 }.
A node q is forced if lambda_q = 0 for every lambda in K, i.e. if coordinate q vanishes on a Z-basis of K.  Since a
lattice is not a finite union of sublattices of smaller rank, some lambda in K has every coordinate nonzero if and
only if no node is forced: this is condition (b) of Corollary B.

Hypothesis (Proj) (Lemma "(Proj) as a linear program" of Section 6.6).  For an interior edge eps = [o, w_j] of the
subdivision of a pentagon or hexagon put
    I_eps = e_{w_j-1} + e_{w_j+1} - a_j e_{w_j} + (a_j - 2) e_o,  where (w_j-1 - o) + (w_j+1 - o) = a_j (w_j - o).
(Proj) holds for sigma if and only if some psi in Q^R has <circ_q, psi> = 0 for all q and <I_eps, psi> > 0 for all
interior edges eps at pentagons and hexagons.  The data give, for every positive statement, an integral vector psi
(a witness), and for every negative statement Farkas certificates: rational y_eps >= 0, not all zero, and nu_q with
    sum_eps y_eps I_eps = sum_q nu_q circ_q,
which make the system infeasible (pair with psi).  A Farkas certificate uses only the rows of the squares, the
pentagons and some hexagons with prescribed subdivisions; it therefore excludes every profile that agrees with these
prescriptions.  Such a pair (partial profile, certificate) is called a nogood below, and a negative statement about a
set of profiles is proved by a list of nogoods covering the set.  This program verifies every certificate and every
covering.

What is checked, for every polytope of the data file:
 (1) the facets and two-faces of Delta are recomputed from its vertices; the stored squares, pentagons (with centre)
     and hexagons (with centre) are exactly the non-triangular two-faces, in a cyclic order; Delta is admissible;
     R is the set of vertices and centres;
 (2) the node parallelograms of every stored profile are recomputed from its subdivisions and compared with the
     stored ones; circuits, a Z-basis of K (integer echelon form of [C | I], which gives a saturated basis), rank K
     and the forced nodes are computed;
 (3) every witness psi and every Farkas certificate is verified, and every list of nogoods is shown to cover the set
     of profiles of its statement;
 (4) the counts of Section 11.4 are recomputed and compared with the numbers printed in the paper, and the per
     polytope results with the results stored in the data file (which record the original computation).
It also checks, on every polytope with one hexagon and all five subdivisions of the hexagon, that the forced nodes
depend on the profile only through its branch (Remark "the tie" of the paper), and on the other hexagonal polytopes
it compares, for each profile used, a second profile with the same branch pattern.
The program exits with status 1 on the first failed check.
"""
from __future__ import annotations

import argparse
import ast
import glob
import gzip
import json
import sys
import time
from itertools import combinations, product
from math import gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]                      # the repository root (paper6/forced_nodes/ -> ../..)
DATA = HERE / "list_polytopes.json.gz"
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:6.1f}s] {msg}", flush=True)


class CheckError(Exception):
    pass


def need(cond, msg):
    if not cond:
        raise CheckError(msg)


# ------------------------------------------------------------------------------------------------------------------
# integer linear algebra

def vsub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def vadd(a, b):
    return tuple(x + y for x, y in zip(a, b))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def vgcd(v):
    g = 0
    for x in v:
        g = gcd(g, x)
    return g


def det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def normal4(d1, d2, d3):
    """A vector orthogonal to d1, d2, d3 in Z^4 (generalised cross product); zero iff they are dependent."""
    rows = [d1, d2, d3]
    n = []
    for k in range(4):
        cols = [j for j in range(4) if j != k]
        m = [[r[j] for j in cols] for r in rows]
        n.append((-1) ** k * det3(m))
    return tuple(n)


def rank_int(rows):
    """Rank over Q of a list of integer vectors (fraction-free elimination)."""
    A = [list(r) for r in rows if any(r)]
    if not A:
        return 0
    ncol = len(A[0])
    r = 0
    for c in range(ncol):
        piv = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        for i in range(r + 1, len(A)):
            if A[i][c]:
                f, g = A[i][c], A[r][c]
                A[i] = [g * x - f * y for x, y in zip(A[i], A[r])]
                h = vgcd(A[i])
                if h > 1:
                    A[i] = [x // h for x in A[i]]
        r += 1
        if r == len(A):
            break
    return r


def left_kernel_basis(C):
    """Z-basis of {x in Z^m : x C = 0} for an integer m x n matrix C (rows).  Integer row echelon form of [C | I] by
    unimodular row operations; the identity part of the rows whose C part vanishes is a basis of the left kernel, and
    it is saturated because the transformation is unimodular."""
    m = len(C)
    n = len(C[0]) if m else 0
    A = [list(C[i]) + [1 if j == i else 0 for j in range(m)] for i in range(m)]
    r = 0
    for c in range(n):
        rows = [i for i in range(r, m) if A[i][c] != 0]
        if not rows:
            continue
        while True:                                   # Euclid on column c among rows r..m-1
            rows = [i for i in range(r, m) if A[i][c] != 0]
            p = min(rows, key=lambda i: abs(A[i][c]))
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
        r += 1
        if r == m:
            break
    basis = [A[i][n:] for i in range(m) if not any(A[i][:n])]
    for b in basis:                                   # sanity: x C = 0
        need(all(sum(b[i] * C[i][j] for i in range(m)) == 0 for j in range(n)), "kernel vector")
    return basis


# ------------------------------------------------------------------------------------------------------------------
# faces of a reflexive 4-polytope

def facets(V):
    """Facets of Delta = conv(V) (0 in the interior): list of (u, frozenset of vertex indices), u integral with
    <u, x> <= 1 on Delta and = 1 on the facet (so u is a vertex of the polar polytope).  Every 4-subset of vertices
    spanning a hyperplane is tested (skipping subsets of facets already found)."""
    n = len(V)
    found = []
    masks = []
    for S in combinations(range(n), 4):
        sm = (1 << S[0]) | (1 << S[1]) | (1 << S[2]) | (1 << S[3])
        if any(sm & ~fm == 0 for fm in masks):
            continue
        p0 = V[S[0]]
        n0, n1, n2, n3 = normal4(vsub(V[S[1]], p0), vsub(V[S[2]], p0), vsub(V[S[3]], p0))
        if not (n0 or n1 or n2 or n3):
            continue
        c = n0 * p0[0] + n1 * p0[1] + n2 * p0[2] + n3 * p0[3]
        lo = hi = False
        for v in V:
            x = n0 * v[0] + n1 * v[1] + n2 * v[2] + n3 * v[3]
            if x > c:
                hi = True
            elif x < c:
                lo = True
            if lo and hi:
                break
        if lo and hi:
            continue
        nv = (n0, n1, n2, n3)
        if hi:
            nv, c = tuple(-x for x in nv), -c
        need(c > 0, "origin not interior")
        need(all(x % c == 0 for x in nv), "facet not at lattice distance one (not reflexive)")
        u = tuple(x // c for x in nv)
        F = frozenset(i for i in range(n) if dot(u, V[i]) == 1)
        found.append((u, F))
        masks.append(sum(1 << i for i in F))
    return found


def affine_rank(P):
    if not P:
        return 0
    p0 = P[0]
    return rank_int([vsub(p, p0) for p in P[1:]]) + 1


def cyclic_order(P):
    """The vertices of a convex lattice polygon in R^4, in cyclic order (exact comparator in a coordinate projection
    that is injective on the plane of the polygon)."""
    k = len(P)
    d1, d2 = None, None
    for i, j in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)):
        a, b = vsub(P[1], P[0]), None
        for p in P[2:]:
            b = vsub(p, P[0])
            if a[i] * b[j] - a[j] * b[i] != 0:
                d1, d2 = i, j
                break
        if d1 is not None:
            break
    need(d1 is not None, "degenerate polygon")
    S = [sum(p[t] for p in P) for t in range(4)]
    pts = [(k * p[d1] - S[d1], k * p[d2] - S[d2]) for p in P]

    def half(v):
        return 0 if (v[1] > 0 or (v[1] == 0 and v[0] > 0)) else 1

    from functools import cmp_to_key

    def cmp(i, j):
        a, b = pts[i], pts[j]
        ha, hb = half(a), half(b)
        if ha != hb:
            return ha - hb
        cr = a[0] * b[1] - a[1] * b[0]
        return -1 if cr > 0 else (1 if cr < 0 else 0)

    order = sorted(range(k), key=cmp_to_key(cmp))
    cyc = [P[i] for i in order]
    # convexity: every consecutive turn has the same strict sign
    q = [pts[i] for i in order]
    turns = [(q[(i + 1) % k][0] - q[i][0]) * (q[(i + 2) % k][1] - q[i][1])
             - (q[(i + 1) % k][1] - q[i][1]) * (q[(i + 2) % k][0] - q[i][0]) for i in range(k)]
    need(all(t > 0 for t in turns), "polygon vertices not in convex position")
    return cyc, (d1, d2)


def normalise_cycle(cyc):
    """Rotate/reflect a cycle so that w_0 is the lexicographically smallest vertex and w_1 < w_{k-1}."""
    k = len(cyc)
    i = cyc.index(min(cyc))
    c = cyc[i:] + cyc[:i]
    if c[1] > c[-1]:
        c = [c[0]] + c[1:][::-1]
    return c


def area2(cyc):
    """Twice the lattice area (normalised area) of a lattice polygon in Z^4, in the lattice of its plane: the sum over
    a fan triangulation of the gcd of the 2x2 minors (the plane lattice is saturated, so its minors are coprime)."""
    tot = 0
    p0 = cyc[0]
    for i in range(1, len(cyc) - 1):
        a, b = vsub(cyc[i], p0), vsub(cyc[i + 1], p0)
        tot += vgcd([a[s] * b[t] - a[t] * b[s] for s, t in combinations(range(4), 2)])
    return tot


def inside_strict(cyc, proj, o):
    d1, d2 = proj
    k = len(cyc)
    q = [(p[d1], p[d2]) for p in cyc]
    oo = (o[d1], o[d2])
    s = [(q[(i + 1) % k][0] - q[i][0]) * (oo[1] - q[i][1]) - (q[(i + 1) % k][1] - q[i][1]) * (oo[0] - q[i][0])
         for i in range(k)]
    return all(x > 0 for x in s) or all(x < 0 for x in s)


def analyse_faces(V):
    """Facets, two-faces and admissibility of Delta.  Returns dict(admissible, squares, pentagons, hexagons, reason)
    with squares as normalised cycles and pentagons/hexagons as dict(cycle, centre)."""
    V = [tuple(v) for v in V]
    F = facets(V)
    faces = {}
    for (ua, Ia), (ub, Ib) in combinations(F, 2):
        I = Ia & Ib
        if len(I) < 3 or I in faces:
            continue
        P = [V[i] for i in sorted(I)]
        if affine_rank(P) != 3:
            continue
        faces[I] = (ua, ub)
    out = dict(admissible=True, squares=[], pentagons=[], hexagons=[], reason=None, nfacets=len(F),
               ntwofaces=len(faces))
    for I, (ua, ub) in faces.items():
        P = [V[i] for i in sorted(I)]
        cyc, proj = cyclic_order(P)
        k = len(cyc)
        if any(vgcd(vsub(cyc[(i + 1) % k], cyc[i])) != 1 for i in range(k)):
            out.update(admissible=False, reason="edge not primitive")
            return out
        A2 = area2(cyc)
        kind = {(3, 1): "triangle", (4, 2): "square", (5, 5): "pentagon", (6, 6): "hexagon"}.get((k, A2))
        if kind is None:
            out.update(admissible=False, reason=f"two-face with {k} vertices and normalised area {A2}")
            return out
        if kind == "triangle":
            continue
        if vgcd(vsub(ua, ub)) != 1:
            out.update(admissible=False, reason="dual edge of length > 1")
            return out
        cyc = normalise_cycle(cyc)
        if kind == "square":
            need(vadd(cyc[0], cyc[2]) == vadd(cyc[1], cyc[3]), "square not a parallelogram")
            out["squares"].append(cyc)
            continue
        cands = {vsub(vadd(cyc[i], cyc[(i + 2) % k]), cyc[(i + 1) % k]) for i in range(k)}
        inner = [o for o in cands if inside_strict(cyc, proj, o)]
        need(len(inner) == 1, "centre of a pentagon or hexagon")
        out[kind + "s"].append(dict(cycle=cyc, centre=inner[0]))
    for key in ("squares",):
        out[key].sort()
    out["pentagons"].sort(key=lambda d: d["cycle"])
    out["hexagons"].sort(key=lambda d: d["cycle"])
    return out


# ------------------------------------------------------------------------------------------------------------------
# profiles, nodes, circuits, (Proj) rows

HEX_TILINGS = ["A1:0", "A1:1", "A1:2", "A2:0", "A2:1"]


def tiling_starts(code, k):
    """Start indices i of the parallelograms {o, w_i, w_i+1, w_i+2} of a hexagon subdivision; the interior edges of
    the subdivision are [o, w_s] for the start indices s of all its tiles."""
    br, j = code.split(":")
    j = int(j)
    if br == "A1":
        return [j, j + 3], [j, j + 2, j + 3, (j + 5) % 6]
    return [j, j + 2, j + 4], [j, j + 2, j + 4]


def pentagon_tiling(face):
    w, o = face["cycle"], face["centre"]
    good = [i for i in range(5) if vadd(w[i], w[(i + 2) % 5]) == vadd(o, w[(i + 1) % 5])]
    starts = [i for i in range(5) if i in good and (i + 2) % 5 in good]
    need(len(starts) == 1, "the pentagon subdivision is unique")
    i = starts[0]
    return [i, (i + 2) % 5], [i, (i + 2) % 5, (i + 4) % 5]


class Polytope:
    """Faces, R and the rows of the (Proj) system of an admissible polytope, as stored in the data file."""

    def __init__(self, rec):
        self.rec = rec
        self.V = [tuple(v) for v in rec["vertices"]]
        self.squares = [[tuple(p) for p in c] for c in rec["squares"]]
        self.pent = [dict(cycle=[tuple(p) for p in f["cycle"]], centre=tuple(f["centre"])) for f in rec["pentagons"]]
        self.hexes = [dict(cycle=[tuple(p) for p in f["cycle"]], centre=tuple(f["centre"])) for f in rec["hexagons"]]
        self.R = [tuple(p) for p in rec["R"]]
        self.ix = {p: i for i, p in enumerate(self.R)}

    def check_faces(self):
        fa = analyse_faces(self.V)
        need(fa["admissible"], f"not admissible: {fa['reason']}")
        need(sorted(fa["squares"]) == sorted(self.squares), "squares")
        need(sorted((f["cycle"], f["centre"]) for f in fa["pentagons"])
             == sorted((f["cycle"], f["centre"]) for f in self.pent), "pentagons")
        need(sorted((f["cycle"], f["centre"]) for f in fa["hexagons"])
             == sorted((f["cycle"], f["centre"]) for f in self.hexes), "hexagons")
        R = sorted(set(self.V) | {f["centre"] for f in self.pent + self.hexes})
        need(R == sorted(self.R) and len(R) == len(self.R), "R = vertices and centres")
        for f in self.hexes:                          # every w_i + w_i+2 = o + w_i+1 in a dP6 hexagon
            w, o = f["cycle"], f["centre"]
            need(all(vadd(w[i], w[(i + 2) % 6]) == vadd(o, w[(i + 1) % 6]) for i in range(6)), "hexagon relations")
        return fa

    # nodes: ids "S<k>", "P<f>:<i>", "H<f>:<i>"; interior edges: "P<f>:e<j>", "H<f>:e<j>"
    def nodes(self, profile):
        """Node parallelograms of a profile (dict hexagon index -> subdivision code), as (id, (a, b, c, d))."""
        out = []
        for k, c in enumerate(self.squares):
            out.append((f"S{k}", (c[0], c[1], c[2], c[3])))
        for f, face in enumerate(self.pent):
            starts, _ = pentagon_tiling(face)
            out += [(f"P{f}:{i}", self.par(face, i)) for i in starts]
        for f, face in enumerate(self.hexes):
            if str(f) in profile:
                starts, _ = tiling_starts(profile[str(f)], 6)
                out += [(f"H{f}:{i % 6}", self.par(face, i)) for i in starts]
        return out

    @staticmethod
    def par(face, i):
        w, o, k = face["cycle"], face["centre"], len(face["cycle"])
        a, b, c, d = o, w[i % k], w[(i + 1) % k], w[(i + 2) % k]
        need(vadd(a, c) == vadd(b, d), "node parallelogram: a + c = b + d")
        return (a, b, c, d)

    def circuit(self, abcd):
        a, b, c, d = abcd
        v = [0] * len(self.R)
        v[self.ix[b]] += 1
        v[self.ix[d]] += 1
        v[self.ix[a]] -= 1
        v[self.ix[c]] -= 1
        return v

    def edge_row(self, face, j):
        w, o, k = face["cycle"], face["centre"], len(face["cycle"])
        s = vadd(vsub(w[(j - 1) % k], o), vsub(w[(j + 1) % k], o))
        t = vsub(w[j % k], o)
        l = next(x for x in t if x != 0)
        a = s[t.index(l)] // l
        need(all(x == a * y for x, y in zip(s, t)), "(w_j-1 - o) + (w_j+1 - o) = a_j (w_j - o)")
        v = [0] * len(self.R)
        v[self.ix[w[(j - 1) % k]]] += 1
        v[self.ix[w[(j + 1) % k]]] += 1
        v[self.ix[w[j % k]]] -= a
        v[self.ix[o]] += a - 2
        return v

    def edges(self, profile):
        """Interior edges of the subdivisions of the pentagons and of the hexagons assigned by `profile`."""
        out = []
        for f, face in enumerate(self.pent):
            _, es = pentagon_tiling(face)
            out += [(f"P{f}:e{j}", self.edge_row(face, j)) for j in es]
        for f, face in enumerate(self.hexes):
            if str(f) in profile:
                _, es = tiling_starts(profile[str(f)], 6)
                out += [(f"H{f}:e{j % 6}", self.edge_row(face, j)) for j in es]
        return out

    def lattice(self, profile):
        """(nodes, Z-basis of K, forced node ids) for a complete profile."""
        nd = self.nodes(profile)
        C = [self.circuit(P) for _, P in nd]
        Kb = left_kernel_basis(C) if C else []
        forced = [nd[q][0] for q in range(len(nd)) if all(b[q] == 0 for b in Kb)]
        return nd, Kb, forced


def branch_pattern(profile, hexes):
    return tuple(profile[str(h)].split(":")[0] for h in hexes)


def owners(ids):
    """Hexagon indices occurring in a list of node or edge ids."""
    return {int(s[1:].split(":")[0]) for s in ids if s.startswith("H")}


# ------------------------------------------------------------------------------------------------------------------
# verification of witnesses, Farkas certificates and coverings

def check_witness(P, profile, psi):
    need(set(profile) == {str(i) for i in range(len(P.hexes))} and all(c in HEX_TILINGS for c in profile.values()), "witness: complete valid profile")
    need(len(psi) == len(P.R) and all(type(x) is int for x in psi), "witness: integral vector on R")
    for _, abcd in P.nodes(profile):
        need(dot(P.circuit(abcd), psi) == 0, "witness: <circ_q, psi> = 0")
    for _, row in P.edges(profile):
        need(dot(row, psi) > 0, "witness: <I_eps, psi> > 0")


def check_farkas(P, nogood):
    """A nogood {"profile": partial profile, "y": {edge id: int >= 0}, "nu": {node id: int}}: sum y I = sum nu circ,
    y >= 0 not all zero, supported on the rows of the squares, pentagons and the hexagons of the partial profile."""
    prof = nogood["profile"]
    need(set(prof) <= {str(i) for i in range(len(P.hexes))} and all(c in HEX_TILINGS for c in prof.values()), "Farkas: valid partial profile")
    E = dict(P.edges(prof))
    N = dict(P.nodes(prof))
    y, nu = nogood["y"], nogood["nu"]
    need(all(e in E for e in y) and all(q in N for q in nu), "Farkas: rows outside the partial profile")
    need(all(type(v) is int and v >= 0 for v in y.values()) and any(v > 0 for v in y.values()),
         "Farkas: y >= 0, not zero")
    need(all(type(v) is int for v in nu.values()), "Farkas: integral coefficients")
    lhs = [0] * len(P.R)
    for e, v in y.items():
        lhs = [s + v * t for s, t in zip(lhs, E[e])]
    rhs = [0] * len(P.R)
    for q, v in nu.items():
        rhs = [s + v * t for s, t in zip(rhs, P.circuit(N[q]))]
    need(lhs == rhs, "Farkas: sum y I_eps = sum nu circ_q")


def check_cover(P, claim, allowed):
    """Every complete profile whose branch pattern (on the hexagons in claim['order']) is allowed agrees with the
    partial profile of some nogood.  Depth-first over the hexagons in the given order, pruning a partial profile as
    soon as it contains the partial profile of a nogood or no allowed pattern extends it.  `allowed(pat, full)`
    decides, for a full pattern, whether it belongs to the statement, and for a partial one whether some extension
    does."""
    order = claim["order"]
    need(sorted(order) == list(range(len(P.hexes))), "order of the hexagons")
    nogoods = [{int(h): c for h, c in ng["profile"].items()} for ng in claim["nogoods"]]
    for ng in claim["nogoods"]:
        check_farkas(P, ng)

    def rec(depth, assign):
        if any(all(assign.get(h) == c for h, c in ng.items()) for ng in nogoods):
            return
        pat = tuple(assign[h].split(":")[0] for h in order[:depth])
        full = depth == len(order)
        if not allowed(pat, full):
            return
        need(not full, "a profile of the statement is not excluded by any nogood")
        h = order[depth]
        for c in HEX_TILINGS:
            assign[h] = c
            rec(depth + 1, assign)
            del assign[h]

    rec(0, {})


def pattern_filter(kind, good=None):
    """allowed(pat, full) for the statements: kind 'A1' (every hexagon on A1), 'A2', 'mixed' (both branches occur)
    or None (any), intersected with a set `good` of full patterns if given."""
    def allowed(pat, full):
        if kind in ("A1", "A2") and any(b != kind for b in pat):
            return False
        if kind == "mixed" and full and len(set(pat)) < 2:
            return False
        if good is not None:
            return (pat in good) if full else any(g[:len(pat)] == pat for g in good)
        return True
    return allowed


# ------------------------------------------------------------------------------------------------------------------
# the list from the repository files (optional)

def list_from_sources():
    """The vertex sets of the list (Section 11.4), rebuilt from the repository: the 77 polytopes of
    paper4/certificates/framework_77.json, Delta_19 and Delta_20 (paper4/certificates/examples.py), the polytope X°
    of Paper 4 (the polar of V_F1 in paper4/certificates/examples.py), the hits of the scan of Paper 3 in
    output/both_sides_v*.json, and the 36-vertex product of two hexagons (src/missing_polytope.py).  Returns the list of
    (source label, vertex list) in this order; equal vertex sets are merged later."""
    portable = HERE / "source_vertices.json"
    if portable.exists():
        return [(label, [tuple(v) for v in V]) for label, V in json.loads(portable.read_text())]
    out = []
    fw = json.load(open(ROOT / "paper4/certificates/framework_77.json"))
    for i, r in enumerate(fw):
        out.append((f"paper4/certificates/framework_77.json[{i}]", [tuple(v) for v in r["V"]]))
    tree = ast.parse(open(ROOT / "paper4/certificates/examples.py").read())
    consts = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            if node.targets[0].id in ("V_19", "V_20", "V_F1"):
                consts[node.targets[0].id] = [tuple(v) for v in ast.literal_eval(node.value)]
    out.append(("paper4/certificates/examples.py:V_19", consts["V_19"]))
    out.append(("paper4/certificates/examples.py:V_20", consts["V_20"]))
    out.append(("polar of paper4/certificates/examples.py:V_F1", [u for u, _ in facets(consts["V_F1"])]))
    n = 0
    for f in sorted(glob.glob(str(ROOT / "output/both_sides_v*.json"))):
        for key, r in json.load(open(f)).items():
            for j, h in enumerate(r["hits"]):
                out.append((f"output/{Path(f).name}:{key}:hit {j} (list number {n})", [tuple(v) for v in h["verts"]]))
                n += 1
    tree = ast.parse(open(ROOT / "src/missing_polytope.py").read())
    for node in tree.body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", None) == "DELTA36":
            out.append((f"src/missing_polytope.py:DELTA36 (list number {n})",
                        [tuple(v) for v in ast.literal_eval(node.value)]))
    return out


def check_sources(D):
    src = list_from_sources()
    merged = {}
    for label, V in src:
        merged.setdefault(frozenset(V), []).append(label)
    log(f"list from the repository files: {len(src)} entries, {len(merged)} distinct vertex sets")
    need(len(merged) == 668, "668 vertex sets")
    adm = {}
    for key, labels in merged.items():
        fa = analyse_faces(sorted(key))
        if fa["admissible"]:
            adm[key] = labels
    log(f"admissible vertex sets: {len(adm)}")
    need(len(adm) == 613, "613 admissible vertex sets")
    stored = {frozenset(tuple(v) for v in r["vertices"]): r for r in D["polytopes"]}
    need(set(stored) == set(adm), "the data file contains exactly the admissible vertex sets of the list")
    for key, r in stored.items():
        need(sorted(r["sources"]) == sorted(adm[key]), f"sources of {r['id']}")
    log("the data file lists exactly the 613 admissible vertex sets of the list, with their sources")


def check_equivalence(rec, D):
    """X° and the polytope it is equivalent to: an explicit matrix g in GL_4(Z) maps one vertex set onto the other."""
    eq = rec.get("equivalent_to")
    if not eq:
        return
    other = next(r for r in D["polytopes"] if r["id"] == eq["id"])
    g = [tuple(r) for r in eq["matrix"]]
    d = sum(g[0][j] * (-1) ** j * det3([[g[i][c] for c in range(4) if c != j] for i in (1, 2, 3)]) for j in range(4))
    need(abs(d) == 1, "equivalence: det = +-1")
    img = {tuple(sum(g[i][j] * v[j] for j in range(4)) for i in range(4)) for v in rec["vertices"]}
    need(img == {tuple(v) for v in other["vertices"]}, "equivalence: g maps the vertex set onto the other")


# ------------------------------------------------------------------------------------------------------------------
# main check

EXPECTED = dict(hexfree=470, hexfree_smooth=8, hexfree_K0=87, hexfree_Knz=375, hexfree_proj=364, hexfree_projfail=11,
                hex=142, hex_proj=83, hex_noproj=59, hex_A1=20, hex_A2=55, hex_mixed=18,
                corB_hexfree=118, corB_hex=43, hex_one=29, corB_hex_one=10, hex_many=54, corB_hex_many=33,
                corB_A1=18, corB_A2=22, corB_mixed=10, scope_all=447, corB_all=161, corB_fail=286)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", default=str(DATA))
    ap.add_argument("--summary", type=Path, help="write verified counts to JSON after a successful check")
    ap.add_argument("--sources", action="store_true", help="also rebuild the list from the repository files")
    A = ap.parse_args()
    D = json.load(gzip.open(A.data, "rt"))
    polys = D["polytopes"]
    log(f"{len(polys)} admissible vertex sets in {Path(A.data).name}")
    if A.sources:
        check_sources(D)

    cnt = {k: 0 for k in EXPECTED}
    ntilecheck = [0, 0]
    ncert = [0, 0]
    for idx, rec in enumerate(polys):
        name = rec["id"]
        try:
            P = Polytope(rec)
            P.check_faces()
            check_equivalence(rec, D)
            st = rec["stored"]
            if rec.get("equivalent_to"):
                # counted once, as the polytope it is equivalent to; its (Proj) statements are still checked
                pass
            if not P.hexes:
                prof = {}
                nd, Kb, forced = P.lattice(prof)
                sing = bool(P.squares or P.pent)
                need(len(nd) == st["nodes"] and len(Kb) == st["rank_K"], "stored nodes / rank K")
                need([n for n, _ in rec["nodes"]] == [n for n, _ in nd]
                     and all(tuple(map(tuple, a)) == b for (_, a), (_, b) in zip(rec["nodes"], nd)),
                     "stored node parallelograms")
                cnt["hexfree"] += 1
                if not sing:
                    cnt["hexfree_smooth"] += 1
                    continue
                if not Kb:
                    cnt["hexfree_K0"] += 1
                    continue
                cnt["hexfree_Knz"] += 1
                pj = rec["proj"]
                if "witness" in pj:
                    check_witness(P, prof, pj["witness"])
                    ncert[0] += 1
                    cnt["hexfree_proj"] += 1
                    ok = not forced
                    need(ok == st["no_forced"] and len(forced) == st["forced"], "stored forced nodes")
                    cnt["corB_hexfree"] += ok
                else:
                    check_farkas(P, {"profile": {}, **pj["farkas"]})
                    ncert[1] += 1
                    cnt["hexfree_projfail"] += 1
                    need(st["proj"] is False, "stored (Proj)")
                continue

            # ---- polytopes with a hexagon
            hx = list(range(len(P.hexes)))
            nh = len(hx)
            main_one = not rec.get("equivalent_to")
            if main_one:
                cnt["hex"] += 1
            feas = {}
            for kind in ("A1", "A2", "mixed"):
                cl = rec["kinds"][kind]
                if "witness" in cl:
                    w = cl["witness"]
                    check_witness(P, w["profile"], w["psi"])
                    ncert[0] += 1
                    pat = branch_pattern(w["profile"], hx)
                    need((kind == "mixed" and len(set(pat)) == 2) or set(pat) == {kind}, "witness of the kind")
                    feas[kind] = w
                else:
                    ncert[1] += len(cl["nogoods"])
                    check_cover(P, cl, pattern_filter(kind))
                    feas[kind] = None
                need(bool(feas[kind]) == st["feasible"][kind], f"stored (Proj) feasibility, {kind}")
            anyproj = any(feas.values())
            if not main_one:
                need(not anyproj, "X° fails (Proj) on every profile")
                continue
            cnt["hex_proj"] += anyproj
            cnt["hex_noproj"] += not anyproj
            for kind in ("A1", "A2", "mixed"):
                cnt["hex_" + kind] += bool(feas[kind])

            # forced nodes depend only on the branch pattern
            if nh == 1:
                per = {}
                for c in HEX_TILINGS:
                    nd, Kb, forced = P.lattice({"0": c})
                    per.setdefault(c[:2], set()).add(
                        (len(nd), len(Kb), tuple(sorted(x for x in forced if not x.startswith("H"))),
                         sum(x.startswith("H") for x in forced)))
                need(all(len(v) == 1 for v in per.values()), "forced nodes depend on the tiling only through the branch")
                ntilecheck[0] += 1

            def forced_of(profile):
                nd, Kb, forced = P.lattice(profile)
                if nh > 1:                        # a second profile with the same branch pattern
                    alt = {h: {"A1:0": "A1:1", "A1:1": "A1:2", "A1:2": "A1:0", "A2:0": "A2:1", "A2:1": "A2:0"}[c]
                           for h, c in profile.items()}
                    nd2, Kb2, f2 = P.lattice(alt)
                    key = lambda n, K, f: (len(n), len(K), sorted(x for x in f if not x.startswith("H")),
                                           sorted(owners([x]).pop() for x in f if x.startswith("H")))
                    need(key(nd, Kb, forced) == key(nd2, Kb2, f2), "forced nodes: second profile, same pattern")
                    ntilecheck[1] += 1
                return nd, Kb, forced

            if not anyproj:
                need(not st["corB"]["any"], "stored condition (b)")
                continue
            # Corollary B, any profile satisfying (Proj)
            cb = rec["corollary_b"]
            if "witness" in cb:
                w = cb["witness"]
                check_witness(P, w["profile"], w["psi"])
                ncert[0] += 1
                nd, Kb, forced = forced_of(w["profile"])
                need(not forced, "condition (b): the witness profile has no forced node")
                need([n for n, _ in w["nodes"]] == [n for n, _ in nd]
                     and all(tuple(map(tuple, a)) == b for (_, a), (_, b) in zip(w["nodes"], nd)),
                     "stored node parallelograms of the witness profile")
                okB = True
            else:
                # every profile satisfying (Proj) has a forced node: the stored (Proj) witnesses have forced nodes,
                # and the profiles with a pattern without forced node are excluded by nogoods
                for kind in ("A1", "A2", "mixed"):
                    if feas[kind]:
                        _, _, fo = forced_of(feas[kind]["profile"])
                        need(fo, "a (Proj) witness profile without forced node")
                good = []
                for pat in product(("A1", "A2"), repeat=nh):
                    prof = {str(h): ("A1:0" if b == "A1" else "A2:0") for h, b in zip(hx, pat)}
                    if not P.lattice(prof)[2]:
                        good.append(pat)
                need(sorted(map(list, good)) == sorted(cb["good_patterns"]), "patterns without forced node")
                gset = [tuple(g) for g in good]
                ordr = cb["order"]
                gset_o = [tuple(g[h] for h in ordr) for g in gset]
                ncert[1] += len(cb["nogoods"])
                check_cover(P, cb, pattern_filter(None, gset_o))
                okB = False
            need(okB == st["corB"]["any"], "stored condition (b), some profile")
            cnt["corB_hex"] += okB
            if nh == 1:
                cnt["hex_one"] += 1
                cnt["corB_hex_one"] += okB
            else:
                cnt["hex_many"] += 1
                cnt["corB_hex_many"] += okB
            # by branch
            for kind in ("A1", "A2"):
                if feas[kind]:
                    _, _, fo = forced_of(feas[kind]["profile"])
                    ok = not fo
                    need(ok == st["corB"][kind], f"stored condition (b), {kind}")
                    cnt["corB_" + kind] += ok
            if feas["mixed"]:
                cm = rec["corollary_b_mixed"]
                if "witness" in cm:
                    w = cm["witness"]
                    check_witness(P, w["profile"], w["psi"])
                    pat = branch_pattern(w["profile"], hx)
                    need(len(set(pat)) == 2, "mixed witness")
                    _, _, fo = forced_of(w["profile"])
                    need(not fo, "mixed witness without forced node")
                    ok = True
                else:
                    _, _, fo = forced_of(feas["mixed"]["profile"])
                    need(fo, "the mixed (Proj) witness has a forced node")
                    good = []
                    for pat in product(("A1", "A2"), repeat=nh):
                        if len(set(pat)) < 2:
                            continue
                        prof = {str(h): ("A1:0" if b == "A1" else "A2:0") for h, b in zip(hx, pat)}
                        if not P.lattice(prof)[2]:
                            good.append(pat)
                    need(sorted(map(list, good)) == sorted(cm["good_patterns"]), "mixed patterns without forced node")
                    ordr = cm["order"]
                    ncert[1] += len(cm["nogoods"])
                    check_cover(P, cm, pattern_filter("mixed", [tuple(g[h] for h in ordr) for g in good]))
                    ok = False
                need(ok == st["corB"]["mixed"], "stored condition (b), mixed")
                cnt["corB_mixed"] += ok
        except CheckError as e:
            print(f"FAILED on {name}: {e}", flush=True)
            sys.exit(1)
        if (idx + 1) % 100 == 0:
            log(f"{idx + 1} polytopes checked")

    cnt["scope_all"] = cnt["hexfree_proj"] + cnt["hex_proj"]
    cnt["corB_all"] = cnt["corB_hexfree"] + cnt["corB_hex"]
    cnt["corB_fail"] = cnt["scope_all"] - cnt["corB_all"]
    log(f"(Proj) witnesses verified: {ncert[0]}; Farkas certificates verified: {ncert[1]}")
    log(f"branch-only dependence of forced nodes: {ntilecheck[0]} polytopes with one hexagon (all five "
        f"subdivisions), {ntilecheck[1]} comparisons of two profiles with the same pattern")
    print()
    print("Table (admissible polytopes of the list, up to GL_4(Z))")
    print(f"  no hexagonal two-face:        {cnt['hexfree']}")
    print(f"    no singular two-face:       {cnt['hexfree_smooth']}")
    print(f"    singular two-face, K = 0:   {cnt['hexfree_K0']}")
    print(f"    K != 0:                     {cnt['hexfree_Knz']}   (Proj): {cnt['hexfree_proj']}, "
          f"fails: {cnt['hexfree_projfail']}")
    print(f"  some hexagonal two-face:      {cnt['hex']}   (Proj) for some profile: {cnt['hex_proj']}, "
          f"for none: {cnt['hex_noproj']}")
    print(f"    every hexagon A1 / every hexagon A2 / mixed: {cnt['hex_A1']} / {cnt['hex_A2']} / {cnt['hex_mixed']}")
    print("Condition (b) of Corollary B (no forced node) for some profile satisfying (Proj)")
    print(f"  hexagon-free:  {cnt['corB_hexfree']} of {cnt['hexfree_proj']}")
    print(f"  hexagonal:     {cnt['corB_hex']} of {cnt['hex_proj']}  (one hexagon {cnt['corB_hex_one']} of "
          f"{cnt['hex_one']}; several {cnt['corB_hex_many']} of {cnt['hex_many']})")
    print(f"  by branch:     A1 {cnt['corB_A1']} of {cnt['hex_A1']}, A2 {cnt['corB_A2']} of {cnt['hex_A2']}, "
          f"mixed {cnt['corB_mixed']} of {cnt['hex_mixed']}")
    print(f"  in all:        {cnt['corB_all']} of {cnt['scope_all']}; every profile with (Proj) has a forced node "
          f"on {cnt['corB_fail']}")
    bad = {k: (cnt[k], v) for k, v in EXPECTED.items() if cnt[k] != v}
    if bad:
        print("MISMATCH with the numbers printed in the paper:", bad)
        sys.exit(1)
    if A.summary:
        import hashlib
        A.summary.write_text(json.dumps(dict(counts=cnt, vertex_sets=len(polys),
            projectivity_witnesses=ncert[0], farkas_certificates=ncert[1],
            data_sha256=hashlib.sha256(Path(A.data).read_bytes()).hexdigest()), indent=2) + "\n")
    log("all counts agree with the paper")


if __name__ == "__main__":
    main()
