#!/usr/bin/env python3
"""
Batyrev Hodge numbers of the MPCP resolution of the generic anticanonical
hypersurface X in P_Delta, for reflexive Delta in Z^4 (N-side / fan polytope):

  h^{1,1}(X^) = l(Delta)  - 5 - sum_{facets G}  l*(G)  + sum_{2-faces F}  l*(F) l*(F*)
  h^{2,1}(X^) = l(Delta*) - 5 - sum_{facets G*} l*(G*) + sum_{2-faces F*} l*(F*) l*(F)

(Batyrev, J. Alg. Geom. 3 (1994); our Delta is his Delta*.)  Anchors: quintic
simplex -> (1, 101); P(1,1,1,3,3) -> (4, 112).  Headline polytopes of the
paper: Delta_A -> (5, 101), chi = -192; Delta_B -> (9, 81), chi = -144.
"""
import os
import sys
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from batyrev_global import (facets, two_faces, dot, vsub, vgcd, affine_rank)


def lattice_points(V, facs):
    los = [min(v[i] for v in V) for i in range(4)]
    his = [max(v[i] for v in V) for i in range(4)]
    pts = []
    for x in range(los[0], his[0] + 1):
        for y in range(los[1], his[1] + 1):
            for z in range(los[2], his[2] + 1):
                for w in range(los[3], his[3] + 1):
                    p = (x, y, z, w)
                    if all(dot(u, p) <= c for u, c, _ in facs):
                        pts.append(p)
    return pts


def face_interior_count(pts, facs, tight_target):
    """#lattice points whose tight facet-set equals tight_target."""
    n = 0
    for p in pts:
        tight = frozenset(i for i, (u, c, _) in enumerate(facs) if dot(u, p) == c)
        if tight == tight_target:
            n += 1
    return n


def hodge_numbers(V, name=""):
    facs = facets(V)
    assert all(c == 1 for _, c, _ in facs), "not reflexive"
    Vd = [u for u, _, _ in facs]                      # dual polytope vertices
    facsd = facets(Vd)
    assert all(c == 1 for _, c, _ in facsd)
    pts, ptsd = lattice_points(V, facs), lattice_points(Vd, facsd)

    def side(V_, facs_, pts_):
        sfac = sum(face_interior_count(pts_, facs_, frozenset([i]))
                   for i in range(len(facs_)))
        s2 = 0
        for I, fpair in two_faces(V_, facs_):
            tight = frozenset(i for i in range(len(facs_)) if I <= facs_[i][2])
            lint = face_interior_count(pts_, facs_, tight)
            u1, u2 = facs_[fpair[0]][0], facs_[fpair[1]][0]
            s2 += lint * (vgcd(vsub(u1, u2)) - 1)     # l* of the dual edge
        return sfac, s2

    sfacN, s2N = side(V, facs, pts)
    sfacM, s2M = side(Vd, facsd, ptsd)
    h11 = len(pts) - 5 - sfacN + s2N
    h21 = len(ptsd) - 5 - sfacM + s2M
    if name:
        print(f"{name}: l(D)={len(pts)}, l(D*)={len(ptsd)}, "
              f"(h11, h21)=({h11}, {h21}), chi={2*(h11-h21)}")
    return h11, h21


if __name__ == "__main__":
    e = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
    assert hodge_numbers(e + [(-1, -1, -1, -1)], "quintic") == (1, 101)
    assert hodge_numbers(e + [(-1, -1, -3, -3)], "P(1,1,1,3,3) X9") == (4, 112)
    F = [(x, y, 1, 0) for (x, y) in [(0, 0), (-2, -1), (-2, -2), (-1, -2)]]
    assert hodge_numbers(F + [(0, 0, 0, 1), (1, 1, -1, -1)], "Delta_A") == (5, 101)
    P = [(x, y, 1, 0) for (x, y) in [(0, 0), (-2, -1), (-3, -2), (-2, -3), (-1, -2)]]
    assert hodge_numbers(P + [(0, 0, 0, 1), (1, 1, -1, -1)], "Delta_B") == (9, 81)
    w3v = [(0, 0, 0, 1), (1, 1, -1, 0), (-2, -2, 1, -1)]
    assert hodge_numbers(F + w3v, "Delta_C (single F1 point)") == (4, 108)
    T = [(x, y, 1, 0) for (x, y) in [(0, 0), (-2, -1), (-1, -2)]]
    assert hodge_numbers(T + w3v, "Delta_C' (single 1/3(1,1,1))") == (3, 111)
    VD = [(1, 0, 0, 0), (0, 1, 0, 0), (0, -1, 0, 0), (0, 0, 1, 0),
          (-3, 1, -1, 0), (0, 1, 1, 0), (0, 0, 0, 1), (-3, 2, 1, -1),
          (-2, 2, 2, -1), (1, 0, 1, 1)]
    assert hodge_numbers(VD, "Delta_D (single def-only point)") == (8, 118)
    print("all Hodge-number asserts pass")
