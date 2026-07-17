#!/usr/bin/env python3
"""
The one polytope the per-vertex-count mirror misses (paper 1, Sec. 6.1).

The HuggingFace mirror calabi-yau-data/polytopes-4d stores the Kreuzer-Skarke
classification in files polytopes-4d-NN-vertices.parquet.  The sweep
(src/ks_sweep.py) covered NN = 05..33, which sums to 473,800,775 polytopes --
one short of the classified 473,800,776.  The gap is the single file beyond
that range: polytopes-4d-36-vertices.parquet, one row, the unique reflexive
4-polytope with more than 33 vertices.

This script proves what that polytope is and what it contributes:

  1. Delta_36 (the row of the parquet, hardcoded below) is reflexive with
     36 vertices, 12 facets, 49 lattice points, 13 dual lattice points.
  2. Delta_36 is GL_4(Z)-equivalent to H x H, the product of two copies of
     the reflexive hexagon H = conv{(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)}
     (the dP6 / "1/3(1,1,1)-dual" hexagon).  The explicit unimodular matrix
     is found by matching the facet-normal sets of the two polytopes.
  3. All 72 edges are unit; the 48 two-faces are 36 unit parallelograms
     (conifold faces) and 12 hexagons GL_2(Z)-equivalent to the dP6 hexagon
     (type (S), i=1, 2 smoothing components); every dual edge has lattice
     length 1.  So the generic anticanonical X has 36 ordinary double points
     and 12 dP6-cone points -- every 2-face is of type (S): Delta_36 is NOT
     non-smoothable-by-Corollary-4.3, carries no def-only face, and (having
     hexagon faces) is not in the Batyrev-Kreuzer all-conifold population.
     Hence none of the sweep's headline counts change; only the total rises
     to the full 473,800,776.
  4. Hodge numbers of the MPCP resolution in our convention: (h11, h21) =
     (44, 8), matching the dataset row's (8, 44) under the dataset's swapped
     (mirror) convention, as validated file-by-file by the sweep.

Run:  python3 src/missing_polytope.py     (~2 s; everything asserted)
"""
from fractions import Fraction
from itertools import permutations
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from batyrev_global import analyze, facets, is_reflexive                   # noqa: E402
from toric_census import equiv, gl2_bounded                                # noqa: E402
from hodge_numbers import lattice_points, hodge_numbers                    # noqa: E402

# The single row of polytopes-4d-36-vertices.parquet (vertices column).
DELTA36 = [
    (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (-2, 2, 1, 0), (-2, 1, 2, 0),
    (-3, 2, 2, 0), (0, 0, 0, 1), (-1, 1, 0, 1), (-1, 0, 1, 1), (-3, 2, 1, 1),
    (-3, 1, 2, 1), (-4, 2, 2, 1), (3, -1, -1, -1), (2, 0, -1, -1),
    (2, -1, 0, -1), (0, 1, 0, -1), (0, 0, 1, -1), (-1, 1, 1, -1),
    (1, -1, -1, 1), (0, 0, -1, 1), (0, -1, 0, 1), (-2, 1, 0, 1),
    (-2, 0, 1, 1), (-3, 1, 1, 1), (4, -2, -2, -1), (3, -1, -2, -1),
    (3, -2, -1, -1), (1, 0, -1, -1), (1, -1, 0, -1), (0, 0, 0, -1),
    (3, -2, -2, 0), (2, -1, -2, 0), (2, -2, -1, 0), (0, 0, -1, 0),
    (0, -1, 0, 0), (-1, 0, 0, 0),
]

HEX = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]     # reflexive hexagon
HEXPROD = [(a, b, c, d) for (a, b) in HEX for (c, d) in HEX]   # H x H, 36 vertices

# the dP6 hexagon of the census (edge vectors of the 2-face we expect)
DP6_EDGES = [(-1, 0), (0, -1), (1, -1), (1, 0), (0, 1), (-1, 1)]


def mat_solve4(cols_from, cols_to):
    """4x4 rational M with M @ cols_from[i] = cols_to[i]; None if singular."""
    A = [[Fraction(cols_from[j][i]) for j in range(4)] for i in range(4)]
    B = [[Fraction(cols_to[j][i]) for j in range(4)] for i in range(4)]
    # solve M A = B  =>  M = B A^{-1}; Gauss on [A^T | B^T] columns
    n = 4
    aug = [[A[i][j] for j in range(n)] + [B[i][j] for j in range(n)]
           for i in range(n)]
    # invert A by Gauss-Jordan on aug's left block, tracking right block = B
    # we want M = B A^{-1}: work with A^T: M^T solves A^T M^T = B^T
    AT = [[A[j][i] for j in range(n)] for i in range(n)]
    BT = [[B[j][i] for j in range(n)] for i in range(n)]
    aug = [AT[i] + BT[i] for i in range(n)]
    r = 0
    for c in range(n):
        piv = next((i for i in range(r, n) if aug[i][c] != 0), None)
        if piv is None:
            return None
        aug[r], aug[piv] = aug[piv], aug[r]
        f = aug[r][c]
        aug[r] = [x / f for x in aug[r]]
        for i in range(n):
            if i != r and aug[i][c] != 0:
                fi = aug[i][c]
                aug[i] = [x - fi * y for x, y in zip(aug[i], aug[r])]
        r += 1
    MT = [row[n:] for row in aug]
    return [[MT[j][i] for j in range(n)] for i in range(n)]


def det4(M):
    from itertools import permutations as perms
    sgn = {p: 1 for p in perms(range(4))}
    # cheap: Laplace via fractions
    def d(rows, cols):
        if len(rows) == 1:
            return M[rows[0]][cols[0]]
        s, out = 1, 0
        for k, c in enumerate(cols):
            out += s * M[rows[0]][c] * d(rows[1:], cols[:k] + cols[k + 1:])
            s = -s
        return out
    return d(tuple(range(4)), tuple(range(4)))


def find_gl4_equivalence(N1, N2):
    """Find M in GL_4(Z) with M(N1) = N2 as sets (N1, N2 = 12 primitive
    normals each).  Returns the integer matrix or None."""
    # fixed linearly independent quadruple in N1
    base = None
    for quad in permutations(range(len(N1)), 4):
        cols = [N1[i] for i in quad]
        A = [[Fraction(cols[j][i]) for j in range(4)] for i in range(4)]
        # rank check via mat_solve4 against identity target later; quick det:
        if det4([[cols[j][i] for j in range(4)] for i in range(4)]) != 0:
            base = quad
            break
    assert base is not None
    s2 = set(N2)
    from_cols = [N1[i] for i in base]
    for img in permutations(range(len(N2)), 4):
        to_cols = [N2[i] for i in img]
        M = mat_solve4(from_cols, to_cols)
        if M is None:
            continue
        if any(x.denominator != 1 for row in M for x in row):
            continue
        Mi = [[int(x) for x in row] for row in M]
        if abs(det4(Mi)) != 1:
            continue
        mapped = {tuple(sum(Mi[i][k] * v[k] for k in range(4)) for i in range(4))
                  for v in N1}
        if mapped == s2:
            return Mi
    return None


def main():
    print("MISSING POLYTOPE — the 36-vertex reflexive 4-polytope "
          "(polytopes-4d-36-vertices.parquet, 1 row)")
    print()

    # ---- 1. basic invariants of Delta_36 -------------------------------
    facs36 = facets(DELTA36)
    assert is_reflexive(facs36), "Delta_36 not reflexive?!"
    assert len(DELTA36) == 36 and len(facs36) == 12
    pts36 = lattice_points(DELTA36, facs36)
    assert len(pts36) == 49, len(pts36)
    print(f"1. Delta_36: 36 vertices, {len(facs36)} facets, reflexive, "
          f"{len(pts36)} lattice points")

    # ---- 2. Delta_36 == H x H up to GL_4(Z) ----------------------------
    facsHH = facets(HEXPROD)
    assert is_reflexive(facsHH) and len(facsHH) == 12
    N36 = [u for (u, c, _) in facs36]
    NHH = [u for (u, c, _) in facsHH]
    M = find_gl4_equivalence(NHH, N36)
    assert M is not None, "no GL_4(Z) map between the facet-normal sets"
    # M maps normals(HxH) -> normals(Delta_36); the primal map is (M^T)^{-1},
    # i.e. M^T maps Delta_36 -> H x H.  Verify on the vertex sets directly.
    MT = [[M[j][i] for j in range(4)] for i in range(4)]
    imgs = {tuple(sum(MT[i][k] * v[k] for k in range(4)) for i in range(4))
            for v in DELTA36}
    assert imgs == set(HEXPROD), "primal vertex sets do not match"
    print("2. Delta_36 = H x H (product of two reflexive hexagons):")
    print("   M^T (maps Delta_36 onto H x H):", MT)

    # ---- 3. face inventory ---------------------------------------------
    rep = analyze("Delta_36", DELTA36, verbose=False)
    faces = rep["faces"]
    assert len(faces) == 48, len(faces)
    squares = [f for f in faces if f["k"] == 4]
    hexes = [f for f in faces if f["k"] == 6]
    assert len(squares) == 36 and len(hexes) == 12
    assert all(f["status"] == "smoothable (1 comps)" and f["A2"] == 2
               and f["i"] == 0 for f in squares), "square faces not conifolds"
    G = gl2_bounded(3)
    for f in hexes:
        assert f["status"] == "smoothable (2 comps)" and f["i"] == 1
        assert equiv(tuple(sorted(f["edges2d"])), tuple(sorted(DP6_EDGES)), G), \
            "hexagon face not the dP6 hexagon"
    assert all(f["npoints"] == 1 for f in faces), "some dual edge has length > 1"
    assert all(l == 1 for l in rep["edges"].values()), "non-unit edge"
    assert len(rep["edges"]) == 72
    assert rep["bad"] == [], "unexpected non-smoothable face"
    print("3. 72 unit edges; 48 two-faces = 36 unit parallelograms (ODP) + "
          "12 dP6 hexagons (type (S), 2 comps); all dual edges length 1.")
    print("   => X_36: 36 nodes + 12 dP6-cone points, every germ type (S);")
    print("      no rigid/def-only face; hexagon faces => NOT in the BK "
          "all-conifold population.")

    # ---- 3b. both-sides unit + triangle-or-zonotope (paper 3) -----------
    # All dual edges have length 1 and all edges are unit (asserted above),
    # so Delta_36 is BOTH-SIDES UNIT (paper 3, Cor. 1.3).  Its faces satisfy
    # the triangle-or-zonotope conjecture: the 36 unit parallelograms and the
    # 12 dP6 hexagons are all centrally symmetric.
    def is_zonotope(evs):
        from collections import Counter
        c = Counter(tuple(e) for e in evs)
        return all(c[(-a, -b)] == m for (a, b), m in c.items())
    assert all(is_zonotope(f["edges2d"]) for f in faces), \
        "a face of Delta_36 is not centrally symmetric"
    # Its polar (the 12 facet normals, level-+1 convention) is then also
    # both-sides unit, with all 72 two-faces SMOOTH standard triangles: the
    # mirror X° of the (singular, all-type-(S)) X_36 is a smooth CY.
    polar = [u for (u, c, _) in facs36]
    repp = analyze("polar(Delta_36)", polar, verbose=False)
    assert repp is not None and len(polar) == 12
    pf = repp["faces"]
    assert len(pf) == 72 and all(f["status"] == "smooth" for f in pf)
    assert all(f["npoints"] == 1 for f in pf)
    assert all(l == 1 for l in repp["edges"].values())
    print("3b. Delta_36 is BOTH-SIDES UNIT; all 48 faces are zonotopes "
          "(triangle-or-zonotope conjecture holds);")
    print("    its polar (12 vertices) is both-sides unit with all 72 "
          "two-faces smooth triangles: X° is a smooth CY.")

    # ---- 4. Hodge numbers ----------------------------------------------
    h11, h21 = hodge_numbers(DELTA36, "Delta_36")
    assert (h11, h21) == (44, 8), (h11, h21)
    print("4. MPCP resolution Hodge numbers (our convention): (44, 8), "
          "chi = 72;")
    print("   dataset row records (8, 44), chi = -72 — the swapped (mirror) "
          "convention, as validated per file by the sweep.")
    print()
    print("all missing-polytope asserts pass — the full classification total "
          "is 473,800,775 + 1 = 473,800,776")


if __name__ == "__main__":
    main()
