#!/usr/bin/env sage
"""A smooth projective crepant fan for the mirror partner X (Paper 4, the
partial-smoothing section).

Same construction as resolution_fan.sage, applied to Delta_{F_1} itself: one
coherent integral lifting of the 24 boundary lattice points (22 vertices and
the interior points of the hexagon and of the F_1 face), with the two
interior points lowered enough to induce the full star subdivision on both
faces.  Every maximal cone is checked to be unimodular (94 cones), the induced
square diagonals are read off, and an explicit strictly convex integral
support function certifies projectivity.

Run:

    sage mirror_partner_fan.sage
    sage mirror_partner_fan.sage --verbose-cones
"""

from collections import Counter
from itertools import combinations
from pathlib import Path
import sys

PAPER4_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = PAPER4_ROOT.parent.parent
sys.path.insert(0, str(PAPER4_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from batyrev_global import dot, facets
from fixed_example import DELTA_F1_VERTICES
from hodge_numbers import lattice_points

HEXAGON = (9, 14, 15, 20, 21, 22)          # Paper 3, Theorem 6.1
F1_FACE = (4, 18, 19, 22)
# plain Python ints: the preparser's Sage Integers do not compare equal to the
# Fraction constants returned by batyrev_global.facets
INTERIOR = {"H6": tuple(int(x) for x in (-1, 0, -1, 0)),
            "F1": tuple(int(x) for x in (-1, 0, 0, 0))}
RAYS = tuple(DELTA_F1_VERTICES) + (INTERIOR["H6"], INTERIOR["F1"])   # 1..24
PERTURBATIONS = (
    tuple(ZZ(2) ** (21 - index) for index in range(22))
    + tuple(-ZZ(2) ** 30 + ZZ(2) ** (1 - index) for index in range(2))
)

facet_data = facets(DELTA_F1_VERTICES)
all_points = lattice_points(DELTA_F1_VERTICES, facet_data)
assert set(all_points) == set(RAYS).union({(0, 0, 0, 0)}), "lattice points"
assert len(all_points) == 25
assert all(gcd(abs(e) for e in ray) == 1 for ray in RAYS)

# ---------------------------------------------------------------- the fan
maximal_cones = set()
facet_triangulations = []
for normal, constant, _ in facet_data:
    ray_indices = tuple(i for i, ray in enumerate(RAYS, 1)
                        if dot(normal, ray) == constant)
    lifted = tuple(tuple(RAYS[i - 1]) + (PERTURBATIONS[i - 1],)
                   for i in ray_indices)
    P = Polyhedron(vertices=lifted)
    if P.dim() == 3:
        assert len(ray_indices) == 4
        tets = (tuple(sorted(ray_indices)),)
    else:
        assert P.dim() == 4
        tets = tuple(sorted(
            tuple(ray_indices[k] for k, pt in enumerate(lifted)
                  if ineq.eval(vector(ZZ, pt)) == 0)
            for ineq in P.inequality_generator() if ineq.A()[-1] > 0))
    assert all(len(t) == 4 for t in tets)
    assert set(ray_indices) == {r for t in tets for r in t}
    maximal_cones.update(tets)
    facet_triangulations.append(tets)
maximal_cones = tuple(sorted(maximal_cones))

# compatibility of the induced triangulations on facet intersections
for left, right in combinations(range(len(facet_data)), 2):
    ln, lc, _ = facet_data[left]
    rn, rc, _ = facet_data[right]
    common = {i + 1 for i, ray in enumerate(RAYS)
              if dot(ln, ray) == lc and dot(rn, ray) == rc}
    if not common:
        continue
    restr = []
    for fi in (left, right):
        cells = {tuple(sorted(set(t) & common))
                 for t in facet_triangulations[fi] if set(t) & common}
        restr.append({c for c in cells
                      if not any(set(c) < set(o) for o in cells)})
    assert restr[0] == restr[1]

dets = Counter(abs(matrix(ZZ, [RAYS[i - 1] for i in c]).det())
               for c in maximal_cones)
used = {i for c in maximal_cones for i in c}
assert used == set(range(1, 25)), sorted(set(range(1, 25)) - used)
assert dict(dets) == {1: 94}, dets

# ---------------------------------------------------------------- stars
def star_cells(face, center):
    face_rays = set(face) | {center}
    cells = {tuple(sorted(set(c) & face_rays)) for c in maximal_cones
             if len(set(c) & face_rays) >= 3}
    return tuple(sorted(c for c in cells
                        if not any(set(c) < set(o) for o in cells)))

star_H6 = star_cells(HEXAGON, 23)
star_F1 = star_cells(F1_FACE, 24)
assert len(star_H6) == 6 and all(23 in c and len(c) == 3 for c in star_H6)
assert len(star_F1) == 4 and all(24 in c and len(c) == 3 for c in star_F1)

# ---------------------------------------------------------------- diagonals
SQUARES = (
    (1, 2, 7, 12), (1, 2, 9, 15), (1, 3, 8, 13), (1, 3, 9, 14),
    (2, 10, 12, 19), (2, 10, 15, 20), (3, 11, 13, 18), (3, 11, 14, 21),
    (5, 11, 17, 21), (6, 10, 16, 20), (7, 9, 12, 15), (7, 12, 17, 19),
    (7, 14, 17, 21), (8, 9, 13, 14), (8, 13, 16, 18), (8, 15, 16, 20),
    (12, 15, 19, 20), (13, 14, 18, 21), (16, 18, 20, 22), (17, 19, 21, 22),
)
diagonals = []
for sq in SQUARES:
    s = set(sq)
    tris = {tuple(sorted(s & set(c))) for c in maximal_cones
            if len(s & set(c)) == 3}
    tris = {t for t in tris if not any(set(t) < set(o) for o in tris)}
    assert len(tris) == 2, (sq, tris)
    d = tuple(sorted(set.intersection(*(set(t) for t in tris))))
    assert len(d) == 2
    # sanity: the diagonal relation v_a + v_c = v_b + v_d
    a, c = d
    b, e = sorted(s - set(d))
    assert tuple(x + y for x, y in zip(RAYS[a - 1], RAYS[c - 1])) == \
        tuple(x + y for x, y in zip(RAYS[b - 1], RAYS[e - 1]))
    diagonals.append(d)

# ---------------------------------------------------------------- projectivity
required = [ZZ(1)]
pieces = []
for cone in maximal_cones:
    containing = [k for k, (n, c0, _) in enumerate(facet_data)
                  if all(dot(n, RAYS[i - 1]) == c0 for i in cone)]
    assert len(containing) == 1
    normal, constant, _ = facet_data[containing[0]]
    assert QQ(constant) == 1
    M = matrix(QQ, [RAYS[i - 1] for i in cone])
    pert = M.solve_right(vector(QQ, [PERTURBATIONS[i - 1] for i in cone]))
    for i, ray in enumerate(RAYS, 1):
        if i in cone:
            continue
        base_coeff = QQ(1 - dot(normal, ray))
        raw = PERTURBATIONS[i - 1] - vector(QQ, ray) * pert
        assert base_coeff >= 0
        if base_coeff == 0:
            assert raw > 0
        else:
            required.append(floor(-raw / base_coeff) + 1)
    pieces.append((cone, M))
SUPPORT_BASELINE = max(required)
heights = tuple(SUPPORT_BASELINE + p for p in PERTURBATIONS)
margins = []
for cone, M in pieces:
    lin = M.solve_right(vector(QQ, [heights[i - 1] for i in cone]))
    for i, ray in enumerate(RAYS, 1):
        if i in cone:
            continue
        m = heights[i - 1] - vector(QQ, ray) * lin
        assert m > 0
        margins.append(m)

print("mirror_partner_fan: all exact assertions passed")
print("boundary lattice points: 24 nonzero rays plus the origin")
print(f"maximal cones: {len(maximal_cones)}; determinant distribution: "
      f"{dict(dets)}")
print(f"support baseline: {SUPPORT_BASELINE}; minimum convexity margin="
      f"{min(margins)}")
print("stars: hexagon", star_H6)
print("       F1     ", star_F1)
print("square diagonals:")
for k, (sq, d) in enumerate(zip(SQUARES, diagonals), 1):
    print(f"  n_{k:02d} {sq}: diagonal {d}")
if "--verbose-cones" in sys.argv:
    print("maximal cones:")
    for c in maximal_cones:
        print(f"  {c}")
