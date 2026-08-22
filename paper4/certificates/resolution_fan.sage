#!/usr/bin/env sage
"""Fix and verify a smooth projective crepant fan for Paper 4.

The construction uses one coherent lifting of all 30 boundary lattice
points.  The four two-face interior points are lowered enough to induce the
full star subdivision on every dP6 and dP7 face.  An explicit strictly
convex integral support function certifies projectivity.

Run:

    sage resolution_fan.sage
    sage resolution_fan.sage --verbose-cones
"""

from collections import Counter
from itertools import combinations
from pathlib import Path
import ast
import sys


PAPER4_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = PAPER4_ROOT.parent.parent
sys.path.insert(0, str(PAPER4_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from batyrev_global import dot, facets
from fixed_example import (
    DP6_FACES,
    DP7_FACES,
    DP_INTERIOR_POINTS,
    NODE_FACES,
    POLAR_VERTICES,
)
from hodge_numbers import lattice_points


EXTRA_RAYS = (
    DP_INTERIOR_POINTS["D_6,1"],
    DP_INTERIOR_POINTS["D_6,2"],
    DP_INTERIOR_POINTS["D_7,1"],
    DP_INTERIOR_POINTS["D_7,2"],
)
RAYS = POLAR_VERTICES + EXTRA_RAYS

REGULAR_PERTURBATIONS = (
    tuple(ZZ(2) ** (25 - index) for index in range(26))
    + tuple(-ZZ(2) ** 30 + ZZ(2) ** (3 - index) for index in range(4))
)
SUPPORT_BASELINE = 2181038071

EXPECTED_MAXIMAL_CONES = (
    (1, 2, 7, 27),
    (1, 2, 7, 29),
    (1, 2, 13, 27),
    (1, 2, 13, 29),
    (1, 3, 7, 27),
    (1, 3, 7, 29),
    (1, 3, 13, 27),
    (1, 3, 13, 29),
    (2, 4, 8, 27),
    (2, 4, 8, 30),
    (2, 4, 14, 27),
    (2, 4, 14, 30),
    (2, 7, 8, 27),
    (2, 7, 8, 29),
    (2, 8, 29, 30),
    (2, 13, 14, 27),
    (2, 13, 14, 29),
    (2, 14, 29, 30),
    (3, 5, 9, 23),
    (3, 5, 9, 27),
    (3, 5, 13, 23),
    (3, 5, 13, 27),
    (3, 7, 9, 27),
    (3, 7, 9, 29),
    (3, 9, 23, 29),
    (3, 13, 23, 29),
    (4, 6, 10, 15),
    (4, 6, 10, 27),
    (4, 6, 15, 27),
    (4, 8, 10, 27),
    (4, 8, 10, 30),
    (4, 10, 15, 30),
    (4, 14, 15, 27),
    (4, 14, 15, 30),
    (5, 6, 11, 25),
    (5, 6, 11, 27),
    (5, 6, 16, 25),
    (5, 6, 16, 27),
    (5, 9, 11, 23),
    (5, 9, 11, 27),
    (5, 11, 23, 25),
    (5, 13, 16, 23),
    (5, 13, 16, 27),
    (5, 16, 23, 25),
    (6, 10, 12, 15),
    (6, 10, 12, 27),
    (6, 11, 12, 25),
    (6, 11, 12, 27),
    (6, 12, 15, 26),
    (6, 12, 25, 26),
    (6, 15, 16, 26),
    (6, 15, 16, 27),
    (6, 16, 25, 26),
    (7, 8, 17, 28),
    (7, 8, 17, 29),
    (7, 8, 27, 28),
    (7, 9, 17, 28),
    (7, 9, 17, 29),
    (7, 9, 27, 28),
    (8, 10, 18, 28),
    (8, 10, 18, 30),
    (8, 10, 27, 28),
    (8, 17, 18, 28),
    (8, 17, 18, 29),
    (8, 18, 29, 30),
    (9, 11, 19, 23),
    (9, 11, 19, 28),
    (9, 11, 27, 28),
    (9, 17, 19, 28),
    (9, 17, 19, 29),
    (9, 19, 23, 29),
    (10, 12, 15, 24),
    (10, 12, 18, 24),
    (10, 12, 18, 28),
    (10, 12, 27, 28),
    (10, 15, 24, 30),
    (10, 18, 24, 30),
    (11, 12, 19, 25),
    (11, 12, 19, 28),
    (11, 12, 27, 28),
    (11, 19, 23, 25),
    (12, 15, 24, 26),
    (12, 18, 20, 24),
    (12, 18, 20, 28),
    (12, 19, 20, 25),
    (12, 19, 20, 28),
    (12, 20, 24, 26),
    (12, 20, 25, 26),
    (13, 14, 16, 22),
    (13, 14, 16, 27),
    (13, 14, 22, 29),
    (13, 16, 22, 23),
    (13, 21, 22, 23),
    (13, 21, 22, 29),
    (13, 21, 23, 29),
    (14, 15, 16, 22),
    (14, 15, 16, 27),
    (14, 15, 22, 30),
    (14, 22, 29, 30),
    (15, 16, 22, 26),
    (15, 22, 24, 26),
    (15, 22, 24, 30),
    (16, 22, 23, 25),
    (16, 22, 25, 26),
    (17, 18, 19, 21),
    (17, 18, 19, 28),
    (17, 18, 21, 29),
    (17, 19, 21, 29),
    (18, 19, 20, 21),
    (18, 19, 20, 28),
    (18, 20, 21, 24),
    (18, 21, 22, 24),
    (18, 21, 22, 30),
    (18, 21, 29, 30),
    (18, 22, 24, 30),
    (19, 20, 21, 25),
    (19, 21, 23, 25),
    (19, 21, 23, 29),
    (20, 21, 24, 26),
    (20, 21, 25, 26),
    (21, 22, 23, 25),
    (21, 22, 24, 26),
    (21, 22, 25, 26),
    (21, 22, 29, 30),
)

EXPECTED_NODE_DIAGONALS = (
    (2, 7),
    (2, 13),
    (3, 7),
    (4, 8),
    (4, 14),
    (5, 9),
    (5, 13),
    (5, 23),
    (9, 23),
    (6, 10),
    (10, 15),
    (6, 11),
    (6, 16),
    (6, 25),
    (8, 17),
    (9, 17),
    (12, 18),
    (12, 24),
    (12, 19),
    (12, 25),
    (14, 16),
    (16, 23),
    (16, 26),
    (18, 19),
    (20, 24),
    (20, 25),
)

EXPECTED_DIVISORIAL_STARS = {
    "D_6,1": (
        (1, 2, 27), (1, 3, 27), (2, 4, 27),
        (3, 5, 27), (4, 6, 27), (5, 6, 27),
    ),
    "D_6,2": (
        (7, 8, 28), (7, 9, 28), (8, 10, 28),
        (9, 11, 28), (10, 12, 28), (11, 12, 28),
    ),
    "D_7,1": (
        (1, 7, 29), (1, 13, 29), (7, 17, 29),
        (13, 21, 29), (17, 21, 29),
    ),
    "D_7,2": (
        (2, 8, 30), (2, 14, 30), (8, 18, 30),
        (14, 22, 30), (18, 22, 30),
    ),
}


def maximal_cones_and_facets():
    """Take the lower triangulation induced by the common integral lift."""
    facet_data = facets(POLAR_VERTICES)
    maximal_cones = set()
    facet_triangulations = []
    for normal, constant, _vertices in facet_data:
        ray_indices = tuple(
            index for index, ray in enumerate(RAYS, 1)
            if dot(normal, ray) == constant
        )
        lifted_points = tuple(
            tuple(RAYS[index - 1]) + (REGULAR_PERTURBATIONS[index - 1],)
            for index in ray_indices
        )
        lifted_polyhedron = Polyhedron(vertices=lifted_points)

        if lifted_polyhedron.dim() == 3:
            # A facet with exactly four rays needs no subdivision; its lift
            # remains a tetrahedron rather than acquiring lower facets.
            assert len(ray_indices) == 4
            tetrahedra = (tuple(sorted(ray_indices)),)
        else:
            assert lifted_polyhedron.dim() == 4
            tetrahedra = tuple(
                sorted(
                    tuple(
                        ray_indices[local_index]
                        for local_index, point in enumerate(lifted_points)
                        if inequality.eval(vector(ZZ, point)) == 0
                    )
                    for inequality in lifted_polyhedron.inequality_generator()
                    if inequality.A()[-1] > 0
                )
            )
        assert all(len(tetrahedron) == 4 for tetrahedron in tetrahedra)
        assert set(ray_indices) == {ray for cone in tetrahedra for ray in cone}
        maximal_cones.update(tetrahedra)
        facet_triangulations.append(tetrahedra)

    maximal_cones = tuple(sorted(maximal_cones))
    assert maximal_cones == EXPECTED_MAXIMAL_CONES
    assert len(maximal_cones) == 124

    # Directly check that the triangulations induced on every intersection of
    # two original facets agree.
    for left, right in combinations(range(len(facet_data)), 2):
        left_normal, left_constant, _ = facet_data[left]
        right_normal, right_constant, _ = facet_data[right]
        common = {
            index + 1 for index, ray in enumerate(RAYS)
            if dot(left_normal, ray) == left_constant
            and dot(right_normal, ray) == right_constant
        }
        if not common:
            continue

        restrictions = []
        for facet_index in (left, right):
            cells = {
                tuple(sorted(set(tetrahedron).intersection(common)))
                for tetrahedron in facet_triangulations[facet_index]
                if set(tetrahedron).intersection(common)
            }
            maximal_cells = {
                cell for cell in cells
                if not any(set(cell) < set(other) for other in cells)
            }
            restrictions.append(maximal_cells)
        assert restrictions[0] == restrictions[1]

    return facet_data, tuple(facet_triangulations), maximal_cones


def divisorial_stars(maximal_cones):
    """Verify the full star subdivision on every dP6 and dP7 face."""
    face_data = {
        "D_6,1": (DP6_FACES[0], 27),
        "D_6,2": (DP6_FACES[1], 28),
        "D_7,1": (DP7_FACES[0], 29),
        "D_7,2": (DP7_FACES[1], 30),
    }
    stars = {}
    for label, (face, center) in face_data.items():
        face_rays = set(face).union({center})
        cells = {
            tuple(sorted(set(cone).intersection(face_rays)))
            for cone in maximal_cones
            if len(set(cone).intersection(face_rays)) >= 3
        }
        maximal_cells = {
            cell for cell in cells
            if not any(set(cell) < set(other) for other in cells)
        }
        stars[label] = tuple(sorted(maximal_cells))
    assert stars == EXPECTED_DIVISORIAL_STARS
    return stars


def node_diagonals(maximal_cones):
    """Read the chosen diagonal from the two triangles in each square."""
    diagonals = []
    for face in NODE_FACES:
        face_set = set(face)
        triangles = {
            tuple(sorted(face_set.intersection(cone)))
            for cone in maximal_cones
            if len(face_set.intersection(cone)) == 3
        }
        maximal_triangles = {
            triangle for triangle in triangles
            if not any(set(triangle) < set(other) for other in triangles)
        }
        assert len(maximal_triangles) == 2
        diagonal = tuple(
            sorted(set.intersection(*(set(item) for item in maximal_triangles)))
        )
        assert len(diagonal) == 2
        diagonals.append(diagonal)
    diagonals = tuple(diagonals)
    assert diagonals == EXPECTED_NODE_DIAGONALS
    return diagonals


def projectivity_certificate(facet_data, maximal_cones):
    """Verify one integral strictly convex support function on the fan."""
    heights = tuple(
        SUPPORT_BASELINE + entry for entry in REGULAR_PERTURBATIONS
    )
    margins = []
    required_baselines = [ZZ(1)]

    for cone in maximal_cones:
        containing_facets = [
            facet_index
            for facet_index, (normal, constant, _vertices) in enumerate(facet_data)
            if all(dot(normal, RAYS[index - 1]) == constant for index in cone)
        ]
        assert len(containing_facets) == 1
        normal, constant, _vertices = facet_data[containing_facets[0]]
        assert QQ(constant) == 1

        ray_matrix = matrix(QQ, [RAYS[index - 1] for index in cone])
        perturbation_piece = ray_matrix.solve_right(
            vector(
                QQ,
                [REGULAR_PERTURBATIONS[index - 1] for index in cone],
            )
        )
        linear_piece = ray_matrix.solve_right(
            vector(QQ, [heights[index - 1] for index in cone])
        )
        for index, ray in enumerate(RAYS, 1):
            if index in cone:
                continue
            baseline_coefficient = QQ(1 - dot(normal, ray))
            raw_margin = (
                REGULAR_PERTURBATIONS[index - 1]
                - vector(QQ, ray) * perturbation_piece
            )
            assert baseline_coefficient >= 0
            if baseline_coefficient == 0:
                assert raw_margin > 0
            else:
                required_baselines.append(
                    floor(-raw_margin / baseline_coefficient) + 1
                )
            margin = heights[index - 1] - vector(QQ, ray) * linear_piece
            assert margin > 0
            margins.append(margin)

    assert max(required_baselines) == SUPPORT_BASELINE
    assert min(margins) == 1
    return heights, tuple(margins)


facet_data = facets(POLAR_VERTICES)
all_points = lattice_points(POLAR_VERTICES, facet_data)
assert set(all_points) == set(RAYS).union({(0, 0, 0, 0)})
assert len(all_points) == 31
assert all(gcd(abs(entry) for entry in ray) == 1 for ray in RAYS)

facet_data, facet_triangulations, maximal_cones = maximal_cones_and_facets()
determinants = tuple(
    abs(matrix(ZZ, [RAYS[index - 1] for index in cone]).det())
    for cone in maximal_cones
)
assert Counter(determinants) == Counter({1: 124})
assert {index for cone in maximal_cones for index in cone} == set(range(1, 31))

stars = divisorial_stars(maximal_cones)
diagonals = node_diagonals(maximal_cones)
heights, margins = projectivity_certificate(facet_data, maximal_cones)

# The research note prints every maximal cone and every node diagonal.  Parse
# those literal blocks so the human-readable record cannot drift from the
# checker.
note_text = (PAPER4_ROOT / "resolution_data.md").read_text()
cone_block = note_text.split("<!-- MAXIMAL_CONES_BEGIN -->", 1)[1].split(
    "<!-- MAXIMAL_CONES_END -->", 1
)[0]
cone_literal = cone_block.split("MAXIMAL_CONES =", 1)[1].split("```", 1)[0]
assert ast.literal_eval(cone_literal.strip()) == EXPECTED_MAXIMAL_CONES

diagonal_block = note_text.split("<!-- NODE_DIAGONALS_BEGIN -->", 1)[1].split(
    "<!-- NODE_DIAGONALS_END -->", 1
)[0]
diagonal_literal = diagonal_block.split("NODE_DIAGONALS =", 1)[1].split(
    "```", 1
)[0]
assert ast.literal_eval(diagonal_literal.strip()) == EXPECTED_NODE_DIAGONALS

print("resolution_fan: all exact assertions passed")
print("boundary lattice points: 30 nonzero rays plus the origin")
print("maximal cones: 124; determinant distribution: {1: 124}")
print(
    f"support baseline: {SUPPORT_BASELINE}; "
    f"minimum convexity margin={min(margins)}"
)
print("divisorial restrictions: full stars on two dP6 and two dP7 faces")
print("node diagonals:")
for index, diagonal in enumerate(diagonals, 1):
    print(f"  N_{index:02d}: {diagonal}")

if "--verbose-cones" in sys.argv:
    print("maximal cones:")
    for cone in maximal_cones:
        print(f"  {cone}")
