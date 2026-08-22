#!/usr/bin/env sage
"""Compute the Cox grading, irrelevant ideal, and singular-chart overlaps.

The fan is the face fan of the fixed 26-vertex polytope Delta^vee.  The
script constructs the divisor-class quotient integrally, rather than using
a rational nullspace, and records every pairwise intersection among the 30
singular face cones.

Run:

    sage cox_data.sage
    sage cox_data.sage --verbose-overlaps
"""

from collections import Counter
from itertools import combinations
from pathlib import Path
import sys


PAPER4_ROOT = globals().get(
    "PAPER4_ROOT_OVERRIDE", Path(__file__).resolve().parent
)
PROJECT_ROOT = PAPER4_ROOT.parent.parent
COX_DATA_QUIET = globals().get("COX_DATA_QUIET", False)
sys.path.insert(0, str(PAPER4_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from batyrev_global import facets
from fixed_example import DP6_FACES, DP7_FACES, NODE_FACES, POLAR_VERTICES


# These four rays span a smooth maximal cone.  Their 4 by 4 ray matrix is
# unimodular, so the classes of the other 22 invariant divisors form a basis
# of the torsion-free class group.
PIVOT_RAYS = (15, 22, 24, 26)

EXPECTED_SMITH_DIAGONAL = (1, 1, 1, 1)
EXPECTED_ANTICANONICAL_DEGREE = (
    3, 2, 3, 1, 2, 1, 3, 2, 3, 1, 2,
    1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1,
)

EXPECTED_FACETS = (
    (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
    (1, 2, 3, 4, 5, 6, 13, 14, 15, 16),
    (1, 2, 7, 8, 13, 14, 17, 18, 21, 22),
    (1, 3, 7, 9, 13, 17, 19, 21, 23),
    (2, 4, 8, 10, 14, 15, 18, 22, 24),
    (3, 5, 9, 11, 19, 23, 25),
    (3, 5, 13, 16, 23, 25),
    (4, 6, 10, 12, 15, 24, 26),
    (5, 6, 11, 12, 25, 26),
    (5, 6, 15, 16, 25, 26),
    (7, 8, 9, 10, 11, 12, 17, 18, 19, 20),
    (10, 12, 18, 20, 24, 26),
    (11, 12, 19, 20, 25, 26),
    (13, 14, 15, 16, 22),
    (13, 16, 21, 22, 23, 25),
    (15, 16, 22, 25, 26),
    (15, 22, 24, 26),
    (17, 18, 19, 20, 21),
    (18, 20, 21, 22, 24, 26),
    (19, 20, 21, 25, 26),
    (19, 21, 23, 25),
    (21, 22, 25, 26),
)

EXPECTED_OVERLAP_COUNTS = {0: 261, 1: 78, 2: 96}


def cox_grading():
    """Return a split integral quotient Z^26 -> Cl(X) = Z^22."""
    ray_matrix = matrix(ZZ, POLAR_VERTICES)
    assert ray_matrix.dimensions() == (26, 4)
    assert ray_matrix.rank() == 4

    smith, left_transform, right_transform = ray_matrix.smith_form()
    assert left_transform * ray_matrix * right_transform == smith
    smith_diagonal = tuple(smith[index, index] for index in range(4))
    assert smith_diagonal == EXPECTED_SMITH_DIAGONAL

    pivot_indices = tuple(index - 1 for index in PIVOT_RAYS)
    free_indices = tuple(
        index for index in range(26) if index not in pivot_indices
    )
    pivot_matrix = ray_matrix.matrix_from_rows(pivot_indices)
    assert pivot_matrix.det() == -1

    # If the nonpivot divisor classes are the standard basis, the relation
    # Q*A=0 forces these four pivot-degree columns.
    pivot_degrees_rational = (
        -ray_matrix.matrix_from_rows(free_indices) * pivot_matrix.inverse()
    )
    assert all(entry in ZZ for entry in pivot_degrees_rational.list())
    pivot_degrees = pivot_degrees_rational.change_ring(ZZ)

    grading = zero_matrix(ZZ, 22, 26)
    for row, column in enumerate(free_indices):
        grading[row, column] = 1
    for row in range(22):
        for local_column, global_column in enumerate(pivot_indices):
            grading[row, global_column] = pivot_degrees[row, local_column]

    assert grading * ray_matrix == 0
    assert grading.rank() == 22
    assert grading.matrix_from_columns(free_indices) == identity_matrix(ZZ, 22)

    anticanonical_degree = tuple(grading * vector(ZZ, [1] * 26))
    assert anticanonical_degree == EXPECTED_ANTICANONICAL_DEGREE

    return {
        "ray_matrix": ray_matrix,
        "smith_diagonal": smith_diagonal,
        "pivot_indices": pivot_indices,
        "free_indices": free_indices,
        "pivot_degrees": pivot_degrees,
        "grading": grading,
        "anticanonical_degree": anticanonical_degree,
    }


def fan_data():
    """Return the maximal cones and irrelevant-ideal generator supports."""
    facet_data = facets(POLAR_VERTICES)
    facet_ray_sets = tuple(
        sorted(
            tuple(sorted(index + 1 for index in face))
            for _normal, _constant, face in facet_data
        )
    )
    assert facet_ray_sets == EXPECTED_FACETS
    all_rays = frozenset(range(1, 27))
    irrelevant_supports = tuple(
        tuple(sorted(all_rays.difference(face))) for face in facet_ray_sets
    )
    assert len(irrelevant_supports) == 22
    assert all(irrelevant_supports)
    assert len(set(irrelevant_supports)) == 22
    return {
        "facets": facet_ray_sets,
        "irrelevant_supports": irrelevant_supports,
    }


def singular_overlap_data():
    """Record all intersections U_F cap U_G = U_{cone(F cap G)}."""
    singular_faces = tuple(
        [(f"N_{index}", face) for index, face in enumerate(NODE_FACES, 1)]
        + [
            ("D_6,1", DP6_FACES[0]),
            ("D_6,2", DP6_FACES[1]),
            ("D_7,1", DP7_FACES[0]),
            ("D_7,2", DP7_FACES[1]),
        ]
    )
    assert len(singular_faces) == 30

    records = []
    for (left_label, left_face), (right_label, right_face) in combinations(
        singular_faces, 2
    ):
        common_rays = tuple(sorted(set(left_face).intersection(right_face)))
        records.append((left_label, right_label, common_rays))

    assert len(records) == binomial(30, 2) == 435
    counts = Counter(len(common) for _left, _right, common in records)
    assert dict(counts) == EXPECTED_OVERLAP_COUNTS
    assert max(counts) == 2

    # The two node-only coloops touch both dP6 cones along a ray, but no dP7
    # cone along a positive-dimensional fan cone.
    by_pair = {
        frozenset((left, right)): common for left, right, common in records
    }
    assert by_pair[frozenset(("N_9", "D_6,1"))] == (3,)
    assert by_pair[frozenset(("N_9", "D_6,2"))] == (9,)
    assert by_pair[frozenset(("N_11", "D_6,1"))] == (4,)
    assert by_pair[frozenset(("N_11", "D_6,2"))] == (10,)
    for node in ("N_9", "N_11"):
        for divisor in ("D_7,1", "D_7,2"):
            assert by_pair[frozenset((node, divisor))] == ()

    return {
        "faces": singular_faces,
        "records": tuple(records),
        "counts": dict(counts),
    }


cox = cox_grading()
fan = fan_data()
overlaps = singular_overlap_data()

if not COX_DATA_QUIET:
    print("cox_data: all exact assertions passed")
    print(f"class group: Z^22; Smith diagonal={cox['smith_diagonal']}")
    print(f"pivot rays: {PIVOT_RAYS}; free divisor rays="
          f"{tuple(index + 1 for index in cox['free_indices'])}")
    for column, ray in enumerate(PIVOT_RAYS):
        degree = tuple(cox["pivot_degrees"].column(column))
        print(f"deg(X_{ray})={degree}")
    print(f"anticanonical degree={cox['anticanonical_degree']}")
    print("irrelevant supports:")
    for index, support in enumerate(fan["irrelevant_supports"], 1):
        print(f"  B_{index:02d}: {support}")
    print(f"overlap counts by number of common rays: {overlaps['counts']}")

    print("divisorial-chart incidences:")
    for divisor in ("D_6,1", "D_6,2", "D_7,1", "D_7,2"):
        adjacent = []
        for left, right, common in overlaps["records"]:
            if common and divisor in (left, right):
                other = right if left == divisor else left
                adjacent.append((other, common))
        print(f"  {divisor}: {tuple(adjacent)}")

    if "--verbose-overlaps" in sys.argv:
        print("all non-torus-only overlaps:")
        for left, right, common in overlaps["records"]:
            if common:
                print(f"  {left}, {right}: {common}")
