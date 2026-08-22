#!/usr/bin/env sage
"""Exact divisor restrictions and exceptional-curve functionals for Paper 4.

The fixed MPCP fan has 30 rays.  This script computes its split Picard
grading, the restriction of every invariant divisor to each of the four
exceptional del Pezzo surfaces, and the divisor-intersection functional of
each of the 26 exceptional node curves.  The final rank calculation is a
numerical input for a possible mixed Friedman--Gross comparison; it is not
itself such a comparison theorem.

Run:  sage restriction_data.sage
"""

from pathlib import Path
import ast
import sys


PAPER4_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = PAPER4_ROOT.parent.parent
sys.path.insert(0, str(PAPER4_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from batyrev_global import facets, int_kernel, solve_int_coords, two_faces, vsub
from fixed_example import (
    DELTA_F1_VERTICES,
    DP6_FACES,
    DP7_FACES,
    DP_INTERIOR_POINTS,
    LOCAL_EDGE_MODELS,
    NODE_FACES,
    POLAR_VERTICES,
)
from paper3_node_relations import diagonal_relation


EXTRA_RAYS = (
    DP_INTERIOR_POINTS["D_6,1"],
    DP_INTERIOR_POINTS["D_6,2"],
    DP_INTERIOR_POINTS["D_7,1"],
    DP_INTERIOR_POINTS["D_7,2"],
)
RAYS = POLAR_VERTICES + EXTRA_RAYS

# These four rays span a maximal cone in the fixed fan and form a unimodular
# basis of the lattice.  The remaining invariant divisors are our Picard
# basis, in the displayed order.
PIVOT_RAYS = (15, 22, 24, 26)
FREE_DIVISOR_RAYS = tuple(
    index for index in range(1, 31) if index not in PIVOT_RAYS
)

EXPECTED_ANTICANONICAL_CLASS = (
    3, 2, 3, 1, 2, 1, 3, 2, 3, 1, 2, 1, 2,
    1, 1, 2, 1, 2, 1, 1, 2, 1, 2, 2, 2, 1,
)

NODE_DIAGONALS = (
    (2, 7), (2, 13), (3, 7), (4, 8), (4, 14), (5, 9),
    (5, 13), (5, 23), (9, 23), (6, 10), (10, 15), (6, 11),
    (6, 16), (6, 25), (8, 17), (9, 17), (12, 18), (12, 24),
    (12, 19), (12, 25), (14, 16), (16, 23), (16, 26),
    (18, 19), (20, 24), (20, 25),
)

SURFACE_MODELS = {
    "D_6,1": {
        "face": DP6_FACES[0],
        "center": 27,
        "cyclic": (1, 3, 5, 6, 4, 2),
        "expected_grading": (
            (-1, 1, 1, 0, 0, 0),
            (1, -1, 0, 1, 0, 0),
            (0, 1, 0, 0, 1, 0),
            (1, 0, 0, 0, 0, 1),
        ),
        "expected_canonical": (-1, -1, -2, -2),
        "expected_intersection": (
            (-1, 0, 1, 0),
            (0, -1, 0, 1),
            (1, 0, -1, 1),
            (0, 1, 1, -1),
        ),
        "expected_self_intersections": (-1, -1, -1, -1, -1, -1),
        "expected_degree": 6,
    },
    "D_6,2": {
        "face": DP6_FACES[1],
        "center": 28,
        "cyclic": (7, 9, 11, 12, 10, 8),
        "expected_grading": (
            (-1, 1, 1, 0, 0, 0),
            (1, -1, 0, 1, 0, 0),
            (0, 1, 0, 0, 1, 0),
            (1, 0, 0, 0, 0, 1),
        ),
        "expected_canonical": (-1, -1, -2, -2),
        "expected_intersection": (
            (-1, 0, 1, 0),
            (0, -1, 0, 1),
            (1, 0, -1, 1),
            (0, 1, 1, -1),
        ),
        "expected_self_intersections": (-1, -1, -1, -1, -1, -1),
        "expected_degree": 6,
    },
    "D_7,1": {
        "face": DP7_FACES[0],
        "center": 29,
        "cyclic": (1, 13, 21, 17, 7),
        "expected_grading": (
            (-1, 1, 1, 0, 0),
            (1, -1, 0, 1, 0),
            (1, 0, 0, 0, 1),
        ),
        "expected_canonical": (-1, -1, -2),
        "expected_intersection": (
            (0, 0, 1),
            (0, -1, 1),
            (1, 1, 0),
        ),
        "expected_self_intersections": (-1, 0, 0, -1, -1),
        "expected_degree": 7,
    },
    "D_7,2": {
        "face": DP7_FACES[1],
        "center": 30,
        "cyclic": (2, 14, 22, 18, 8),
        "expected_grading": (
            (-1, 1, 1, 0, 0),
            (1, -1, 0, 1, 0),
            (0, 1, 0, 0, 1),
        ),
        "expected_canonical": (-1, -1, -2),
        "expected_intersection": (
            (-1, 0, 1),
            (0, 0, 1),
            (1, 1, 0),
        ),
        "expected_self_intersections": (-1, -1, 0, 0, -1),
        "expected_degree": 7,
    },
}


def smith_nonzero_diagonal(input_matrix):
    """Return the positive nonzero Smith diagonal, checking transforms."""
    smith, left, right = input_matrix.smith_form()
    assert left * input_matrix * right == smith
    return tuple(
        abs(smith[index, index])
        for index in range(min(smith.dimensions()))
        if smith[index, index]
    )


def ambient_grading():
    """Return the split quotient Div_T -> Pic of the resolved ambient fan."""
    ray_matrix = matrix(ZZ, RAYS)
    assert ray_matrix.dimensions() == (30, 4)
    assert ray_matrix.rank() == 4
    smith, left, right = ray_matrix.smith_form()
    assert left * ray_matrix * right == smith
    assert tuple(smith[index, index] for index in range(4)) == (1, 1, 1, 1)

    pivot_indices = tuple(index - 1 for index in PIVOT_RAYS)
    free_indices = tuple(index - 1 for index in FREE_DIVISOR_RAYS)
    pivot_matrix = ray_matrix.matrix_from_rows(pivot_indices)
    assert pivot_matrix.det() == -1
    pivot_degrees = (
        -ray_matrix.matrix_from_rows(free_indices) * pivot_matrix.inverse()
    )
    assert all(entry in ZZ for entry in pivot_degrees.list())
    pivot_degrees = pivot_degrees.change_ring(ZZ)

    grading = zero_matrix(ZZ, 26, 30)
    for row, column in enumerate(free_indices):
        grading[row, column] = 1
    for row in range(26):
        for local_column, global_column in enumerate(pivot_indices):
            grading[row, global_column] = pivot_degrees[row, local_column]

    assert grading * ray_matrix == 0
    assert grading.matrix_from_columns(free_indices) == identity_matrix(ZZ, 26)
    anticanonical = tuple(grading * vector(ZZ, [1] * 30))
    assert anticanonical == EXPECTED_ANTICANONICAL_CLASS
    return ray_matrix, grading, free_indices, anticanonical


def canonical_face_coordinates(label, model):
    """Map labeled face rays and its center to Altmann's polygon lattice."""
    facet_data = facets(POLAR_VERTICES)
    face_data = two_faces(POLAR_VERTICES, facet_data)
    face = frozenset(index - 1 for index in model["face"])
    containing = next(
        containing for candidate, containing in face_data if candidate == face
    )
    difference_basis = int_kernel(
        [list(facet_data[index][0]) for index in containing]
    )
    base_vertex = POLAR_VERTICES[min(face)]
    linear_map = matrix(ZZ, LOCAL_EDGE_MODELS[label]["matrix"])
    translation = vector(ZZ, LOCAL_EDGE_MODELS[label]["translation"])

    def coordinate(point):
        local = vector(
            ZZ,
            solve_int_coords(difference_basis, vsub(point, base_vertex)),
        )
        return tuple(linear_map * local + translation)

    coordinates = {
        index: coordinate(RAYS[index - 1])
        for index in tuple(model["face"]) + (model["center"],)
    }
    center = vector(ZZ, coordinates[model["center"]])
    surface_rays = {
        index: vector(ZZ, coordinates[index]) - center
        for index in model["face"]
    }
    assert all(
        matrix(
            ZZ,
            [
                surface_rays[model["cyclic"][index]],
                surface_rays[
                    model["cyclic"][(index + 1) % len(model["cyclic"])]
                ],
            ],
        ).det() == 1
        for index in range(len(model["cyclic"]))
    )
    return coordinates, surface_rays


def surface_record(label, model, ambient_data):
    """Compute Pic(E), its intersection form, and Pic(Y) -> Pic(E)."""
    ray_matrix, ambient_pic_grading, free_indices, _anticanonical = ambient_data
    coordinates, surface_rays = canonical_face_coordinates(label, model)
    boundary = model["face"]
    boundary_matrix = matrix(ZZ, [surface_rays[index] for index in boundary])
    assert boundary_matrix.rank() == 2

    pivot_matrix = boundary_matrix.matrix_from_rows((0, 1))
    assert abs(pivot_matrix.det()) == 1
    free_positions = tuple(range(2, len(boundary)))
    pivot_degrees = (
        -boundary_matrix.matrix_from_rows(free_positions)
        * pivot_matrix.inverse()
    ).change_ring(ZZ)
    picard_rank = len(boundary) - 2
    surface_grading = zero_matrix(ZZ, picard_rank, len(boundary))
    for row, column in enumerate(free_positions):
        surface_grading[row, column] = 1
    for row in range(picard_rank):
        for column in range(2):
            surface_grading[row, column] = pivot_degrees[row, column]
    assert surface_grading * boundary_matrix == 0
    assert tuple(tuple(row) for row in surface_grading.rows()) == model[
        "expected_grading"
    ]

    # Toric boundary intersections in the verified counterclockwise order.
    boundary_position = {ray: position for position, ray in enumerate(boundary)}
    intersection = zero_matrix(ZZ, len(boundary), len(boundary))
    self_intersections = []
    cyclic = model["cyclic"]
    for cyclic_position, ray in enumerate(cyclic):
        previous_ray = cyclic[cyclic_position - 1]
        next_ray = cyclic[(cyclic_position + 1) % len(cyclic)]
        previous_vector = surface_rays[previous_ray]
        ray_vector = surface_rays[ray]
        next_vector = surface_rays[next_ray]
        self_intersection = -matrix(ZZ, [previous_vector, next_vector]).det()
        assert previous_vector + next_vector + self_intersection * ray_vector == 0
        self_intersections.append(self_intersection)
        row = boundary_position[ray]
        intersection[row, row] = self_intersection
        intersection[row, boundary_position[previous_ray]] = 1
        intersection[row, boundary_position[next_ray]] = 1
    assert tuple(self_intersections) == model["expected_self_intersections"]
    assert intersection * boundary_matrix == 0

    picard_intersection = intersection.matrix_from_rows_and_columns(
        free_positions, free_positions
    )
    assert intersection == (
        surface_grading.transpose() * picard_intersection * surface_grading
    )
    assert tuple(tuple(row) for row in picard_intersection.rows()) == model[
        "expected_intersection"
    ]

    canonical = -surface_grading * vector(ZZ, [1] * len(boundary))
    assert tuple(canonical) == model["expected_canonical"]
    degree = canonical * picard_intersection * canonical
    assert degree == model["expected_degree"]

    # At the generic point of the singular orbit, only boundary divisors meet
    # the exceptional fiber.  Adjunction gives D_center|_E = K_E.
    divisor_restriction = zero_matrix(ZZ, picard_rank, 30)
    for local_column, ray in enumerate(boundary):
        divisor_restriction.set_column(
            ray - 1, surface_grading.column(local_column)
        )
    divisor_restriction.set_column(model["center"] - 1, canonical)
    assert divisor_restriction * ray_matrix == 0

    picard_restriction = divisor_restriction.matrix_from_columns(free_indices)
    assert divisor_restriction == picard_restriction * ambient_pic_grading
    assert picard_restriction.rank() == picard_rank
    assert smith_nonzero_diagonal(picard_restriction) == (1,) * picard_rank

    return {
        "label": label,
        "boundary": boundary,
        "center": model["center"],
        "cyclic": cyclic,
        "coordinates": coordinates,
        "surface_rays": surface_rays,
        "surface_grading": surface_grading,
        "canonical": canonical,
        "intersection": picard_intersection,
        "self_intersections": tuple(self_intersections),
        "degree": degree,
        "divisor_restriction": divisor_restriction,
        "picard_restriction": picard_restriction,
    }


def node_curve_matrix(ambient_data):
    """Return the 26 chosen exceptional-curve functionals on Pic(Y)."""
    ray_matrix, grading, free_indices, _anticanonical = ambient_data
    rows = []
    divisor_rows = []
    for face, diagonal in zip(NODE_FACES, NODE_DIAGONALS):
        derived_face, raw_row = diagonal_relation(
            POLAR_VERTICES, frozenset(index - 1 for index in face)
        )
        assert tuple(index + 1 for index in derived_face) == face
        assert set(diagonal).issubset(face)
        diagonal_coefficients = {raw_row[index - 1] for index in diagonal}
        assert len(diagonal_coefficients) == 1
        sign = 1 if diagonal_coefficients == {-1} else -1
        divisor_row = tuple(sign * entry for entry in raw_row) + (0, 0, 0, 0)
        assert all(divisor_row[index - 1] == -1 for index in diagonal)
        assert vector(ZZ, divisor_row) * ray_matrix == 0
        picard_row = tuple(divisor_row[index] for index in free_indices)
        assert vector(ZZ, picard_row) * grading == vector(ZZ, divisor_row)
        divisor_rows.append(divisor_row)
        rows.append(picard_row)
    curve_matrix = matrix(ZZ, rows)
    assert curve_matrix.dimensions() == (26, 26)
    assert curve_matrix.rank() == 18
    return curve_matrix, tuple(divisor_rows)


# The diagonal choices are parsed from the fan note so the two exact scripts
# cannot silently drift apart.
resolution_note = (PAPER4_ROOT / "resolution_data.md").read_text()
diagonal_block = resolution_note.split("<!-- NODE_DIAGONALS_BEGIN -->", 1)[1].split(
    "<!-- NODE_DIAGONALS_END -->", 1
)[0]
diagonal_literal = diagonal_block.split("NODE_DIAGONALS =", 1)[1].split("```", 1)[0]
assert ast.literal_eval(diagonal_literal.strip()) == NODE_DIAGONALS

ambient = ambient_grading()
surfaces = tuple(
    surface_record(label, model, ambient)
    for label, model in SURFACE_MODELS.items()
)
surface_map = matrix(
    ZZ,
    [
        list(row)
        for record in surfaces
        for row in record["picard_restriction"].rows()
    ],
)
assert surface_map.dimensions() == (14, 26)
assert surface_map.rank() == 14
assert smith_nonzero_diagonal(surface_map) == (1,) * 14

curve_map, node_divisor_rows = node_curve_matrix(ambient)
combined_map = curve_map.stack(surface_map)
assert combined_map.dimensions() == (40, 26)
assert combined_map.rank() == 25
assert smith_nonzero_diagonal(combined_map) == (1,) * 25

anticanonical = vector(ZZ, ambient[3])
assert gcd(abs(entry) for entry in anticanonical) == 1
assert combined_map * anticanonical == 0
kernel_basis = combined_map.right_kernel_matrix()
assert kernel_basis.dimensions() == (1, 26)
kernel_generator = vector(ZZ, kernel_basis.row(0))
if kernel_generator[0] < 0:
    kernel_generator = -kernel_generator
assert kernel_generator == anticanonical

# The anticanonical moment polytope of the resolved toric fourfold is the
# original reflexive polytope.  Its normalized volume is the top
# anticanonical intersection.  This is used to prove that the one ambient
# class invisible on the exceptional strata does not restrict trivially to
# the anticanonical Calabi--Yau hypersurface.
anticanonical_top_intersection = (
    factorial(4) * Polyhedron(vertices=DELTA_F1_VERTICES).volume()
)
assert anticanonical_top_intersection == 94

def report():
    """Print the compact human-readable verification summary."""
    print("restriction_data: all exact assertions passed")
    print(f"ambient Picard basis rays: {FREE_DIVISOR_RAYS}")
    print(f"ambient anticanonical class: {tuple(anticanonical)}")
    print(f"ambient anticanonical fourth power: {anticanonical_top_intersection}")
    for record in surfaces:
        print(
            f"{record['label']}: "
            f"Picard rank={record['picard_restriction'].nrows()}, "
            f"cyclic rays={record['cyclic']}, "
            f"boundary self-intersections={record['self_intersections']}, "
            f"K^2={record['degree']}"
        )
        print(f"  K={tuple(record['canonical'])}")
        print(
            "  boundary grading="
            f"{tuple(tuple(row) for row in record['surface_grading'].rows())}"
        )
        print(
            "  intersection form="
            f"{tuple(tuple(row) for row in record['intersection'].rows())}"
        )
    print(
        "surface restriction map: 14 x 26, rank 14, "
        "nonzero Smith invariants all 1"
    )
    print("node curve map: 26 x 26, rank 18")
    print(
        "combined numerical map: 40 x 26, rank 25, "
        "nonzero Smith invariants all 1"
    )
    print("combined kernel: primitive anticanonical class")


if not globals().get("RESTRICTION_DATA_QUIET", False):
    report()
