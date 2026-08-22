#!/usr/bin/env sage
"""Fix the generic anticanonical equation used by Paper 4.

The script enumerates the 25 lattice monomials of Delta, normalizes the
coefficients of 1,x1,x2,x3,x4 by scalar and dense-torus actions, and builds
the remaining 20-parameter Laurent polynomial over a rational function
field.  It also homogenizes every monomial in the 26 Cox variables and uses
exact facewise specializations to prove that every one of the 189 face
nondegeneracy opens is nonempty.  Hence their finite intersection contains
the generic point of the irreducible coefficient chart.

Run:

    sage global_section.sage
    sage global_section.sage --verify-common-witness

The optional common-witness gate is substantially more expensive than the
standard facewise proof: it checks that the distinct-prime coefficient
specialization is nondegenerate on the full-dimensional face as well as on
every proper face.
"""

from itertools import product
from pathlib import Path
import sys


PAPER4_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = PAPER4_ROOT.parent.parent
GLOBAL_SECTION_QUIET = globals().get("GLOBAL_SECTION_QUIET", False)
GLOBAL_SECTION_VERIFY_FACES = globals().get(
    "GLOBAL_SECTION_VERIFY_FACES", True
)
GLOBAL_SECTION_VERIFY_COMMON_WITNESS = globals().get(
    "GLOBAL_SECTION_VERIFY_COMMON_WITNESS",
    "--verify-common-witness" in sys.argv[1:],
)
sys.path.insert(0, str(PAPER4_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from batyrev_global import facets
from fixed_example import (
    DELTA_F1_VERTICES,
    DP6_FACES,
    DP7_FACES,
    NODE_FACES,
    POLAR_VERTICES,
)


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def lattice_points_of_delta():
    """Enumerate Delta from its 26 polar inequalities."""
    bounds = tuple(
        (
            min(vertex[index] for vertex in DELTA_F1_VERTICES),
            max(vertex[index] for vertex in DELTA_F1_VERTICES),
        )
        for index in range(4)
    )
    points = tuple(
        point
        for point in product(*(range(lower, upper + 1) for lower, upper in bounds))
        if all(dot(point, ray) >= -1 for ray in POLAR_VERTICES)
    )
    assert len(points) == 25
    assert set(DELTA_F1_VERTICES).issubset(points)
    return points


def nonempty_face_vertex_sets():
    """Return every nonempty face as a set of Delta vertex indices."""
    facet_records = facets(DELTA_F1_VERTICES)
    facet_sets = tuple(
        frozenset(indices) for _normal, _constant, indices in facet_records
    )
    full_face = frozenset(range(len(DELTA_F1_VERTICES)))
    faces = {full_face}
    frontier = [full_face]
    while frontier:
        face = frontier.pop()
        for facet in facet_sets:
            intersection = face.intersection(facet)
            if intersection and intersection not in faces:
                faces.add(intersection)
                frontier.append(intersection)
    ordered = tuple(sorted(faces, key=lambda face: (len(face), tuple(face))))
    dimensions = {}
    for face in ordered:
        vertices = [vector(QQ, DELTA_F1_VERTICES[index]) for index in face]
        base = vertices[0]
        dimension = matrix(QQ, [vertex - base for vertex in vertices[1:]]).rank()
        dimensions[face] = dimension
    dimension_counts = {
        dimension: sum(value == dimension for value in dimensions.values())
        for dimension in range(5)
    }
    assert dimension_counts == {0: 22, 1: 68, 2: 72, 3: 26, 4: 1}
    assert len(ordered) == 189
    return facet_records, ordered, dimensions


DELTA_LATTICE_POINTS = lattice_points_of_delta()
GAUGE_POINTS = (
    (0, 0, 0, 0),
    (1, 0, 0, 0),
    (0, 1, 0, 0),
    (0, 0, 1, 0),
    (0, 0, 0, 1),
)
assert all(point in DELTA_LATTICE_POINTS for point in GAUGE_POINTS)
PARAMETER_POINTS = tuple(
    point for point in DELTA_LATTICE_POINTS if point not in GAUGE_POINTS
)
assert len(PARAMETER_POINTS) == 20

parameter_names = tuple(f"a_{index:02d}" for index in range(1, 21))
coefficient_ring = PolynomialRing(QQ, names=parameter_names)
coefficient_field = coefficient_ring.fraction_field()
laurent_ring = LaurentPolynomialRing(
    coefficient_field, names=("x1", "x2", "x3", "x4")
)
x = laurent_ring.gens()


def laurent_monomial(ring_generators, exponent):
    return prod(
        ring_generators[index] ** exponent[index] for index in range(4)
    )


generic_coefficients = {
    point: coefficient_field(1) for point in GAUGE_POINTS
}
generic_coefficients.update(
    {
        point: coefficient_field.gen(index)
        for index, point in enumerate(PARAMETER_POINTS)
    }
)
generic_laurent_polynomial = sum(
    generic_coefficients[point] * laurent_monomial(x, point)
    for point in DELTA_LATTICE_POINTS
)
assert len(generic_laurent_polynomial.dict()) == 25

coefficient_labels = {point: "1" for point in GAUGE_POINTS}
coefficient_labels.update(dict(zip(PARAMETER_POINTS, parameter_names)))

# A character x^m homogenizes to prod_j X_j^(<m,q_j>+1).
cox_exponents = {
    point: tuple(dot(point, ray) + 1 for ray in POLAR_VERTICES)
    for point in DELTA_LATTICE_POINTS
}
assert all(min(exponents) >= 0 for exponents in cox_exponents.values())

PAPER4_ROOT_OVERRIDE = PAPER4_ROOT
COX_DATA_QUIET = True
load(str(PAPER4_ROOT / "cox_data.sage"))
for exponents in cox_exponents.values():
    assert tuple(cox["grading"] * vector(ZZ, exponents)) == tuple(
        cox["anticanonical_degree"]
    )


singular_faces = tuple(
    [(f"N_{index}", face) for index, face in enumerate(NODE_FACES, 1)]
    + [
        ("D_6,1", DP6_FACES[0]),
        ("D_6,2", DP6_FACES[1]),
        ("D_7,1", DP7_FACES[0]),
        ("D_7,2", DP7_FACES[1]),
    ]
)
singular_orbit_records = []
for label, face in singular_faces:
    dual_edge = tuple(
        point
        for point in DELTA_LATTICE_POINTS
        if all(dot(point, POLAR_VERTICES[index - 1]) == -1 for index in face)
    )
    assert len(dual_edge) == 2
    left, right = dual_edge
    direction = tuple(right[index] - left[index] for index in range(4))
    assert gcd(abs(entry) for entry in direction) == 1
    left_coefficient = generic_coefficients[left]
    right_coefficient = generic_coefficients[right]
    orbit_root = -left_coefficient / right_coefficient
    assert orbit_root != 0
    singular_orbit_records.append(
        {
            "label": label,
            "face": face,
            "dual_edge": dual_edge,
            "direction": direction,
            "coefficient_labels": (
                coefficient_labels[left], coefficient_labels[right]
            ),
            "orbit_root": orbit_root,
        }
    )
assert len(singular_orbit_records) == 30
assert len({record["dual_edge"] for record in singular_orbit_records}) == 30


def face_lattice_points(face, facet_records):
    """Return the lattice points on a face from its supporting equalities."""
    containing = tuple(
        index
        for index, (_normal, _constant, indices) in enumerate(facet_records)
        if face.issubset(indices)
    )
    return tuple(
        point
        for point in DELTA_LATTICE_POINTS
        if all(
            dot(point, facet_records[index][0]) == facet_records[index][1]
            for index in containing
        )
    )


facet_records, faces, face_dimensions = nonempty_face_vertex_sets()
face_point_sets = {
    face: face_lattice_points(face, facet_records) for face in faces
}
for face, points in face_point_sets.items():
    vertex_points = {
        DELTA_F1_VERTICES[index] for index in face
    }
    assert vertex_points.issubset(points)


WITNESS_PARAMETER_VALUES = (
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
)
assert len(WITNESS_PARAMETER_VALUES) == len(PARAMETER_POINTS)
witness_coefficients = {point: QQ(1) for point in GAUGE_POINTS}
witness_coefficients.update(
    dict(zip(PARAMETER_POINTS, WITNESS_PARAMETER_VALUES))
)


def verify_generic_face_nondegeneracy():
    """Exhibit a nondegenerate coefficient witness for every face."""
    rational_laurent_ring = LaurentPolynomialRing(
        QQ, names=("y1", "y2", "y3", "y4")
    )
    y = rational_laurent_ring.gens()
    checked_by_dimension = {dimension: 0 for dimension in range(5)}
    for face in faces:
        # The all-one full polynomial has a quick exact dense-torus check.
        # Distinct-prime coefficients avoid its special degeneracies on
        # proper faces.  A separate witness on each face suffices because
        # the coefficient chart is irreducible and there are finitely many
        # face-nondegeneracy opens.
        coefficients = (
            {point: QQ(1) for point in DELTA_LATTICE_POINTS}
            if face_dimensions[face] == 4
            else witness_coefficients
        )
        face_polynomial = sum(
            coefficients[point] * laurent_monomial(y, point)
            for point in face_point_sets[face]
        )
        logarithmic_jacobian = [
            face_polynomial,
            *(
                y[index] * face_polynomial.derivative(y[index])
                for index in range(4)
            ),
        ]
        if not rational_laurent_ring.ideal(logarithmic_jacobian).is_one():
            raise AssertionError(
                "facewise witness is degenerate on face "
                f"{tuple(index + 1 for index in sorted(face))}"
            )
        checked_by_dimension[face_dimensions[face]] += 1
    assert checked_by_dimension == {0: 22, 1: 68, 2: 72, 3: 26, 4: 1}
    return checked_by_dimension


def verify_common_witness_full_face_nondegeneracy():
    """Certify the distinct-prime witness on the full-dimensional face."""
    rational_laurent_ring = LaurentPolynomialRing(
        QQ, names=("u1", "u2", "u3", "u4")
    )
    u = rational_laurent_ring.gens()
    witness_polynomial = sum(
        witness_coefficients[point] * laurent_monomial(u, point)
        for point in DELTA_LATTICE_POINTS
    )
    logarithmic_jacobian = [
        witness_polynomial,
        *(
            u[index] * witness_polynomial.derivative(u[index])
            for index in range(4)
        ),
    ]
    if not rational_laurent_ring.ideal(logarithmic_jacobian).is_one():
        raise AssertionError(
            "the distinct-prime coefficient witness is degenerate on "
            "the full-dimensional face"
        )
    return True


nondegenerate_face_counts = (
    verify_generic_face_nondegeneracy()
    if GLOBAL_SECTION_VERIFY_FACES
    else None
)
common_witness_full_face_nondegenerate = (
    verify_common_witness_full_face_nondegeneracy()
    if GLOBAL_SECTION_VERIFY_COMMON_WITNESS
    else None
)

if not GLOBAL_SECTION_QUIET:
    print("global_section: all exact assertions passed")
    print("Newton lattice points: 25")
    print(f"normalized coefficient points: {GAUGE_POINTS}")
    print("generic coefficient parameters: 20")
    if GLOBAL_SECTION_VERIFY_FACES:
        print(
            "face counts checked for Delta-regularity: "
            f"{nondegenerate_face_counts}"
        )
    if GLOBAL_SECTION_VERIFY_COMMON_WITNESS:
        print(
            "distinct-prime common witness: full-face nondegeneracy passed"
        )
    print("parameter-to-exponent dictionary:")
    for name, point in zip(parameter_names, PARAMETER_POINTS):
        print(f"  {name}: {point}")
    print("singular-orbit dual edges and roots:")
    for record in singular_orbit_records:
        left, right = record["dual_edge"]
        left_label, right_label = record["coefficient_labels"]
        print(
            f"  {record['label']}: edge={left}->{right}, "
            f"direction={record['direction']}, "
            f"coefficients=({left_label},{right_label}), "
            f"root=-({left_label})/({right_label})"
        )
