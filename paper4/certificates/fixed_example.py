#!/usr/bin/env python3
"""Freeze and verify the exact geometric input for Paper 4.

The script owns a literal copy of the two vertex lists used in Paper 4. It
then re-derives reflexivity, polarity, face inventories, Hodge numbers, all
30 singular faces, the 26 nodal diagonal relations, and the two coloop rows.
Every calculation is integral or rational; there is no floating-point path.

Run:  python3 fixed_example.py
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import product
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from batyrev_global import (  # noqa: E402
    classify_polygon,
    dual_edge_length,
    face_lattice_polygon,
    facets,
    is_reflexive,
    two_faces,
)
from hodge_numbers import hodge_numbers  # noqa: E402
from paper3_node_relations import diagonal_relation, rational_rank  # noqa: E402
from toric_census import smoothing_components  # noqa: E402


# Paper 3, Theorem 6.1. This ordering is authoritative for Delta.
DELTA_F1_VERTICES = (
    (1, 0, 0, 0),
    (0, 1, 0, 0),
    (1, -1, 0, 0),
    (0, 0, 1, 0),
    (0, 0, 0, 1),
    (0, 0, 1, -1),
    (0, 0, -1, 1),
    (0, 0, 0, -1),
    (0, 0, -1, 0),
    (-1, 1, 0, 0),
    (0, -1, 0, 0),
    (-1, 1, -1, 1),
    (0, -1, 0, -1),
    (0, -1, -1, 0),
    (-1, 1, -1, 0),
    (-1, 0, 0, -1),
    (-1, 0, -1, 1),
    (-1, -1, 0, -1),
    (-2, 1, -1, 1),
    (-2, 1, -1, 0),
    (-1, -1, -1, 0),
    (-2, 0, -1, 0),
)


# Paper 4's fixed q_1,...,q_26 ordering for Delta^circ. The set is re-derived
# below from the primitive facet normals of DELTA_F1_VERTICES.
POLAR_VERTICES = (
    (-1, -1, -1, -1),  # q1
    (-1, -1, -1, 0),   # q2
    (-1, -1, 0, -1),   # q3
    (-1, -1, 0, 1),    # q4
    (-1, -1, 1, 0),    # q5
    (-1, -1, 1, 1),    # q6
    (-1, 0, -1, -1),   # q7
    (-1, 0, -1, 0),    # q8
    (-1, 0, 0, -1),    # q9
    (-1, 0, 0, 1),     # q10
    (-1, 0, 1, 0),     # q11
    (-1, 0, 1, 1),     # q12
    (0, -1, -1, -1),   # q13
    (0, -1, -1, 0),    # q14
    (0, -1, 0, 1),     # q15
    (0, -1, 0, 0),     # q16
    (0, 1, -1, -1),    # q17
    (0, 1, -1, 0),     # q18
    (0, 1, 0, -1),     # q19
    (0, 1, 0, 0),      # q20
    (1, 1, -1, -1),    # q21
    (1, 0, -1, 0),     # q22
    (0, 0, 0, -1),     # q23
    (0, 0, 0, 1),      # q24
    (0, 0, 1, 0),      # q25
    (0, 0, 1, 1),      # q26
)


# Lexicographic order is the stable Paper 4 labeling N_1,...,N_26.
NODE_FACES = (
    (1, 2, 7, 8),
    (1, 2, 13, 14),
    (1, 3, 7, 9),
    (2, 4, 8, 10),
    (2, 4, 14, 15),
    (3, 5, 9, 11),
    (3, 5, 13, 16),
    (3, 5, 23, 25),
    (3, 9, 19, 23),    # N9, coloop
    (4, 6, 10, 12),
    (4, 10, 15, 24),   # N11, coloop
    (5, 6, 11, 12),
    (5, 6, 15, 16),
    (5, 6, 25, 26),
    (7, 8, 17, 18),
    (7, 9, 17, 19),
    (10, 12, 18, 20),
    (10, 12, 24, 26),
    (11, 12, 19, 20),
    (11, 12, 25, 26),
    (13, 14, 15, 16),
    (13, 16, 23, 25),
    (15, 16, 25, 26),
    (17, 18, 19, 20),
    (18, 20, 24, 26),
    (19, 20, 25, 26),
)

DP6_FACES = (
    (1, 2, 3, 4, 5, 6),       # D_{6,1}
    (7, 8, 9, 10, 11, 12),    # D_{6,2}
)

DP7_FACES = (
    (1, 7, 13, 17, 21),       # D_{7,1}
    (2, 8, 14, 18, 22),       # D_{7,2}
)

DP_INTERIOR_POINTS = {
    "D_6,1": (-1, -1, 0, 0),
    "D_6,2": (-1, 0, 0, 0),
    "D_7,1": (0, 0, -1, -1),
    "D_7,2": (0, 0, -1, 0),
}

# The face-lattice coordinates returned by face_lattice_polygon are matched
# to the polygons in Altmann's dP6 and dP7 calculations. Affine maps act on
# column vectors as x |-> Mx+b. Matching polygon vertices is essential:
# matching only an unordered edge-vector set is too weak.
LOCAL_EDGE_MODELS = {
    "D_6,1": {
        "face": DP6_FACES[0],
        "matrix": ((1, 0), (0, 1)),
        "translation": (0, 0),
        "canonical_vertices": frozenset(
            ((0, 0), (1, 0), (2, 1), (2, 2), (1, 2), (0, 1))
        ),
    },
    "D_6,2": {
        "face": DP6_FACES[1],
        "matrix": ((1, 0), (0, 1)),
        "translation": (0, 0),
        "canonical_vertices": frozenset(
            ((0, 0), (1, 0), (2, 1), (2, 2), (1, 2), (0, 1))
        ),
    },
    "D_7,1": {
        "face": DP7_FACES[0],
        "matrix": ((1, 0), (-1, 1)),
        "translation": (-1, 1),
        "canonical_vertices": frozenset(
            ((0, 0), (1, 1), (0, 2), (-1, 2), (-1, 1))
        ),
    },
    "D_7,2": {
        "face": DP7_FACES[1],
        "matrix": ((0, 1), (-1, 0)),
        "translation": (-1, 2),
        "canonical_vertices": frozenset(
            ((0, 0), (1, 1), (0, 2), (-1, 2), (-1, 1))
        ),
    },
}

COLOOP_FACES = frozenset((NODE_FACES[8], NODE_FACES[10]))


def _apply_matrix(matrix, vector):
    """Apply a two-by-two integral matrix to a column vector."""
    return tuple(
        sum(matrix[row][column] * vector[column] for column in range(2))
        for row in range(2)
    )


def _face_records(vertices):
    """Return exact records for all two-faces of a reflexive 4-polytope."""
    facet_data = facets(vertices)
    records = []
    for face, containing in two_faces(vertices, facet_data):
        assert len(containing) == 2
        u1 = facet_data[containing[0]][0]
        u2 = facet_data[containing[1]][0]
        coordinates, edge_vectors, edge_lengths = face_lattice_polygon(
            vertices, face, u1, u2
        )
        classification = classify_polygon(edge_vectors, edge_lengths)
        records.append(
            {
                "face": tuple(i + 1 for i in sorted(face)),
                "vertices_2d": tuple(coordinates),
                "edge_vectors": tuple(edge_vectors),
                "edge_lengths": tuple(edge_lengths),
                "dual_length": dual_edge_length(u1, u2),
                **classification,
            }
        )
    return sorted(records, key=lambda record: record["face"])


def _kind(record):
    if record["status"] == "smooth":
        return "smooth triangle"
    if (
        record["k"],
        record["i"],
        record["sc"],
    ) == (4, 0, 1):
        return "node"
    if (
        record["k"],
        record["i"],
        record["sc"],
    ) == (5, 1, 1):
        return "dP7"
    if (
        record["k"],
        record["i"],
        record["sc"],
    ) == (6, 1, 2):
        return "dP6"
    if (
        record["k"] == 4
        and record["i"] == 1
        and record["status"].startswith("RIGID")
    ):
        return "F1"
    return "other"


def _lattice_points_on_face(vertices, one_based_face):
    """Enumerate all lattice points of a face using exact facet inequalities."""
    zero_based_face = frozenset(index - 1 for index in one_based_face)
    facet_data = facets(vertices)
    containing = next(
        containing
        for face, containing in two_faces(vertices, facet_data)
        if face == zero_based_face
    )
    face_vertices = [vertices[index] for index in zero_based_face]
    coordinate_ranges = [
        range(
            min(vertex[coordinate] for vertex in face_vertices),
            max(vertex[coordinate] for vertex in face_vertices) + 1,
        )
        for coordinate in range(4)
    ]
    points = []
    for point in product(*coordinate_ranges):
        evaluations = [
            sum(normal[i] * point[i] for i in range(4))
            for normal, _constant, _indices in facet_data
        ]
        if not all(
            evaluation <= constant
            for evaluation, (_normal, constant, _indices) in zip(
                evaluations, facet_data
            )
        ):
            continue
        if not all(
            Fraction(evaluations[facet_index]) == facet_data[facet_index][1]
            for facet_index in containing
        ):
            continue
        points.append(point)
    return frozenset(points)


def verify():
    """Run every Paper 4 fixed-input assertion and return the node rows."""
    delta_facets = facets(DELTA_F1_VERTICES)
    polar_facets = facets(POLAR_VERTICES)
    assert len(DELTA_F1_VERTICES) == 22
    assert len(delta_facets) == 26
    assert is_reflexive(delta_facets)
    assert len(POLAR_VERTICES) == 26
    assert len(polar_facets) == 22
    assert is_reflexive(polar_facets)

    # facets() uses <u,v> <= 1, while the papers use the polar convention
    # <q,v> >= -1. Hence polar vertices are the negatives of facet normals.
    derived_polar = {
        tuple(-coordinate for coordinate in normal)
        for normal, constant, _indices in delta_facets
        if constant == 1
    }
    assert derived_polar == set(POLAR_VERTICES)
    derived_original = {
        tuple(-coordinate for coordinate in normal)
        for normal, constant, _indices in polar_facets
        if constant == 1
    }
    assert derived_original == set(DELTA_F1_VERTICES)

    delta_records = _face_records(DELTA_F1_VERTICES)
    polar_records = _face_records(POLAR_VERTICES)
    assert len(delta_records) == 72
    assert len(polar_records) == 68
    assert all(record["dual_length"] == 1 for record in delta_records)
    assert all(record["dual_length"] == 1 for record in polar_records)

    assert Counter(_kind(record) for record in delta_records) == Counter(
        {
            "smooth triangle": 50,
            "node": 20,
            "dP6": 1,
            "F1": 1,
        }
    )
    assert Counter(_kind(record) for record in polar_records) == Counter(
        {
            "smooth triangle": 38,
            "node": 26,
            "dP7": 2,
            "dP6": 2,
        }
    )

    by_kind = {}
    for record in polar_records:
        by_kind.setdefault(_kind(record), []).append(record)
    assert tuple(record["face"] for record in by_kind["node"]) == NODE_FACES
    assert tuple(record["face"] for record in by_kind["dP6"]) == DP6_FACES
    assert tuple(record["face"] for record in by_kind["dP7"]) == DP7_FACES
    assert all(
        smoothing_components(record["edge_vectors"]) == 2
        for record in by_kind["dP6"]
    )
    assert all(
        smoothing_components(record["edge_vectors"]) == 1
        for record in by_kind["dP7"]
    )
    records_by_face = {record["face"]: record for record in polar_records}
    for model in LOCAL_EDGE_MODELS.values():
        matrix = model["matrix"]
        determinant = (
            matrix[0][0] * matrix[1][1]
            - matrix[0][1] * matrix[1][0]
        )
        assert abs(determinant) == 1
        transformed_vertices = {
            tuple(
                coordinate + shift
                for coordinate, shift in zip(
                    _apply_matrix(matrix, vertex),
                    model["translation"],
                )
            )
            for vertex in records_by_face[model["face"]]["vertices_2d"]
        }
        assert transformed_vertices == model["canonical_vertices"]

    divisorial_faces = {
        "D_6,1": DP6_FACES[0],
        "D_6,2": DP6_FACES[1],
        "D_7,1": DP7_FACES[0],
        "D_7,2": DP7_FACES[1],
    }
    for label, face in divisorial_faces.items():
        lattice_points = _lattice_points_on_face(POLAR_VERTICES, face)
        face_vertices = {POLAR_VERTICES[index - 1] for index in face}
        assert lattice_points - face_vertices == {DP_INTERIOR_POINTS[label]}

    rows = []
    terms = []
    for face in NODE_FACES:
        zero_based = frozenset(index - 1 for index in face)
        derived_face, row = diagonal_relation(POLAR_VERTICES, zero_based)
        assert tuple(index + 1 for index in derived_face) == face
        rows.append(row)
        terms.append(
            tuple(
                (index + 1, coefficient)
                for index, coefficient in enumerate(row)
                if coefficient
            )
        )

    rank = rational_rank(rows)
    assert rank == 18
    deletion_ranks = [
        rational_rank(rows[:index] + rows[index + 1 :])
        for index in range(len(rows))
    ]
    derived_coloops = frozenset(
        NODE_FACES[index]
        for index, deletion_rank in enumerate(deletion_ranks)
        if deletion_rank < rank
    )
    assert derived_coloops == COLOOP_FACES
    assert deletion_ranks[8] == deletion_ranks[10] == 17
    assert all(
        deletion_rank == 18
        for index, deletion_rank in enumerate(deletion_ranks)
        if index not in (8, 10)
    )
    assert rational_rank(
        [row for index, row in enumerate(rows) if index not in (8, 10)]
    ) == 16

    assert hodge_numbers(DELTA_F1_VERTICES) == (20, 26)
    assert hodge_numbers(POLAR_VERTICES) == (26, 20)
    return rows, terms


def main():
    _rows, terms = verify()
    print("fixed_example: all exact assertions passed")
    print("Delta_F1: 22 vertices, 26 facets, 72 two-faces, Hodge (20,26)")
    print("Delta_F1^circ: 26 vertices, 22 facets, 68 two-faces, Hodge (26,20)")
    print("X^circ singularities: 26 nodes + 2 dP7 cones + 2 dP6 cones")
    print("node matrix: rank 18; N9 and N11 are the two coloops")
    print()
    print("nodal labels and diagonal relations:")
    for index, (face, relation_terms) in enumerate(zip(NODE_FACES, terms), 1):
        marker = " [COLOOP]" if face in COLOOP_FACES else ""
        print(f"  N{index:02d}: face={face}, rho={relation_terms}{marker}")
    print()
    print(f"D_6,1={DP6_FACES[0]}; D_6,2={DP6_FACES[1]}")
    print(f"D_7,1={DP7_FACES[0]}; D_7,2={DP7_FACES[1]}")
    print("divisorial interior points:")
    for label, point in DP_INTERIOR_POINTS.items():
        print(f"  {label}: {point}")
    print("Altmann edge-coordinate maps:")
    for label, model in LOCAL_EDGE_MODELS.items():
        matrix = model["matrix"]
        determinant = (
            matrix[0][0] * matrix[1][1]
            - matrix[0][1] * matrix[1][0]
        )
        print(
            f"  {label}: M={matrix}, b={model['translation']}, "
            f"det={determinant}"
        )


if __name__ == "__main__":
    main()
