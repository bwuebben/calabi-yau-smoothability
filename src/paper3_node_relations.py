#!/usr/bin/env python3
"""Exact node-relation check for the mirror in paper 3, Theorem 6.1.

The mirror of the printed 22-vertex Delta_{F_1} has 26 square faces, hence
26 ordinary double points, as well as four non-nodal del Pezzo-cone points.
For every square face F, this script forms its diagonal relation rho_F in
Q^26 and computes the row matroid exactly (using fractions, not floating
point arithmetic).

The output proves a statement about the *node subsystem*: two square-face
relations are coloops, so there is no dependence of the 26 diagonal rows in
which every coefficient is nonzero.  Friedman--Batyrev--Kreuzer's criterion
applies directly only to an all-nodal threefold with a projective small
resolution.  The computation therefore does not decide smoothability of the
actual mixed-singularity mirror; the four del Pezzo deformation directions
may participate in its local-to-global obstruction map.

Run from the repository root:

    python3 src/paper3_node_relations.py
"""

from fractions import Fraction
from itertools import combinations

from batyrev_global import facets, two_faces
from both_sides_census import DF1_PAPER


def rational_rank(rows):
    """Rank over Q by exact reduced row elimination."""
    matrix = [[Fraction(entry) for entry in row] for row in rows]
    if not matrix:
        return 0
    nrows, ncols = len(matrix), len(matrix[0])
    rank = 0
    for col in range(ncols):
        pivot = next((i for i in range(rank, nrows) if matrix[i][col]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        value = matrix[rank][col]
        matrix[rank] = [entry / value for entry in matrix[rank]]
        for i in range(nrows):
            if i == rank or not matrix[i][col]:
                continue
            value = matrix[i][col]
            matrix[i] = [
                matrix[i][j] - value * matrix[rank][j]
                for j in range(ncols)
            ]
        rank += 1
    return rank


def diagonal_relation(vertices, face):
    """Return the coefficient vector of the parallelogram diagonal relation."""
    face = tuple(sorted(face))
    points = [vertices[i] for i in face]
    pairings = (
        ((0, 1), (2, 3)),
        ((0, 2), (1, 3)),
        ((0, 3), (1, 2)),
    )
    matches = []
    for left, right in pairings:
        left_sum = tuple(points[left[0]][j] + points[left[1]][j] for j in range(4))
        right_sum = tuple(points[right[0]][j] + points[right[1]][j] for j in range(4))
        if left_sum == right_sum:
            matches.append((left, right))
    assert len(matches) == 1, (face, matches)

    left, right = matches[0]
    relation = [0] * len(vertices)
    for local_index in left:
        relation[face[local_index]] = 1
    for local_index in right:
        relation[face[local_index]] = -1
    assert all(
        sum(relation[i] * vertices[i][j] for i in range(len(vertices))) == 0
        for j in range(4)
    )
    return face, relation


def main():
    # facets() uses the sign convention <p,v> <= 1.  Paper 3 uses the polar
    # convention <q,v> >= -1, so q=-p.  Global negation preserves faces and
    # all diagonal relations.
    polar_vertices = [
        tuple(-coordinate for coordinate in normal)
        for normal, _constant, _indices in facets(DF1_PAPER)
    ]
    assert len(polar_vertices) == 26

    polar_facets = facets(polar_vertices)
    square_faces = sorted(
        tuple(sorted(face))
        for face, _containing_facets in two_faces(polar_vertices, polar_facets)
        if len(face) == 4
    )
    assert len(square_faces) == 26

    records = [diagonal_relation(polar_vertices, face) for face in square_faces]
    rows = [row for _face, row in records]
    full_rank = rational_rank(rows)
    assert full_rank == 18

    deletion_ranks = []
    for i in range(len(rows)):
        deletion_ranks.append(rational_rank(rows[:i] + rows[i + 1 :]))
    coloop_indices = [i for i, rank in enumerate(deletion_ranks) if rank < full_rank]
    assert len(coloop_indices) == 2
    assert [deletion_ranks[i] for i in coloop_indices] == [17, 17]
    assert rational_rank(
        [row for i, row in enumerate(rows) if i not in coloop_indices]
    ) == 16

    coloop_faces = [tuple(i + 1 for i in records[j][0]) for j in coloop_indices]
    assert set(coloop_faces) == {(3, 9, 19, 23), (4, 10, 15, 24)}

    print(f"polar vertices: {len(polar_vertices)}")
    print(f"square faces / node relations: {len(rows)}")
    print(f"rank over Q: {full_rank}")
    print("coloop square faces (one-based polar-vertex indices):")
    for j in coloop_indices:
        face, row = records[j]
        one_based_face = tuple(i + 1 for i in face)
        terms = [(i + 1, coefficient) for i, coefficient in enumerate(row) if coefficient]
        coords = [polar_vertices[i] for i in face]
        print(f"  {one_based_face}: {terms}; vertices={coords}; deletion rank=17")
    print("rank after deleting both coloops: 16")
    print("no linear dependence of all 26 rows has every coefficient nonzero")


if __name__ == "__main__":
    main()
