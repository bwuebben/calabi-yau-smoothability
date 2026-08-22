#!/usr/bin/env sage
"""Pull the canonical dP6 and dP7 characters to the four ambient charts.

For each divisorial face, the script derives a basis of the saturated
rank-three cone lattice, extends the two-dimensional edge-coordinate map to
an affine unimodular map of Gorenstein cones, dualizes it, and chooses a
small integral representative in the ambient character lattice Z^4.

An ambient representative is unique only modulo the primitive character
annihilating the face cone. The script prints that unit direction and checks
that every representative has the required restriction and is regular on
the face chart.

Run:  sage chart_characters.sage
"""

from pathlib import Path
import sys


PAPER4_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = PAPER4_ROOT.parent.parent
sys.path.insert(0, str(PAPER4_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from batyrev_global import facets, int_kernel, solve_int_coords, two_faces
from fixed_example import DP6_FACES, DP7_FACES, POLAR_VERTICES


DP6_CANONICAL_VERTICES = frozenset(
    ((0, 0), (1, 0), (2, 1), (2, 2), (1, 2), (0, 1))
)
DP7_CANONICAL_VERTICES = frozenset(
    ((0, 0), (1, 1), (0, 2), (-1, 2), (-1, 1))
)

DP6_CHARACTERS = (
    ("z1", (0, 1, 0)),
    ("z2", (1, 0, 0)),
    ("z3", (1, -1, 1)),
    ("z4", (0, -1, 2)),
    ("z5", (-1, 0, 2)),
    ("z6", (-1, 1, 1)),
    ("tau", (0, 0, 1)),
)

DP7_CHARACTERS = (
    ("y1", (1, 0, 1)),
    ("y2", (1, 1, 0)),
    ("y3", (0, -1, 2)),
    ("y4", (-1, -1, 2)),
    ("y5", (-1, 1, 0)),
    ("y6", (0, 1, 0)),
    ("y7", (-1, 0, 1)),
    ("tau", (0, 0, 1)),
)

FACE_MODELS = {
    "D_6,1": {
        "face": DP6_FACES[0],
        "linear_map": ((1, 0), (0, 1)),
        "canonical_vertices": DP6_CANONICAL_VERTICES,
        "characters": DP6_CHARACTERS,
    },
    "D_6,2": {
        "face": DP6_FACES[1],
        "linear_map": ((1, 0), (0, 1)),
        "canonical_vertices": DP6_CANONICAL_VERTICES,
        "characters": DP6_CHARACTERS,
    },
    "D_7,1": {
        "face": DP7_FACES[0],
        "linear_map": ((1, 0), (-1, 1)),
        "canonical_vertices": DP7_CANONICAL_VERTICES,
        "characters": DP7_CHARACTERS,
    },
    "D_7,2": {
        "face": DP7_FACES[1],
        "linear_map": ((0, 1), (-1, 0)),
        "canonical_vertices": DP7_CANONICAL_VERTICES,
        "characters": DP7_CHARACTERS,
    },
}

# These values are deliberately duplicated from the derived output below.
# They turn the script into a regression test for the chosen face bases,
# affine polygon maps, and ambient representatives rather than a printer whose
# output could change silently after an upstream edit.
EXPECTED_CHART_DATA = {
    "D_6,1": {
        "translation": (0, 0),
        "unit_direction": (-1, 1, 0, 0),
        "ambient_characters": (
            ("z1", (-1, 0, 0, 1)),
            ("z2", (-1, 0, 1, 0)),
            ("z3", (-1, 0, 1, -1)),
            ("z4", (-1, 0, 0, -1)),
            ("z5", (-1, 0, -1, 0)),
            ("z6", (-1, 0, -1, 1)),
            ("tau", (-1, 0, 0, 0)),
        ),
    },
    "D_6,2": {
        "translation": (0, 0),
        "unit_direction": (0, -1, 0, 0),
        "ambient_characters": (
            ("z1", (-1, 0, 0, 1)),
            ("z2", (-1, 0, 1, 0)),
            ("z3", (-1, 0, 1, -1)),
            ("z4", (-1, 0, 0, -1)),
            ("z5", (-1, 0, -1, 0)),
            ("z6", (-1, 0, -1, 1)),
            ("tau", (-1, 0, 0, 0)),
        ),
    },
    "D_7,1": {
        "translation": (-1, 1),
        "unit_direction": (0, 0, -1, 1),
        "ambient_characters": (
            ("y1", (1, 0, -1, 0)),
            ("y2", (0, 1, -1, 0)),
            ("y3", (1, -1, -1, 0)),
            ("y4", (0, -1, -1, 0)),
            ("y5", (-2, 1, -1, 0)),
            ("y6", (-1, 1, -1, 0)),
            ("y7", (-1, 0, -1, 0)),
            ("tau", (0, 0, -1, 0)),
        ),
    },
    "D_7,2": {
        "translation": (-1, 2),
        "unit_direction": (0, 0, 0, -1),
        "ambient_characters": (
            ("y1", (0, 1, -1, 0)),
            ("y2", (-1, 1, -1, 0)),
            ("y3", (1, 0, -1, 0)),
            ("y4", (1, -1, -1, 0)),
            ("y5", (-1, -1, -1, 0)),
            ("y6", (-1, 0, -1, 0)),
            ("y7", (0, -1, -1, 0)),
            ("tau", (0, 0, -1, 0)),
        ),
    },
}


def primitive(vector_entries):
    """Return the primitive integral vector on the same positive ray."""
    entries = tuple(ZZ(entry) for entry in vector_entries)
    divisor = gcd(abs(entry) for entry in entries)
    assert divisor > 0
    return tuple(entry // divisor for entry in entries)


def ambient_solution(basis_matrix, local_character):
    """Choose a small integral lift of a rank-three local character."""
    restriction = basis_matrix.transpose()
    candidates = []
    for free_coordinate in range(4):
        retained = [index for index in range(4) if index != free_coordinate]
        square = restriction.matrix_from_columns(retained)
        if square.det() == 0:
            continue
        solution = square.solve_right(vector(QQ, local_character))
        if not all(entry in ZZ for entry in solution):
            continue
        ambient = [ZZ(0)] * 4
        for index, entry in zip(retained, solution):
            ambient[index] = ZZ(entry)
        candidates.append(tuple(ambient))
    assert candidates
    return min(
        set(candidates),
        key=lambda candidate: (
            sum(abs(entry) for entry in candidate),
            candidate,
        ),
    )


def chart_record(label, model, facet_data, face_data):
    """Derive one complete face-coordinate and character record."""
    one_based_face = model["face"]
    face = frozenset(index - 1 for index in one_based_face)
    containing = next(containing for candidate, containing in face_data if candidate == face)
    assert len(containing) == 2
    normal_1 = vector(ZZ, facet_data[containing[0]][0])
    normal_2 = vector(ZZ, facet_data[containing[1]][0])

    difference_basis_raw = int_kernel(
        [
            [int(entry) for entry in normal_1],
            [int(entry) for entry in normal_2],
        ]
    )
    difference_basis = [vector(ZZ, entry) for entry in difference_basis_raw]
    assert len(difference_basis) == 2
    base_vertex = vector(ZZ, POLAR_VERTICES[min(face)])
    cone_basis = matrix(
        ZZ,
        4,
        3,
        lambda row, column: (
            difference_basis[0],
            difference_basis[1],
            base_vertex,
        )[column][row],
    )
    assert gcd(abs(minor) for minor in cone_basis.minors(3)) == 1

    computed_vertices = []
    for vertex_index in sorted(face):
        difference = vector(ZZ, POLAR_VERTICES[vertex_index]) - base_vertex
        coordinates = vector(
            ZZ,
            solve_int_coords(
                difference_basis_raw,
                tuple(int(entry) for entry in difference),
            ),
        )
        computed_vertices.append((vertex_index + 1, coordinates))

    linear_map = matrix(ZZ, model["linear_map"])
    candidate_translations = set()
    for _index, computed in computed_vertices:
        for canonical in model["canonical_vertices"]:
            translation = vector(ZZ, canonical) - linear_map * computed
            transformed = {
                tuple(linear_map * coordinates + translation)
                for _vertex_index, coordinates in computed_vertices
            }
            if transformed == model["canonical_vertices"]:
                candidate_translations.add(tuple(translation))
    if len(candidate_translations) != 1:
        print(
            "translation diagnostic",
            label,
            computed_vertices,
            model["linear_map"],
            candidate_translations,
        )
    assert len(candidate_translations) == 1
    translation = vector(ZZ, candidate_translations.pop())

    cone_map = block_matrix(
        ZZ,
        [
            [linear_map, matrix(ZZ, 2, 1, list(translation))],
            [matrix(ZZ, 1, 2, [0, 0]), matrix(ZZ, 1, 1, [1])],
        ],
    )
    assert abs(cone_map.det()) == 1
    for _index, computed in computed_vertices:
        image = cone_map * vector(ZZ, tuple(computed) + (1,))
        assert tuple(image[:2]) in model["canonical_vertices"]
        assert image[2] == 1

    unit_direction = primitive(normal_1 - normal_2)
    assert cone_basis.transpose() * vector(ZZ, unit_direction) == 0

    ambient_characters = []
    for name, canonical_character in model["characters"]:
        local_character = cone_map.transpose() * vector(ZZ, canonical_character)
        ambient = ambient_solution(cone_basis, local_character)
        assert cone_basis.transpose() * vector(ZZ, ambient) == local_character

        face_values = tuple(
            sum(
                ambient[coordinate] * POLAR_VERTICES[index - 1][coordinate]
                for coordinate in range(4)
            )
            for index in one_based_face
        )
        assert min(face_values) >= 0
        ambient_characters.append(
            {
                "name": name,
                "canonical": canonical_character,
                "local": tuple(local_character),
                "ambient": ambient,
                "face_values": face_values,
            }
        )

    return {
        "label": label,
        "face": one_based_face,
        "normals": (tuple(normal_1), tuple(normal_2)),
        "difference_basis": tuple(tuple(entry) for entry in difference_basis),
        "base_vertex": tuple(base_vertex),
        "translation": tuple(translation),
        "cone_map": cone_map,
        "unit_direction": unit_direction,
        "characters": ambient_characters,
    }


facet_data = facets(POLAR_VERTICES)
face_data = two_faces(POLAR_VERTICES, facet_data)
records = [
    chart_record(label, model, facet_data, face_data)
    for label, model in FACE_MODELS.items()
]

for record in records:
    expected = EXPECTED_CHART_DATA[record["label"]]
    assert record["translation"] == expected["translation"]
    assert record["unit_direction"] == expected["unit_direction"]
    assert tuple(
        (character["name"], character["ambient"])
        for character in record["characters"]
    ) == expected["ambient_characters"]

print("chart_characters: all exact assertions passed")
for record in records:
    print()
    print(
        f"{record['label']}: face={record['face']}, "
        f"translation={record['translation']}, "
        f"unit direction={record['unit_direction']}"
    )
    print(f"  difference basis={record['difference_basis']}")
    print(f"  base vertex={record['base_vertex']}")
    for character in record["characters"]:
        print(
            f"  {character['name']}: canonical={character['canonical']}, "
            f"ambient={character['ambient']}"
        )
