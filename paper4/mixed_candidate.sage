#!/usr/bin/env sage
"""Assemble and test the candidate mixed Friedman--Gross relation matrix.

The construction combines the 26 exceptional-node curve classes with the
ten functionals on exceptional-surface Picard groups dual to Altmann's local
T^1 parameters.  It then tests the four reduced dP6 line/plane smoothing
profiles in the kernel of the transpose matrix.

This is the exact numerical matrix in cyclic-difference coordinates.  The
identification of its branch-avoidance result with the intrinsic statement
on the singularity links is Section 2.4 and Lemma 4.2 of the paper.

Run from cy_smoothing or paper4:

    sage paper4/mixed_candidate.sage
    sage mixed_candidate.sage
"""

from itertools import product
from pathlib import Path


PAPER4_ROOT = Path(
    globals().get("PAPER4_ROOT_OVERRIDE", Path(__file__).resolve().parent)
)
MIXED_CANDIDATE_QUIET = globals().get("MIXED_CANDIDATE_QUIET", False)
RESTRICTION_DATA_QUIET = True
load(str(PAPER4_ROOT / "restriction_data.sage"))


DP6_EDGE_VECTORS = (
    (1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1),
)
DP7_EDGE_VECTORS = (
    (1, 1), (-1, 1), (-1, 0), (0, -1), (1, -1),
)

# These orders start at the first vertex in the canonical polygons printed in
# local_deformations.md, so edge i runs from vertex i to vertex i+1.
CANONICAL_BOUNDARY_ORDERS = {
    "D_6,1": (1, 3, 5, 6, 4, 2),
    "D_6,2": (7, 9, 11, 12, 10, 8),
    "D_7,1": (13, 21, 17, 7, 1),
    "D_7,2": (22, 18, 8, 2, 14),
}

# Derivatives of Altmann's t_i with respect to the displayed Kuranishi
# coordinates.  A common additive t-direction has already been removed.
LOCAL_T_VECTORS = {
    "D_6,1": (
        ("s1", (0, 0, -1, 0, 0, -1)),
        ("s2", (0, -1, 0, 0, -1, 0)),
        ("s3", (0, -1, 0, -1, 0, -1)),
    ),
    "D_6,2": (
        ("s1", (0, 0, -1, 0, 0, -1)),
        ("s2", (0, -1, 0, 0, -1, 0)),
        ("s3", (0, -1, 0, -1, 0, -1)),
    ),
    "D_7,1": (
        ("alpha", (0, 0, -1, 1, -1)),
        ("beta", (0, -1, 0, 0, -1)),
    ),
    "D_7,2": (
        ("alpha", (0, 0, -1, 1, -1)),
        ("beta", (0, -1, 0, 0, -1)),
    ),
}


def edge_vectors(record):
    """Derive the canonical cyclic edge vectors from the labeled vertices."""
    order = CANONICAL_BOUNDARY_ORDERS[record["label"]]
    vertices = [vector(ZZ, record["coordinates"][ray]) for ray in order]
    return tuple(
        tuple(vertices[(index + 1) % len(vertices)] - vertices[index])
        for index in range(len(vertices))
    )


def parameter_records(surface_records):
    """Return the ten local-parameter rows in the ambient Picard dual."""
    records = []
    for surface in surface_records:
        first_record = len(records)
        label = surface["label"]
        order = CANONICAL_BOUNDARY_ORDERS[label]
        derived_edges = edge_vectors(surface)
        expected_edges = (
            DP6_EDGE_VECTORS if label.startswith("D_6") else DP7_EDGE_VECTORS
        )
        assert derived_edges == expected_edges

        boundary = surface["boundary"]
        boundary_position = {
            ray: position for position, ray in enumerate(boundary)
        }
        assert set(order) == set(boundary)
        t_rows = LOCAL_T_VECTORS[label]
        assert matrix(ZZ, [row for _name, row in t_rows]).rank() == len(t_rows)

        for parameter, t_vector in t_rows:
            assert sum(
                t_vector[index] * vector(ZZ, derived_edges[index])
                for index in range(len(t_vector))
            ) == 0

            # If d_i=v_{i+1}-v_i, then lambda_i=t_i-t_{i-1} satisfies
            # sum lambda_i=sum lambda_i v_i=0.  It is therefore a functional
            # on Pic(E) that kills K_E.
            lambdas_cyclic = tuple(
                t_vector[index] - t_vector[index - 1]
                for index in range(len(t_vector))
            )
            assert sum(lambdas_cyclic) == 0
            canonical_vertices = [
                vector(ZZ, surface["coordinates"][ray]) for ray in order
            ]
            assert sum(
                lambdas_cyclic[index] * canonical_vertices[index]
                for index in range(len(order))
            ) == 0

            lambda_by_ray = dict(zip(order, lambdas_cyclic))
            lambdas_boundary = vector(
                ZZ, [lambda_by_ray[ray] for ray in boundary]
            )
            picard_functional = vector(
                ZZ,
                [lambdas_boundary[index] for index in range(2, len(boundary))],
            )
            assert (
                picard_functional * surface["surface_grading"]
                == lambdas_boundary
            )
            assert picard_functional * surface["canonical"] == 0

            ambient_functional = (
                picard_functional * surface["picard_restriction"]
            )
            assert ambient_functional * vector(ZZ, ambient[3]) == 0
            records.append(
                {
                    "label": f"{label}:{parameter}",
                    "surface": label,
                    "parameter": parameter,
                    "t": tuple(t_vector),
                    "lambda_cyclic": lambdas_cyclic,
                    "lambda_boundary": tuple(lambdas_boundary),
                    "picard_functional": tuple(picard_functional),
                    "ambient_functional": tuple(ambient_functional),
                }
            )
        local_functionals = matrix(
            ZZ,
            [
                record["picard_functional"]
                for record in records[first_record:]
            ],
        )
        expected_rank = surface["picard_restriction"].nrows() - 1
        assert local_functionals.rank() == expected_rank
        assert smith_nonzero_diagonal(local_functionals) == (1,) * expected_rank
        local_kernel = vector(ZZ, local_functionals.right_kernel_matrix().row(0))
        canonical = vector(ZZ, surface["canonical"])
        if local_kernel[0] * canonical[0] < 0:
            local_kernel = -local_kernel
        assert local_kernel == canonical
    assert len(records) == 10
    return tuple(records)


def normalize(vector_input):
    """Return a primitive integer vector with first nonzero entry positive."""
    result = vector(ZZ, vector_input)
    divisor = gcd(abs(entry) for entry in result)
    assert divisor > 0
    result = result / divisor
    first = next(entry for entry in result if entry)
    if first < 0:
        result = -result
    return vector(ZZ, result)


def find_witness(kernel_basis, required_functionals):
    """Find a small integral kernel vector off all required hyperplanes."""
    assert kernel_basis.nrows() > 0
    for bound in (1, 2):
        coefficient_range = tuple(range(-bound, bound + 1))
        best = None
        for coefficients in product(
            coefficient_range, repeat=kernel_basis.nrows()
        ):
            if not any(coefficients):
                continue
            candidate = vector(ZZ, coefficients) * kernel_basis
            if any(functional * candidate == 0 for functional in required_functionals):
                continue
            candidate = normalize(candidate)
            key = (
                max(abs(entry) for entry in candidate),
                sum(abs(entry) for entry in candidate),
                tuple(candidate),
            )
            if best is None or key < best[0]:
                best = (key, candidate)
        if best is not None:
            return best[1]
    raise AssertionError("no small witness found despite nonzero functionals")


def profile_record(left_choice, right_choice, relation_matrix, labels):
    """Test one of the four reduced dP6 component profiles."""
    choices = {"D_6,1": left_choice, "D_6,2": right_choice}
    forbidden = {"D_7,1:alpha", "D_7,2:alpha"}
    required = [f"N_{index}" for index in range(1, 27)]
    required.extend(("D_7,1:beta", "D_7,2:beta"))
    difference_requirements = []
    for surface, choice in choices.items():
        if choice == "L":
            forbidden.update((f"{surface}:s1", f"{surface}:s2"))
            required.append(f"{surface}:s3")
        else:
            assert choice == "P"
            forbidden.add(f"{surface}:s3")
            required.extend((f"{surface}:s1", f"{surface}:s2"))
            difference_requirements.append(
                (f"{surface}:s1-s2", f"{surface}:s1", f"{surface}:s2")
            )

    label_position = {label: index for index, label in enumerate(labels)}
    allowed_indices = tuple(
        index for index, label in enumerate(labels) if label not in forbidden
    )
    allowed_position = {
        global_index: local_index
        for local_index, global_index in enumerate(allowed_indices)
    }
    restricted_relations = relation_matrix.transpose().matrix_from_columns(
        allowed_indices
    )
    kernel_basis = restricted_relations.right_kernel_matrix()
    assert restricted_relations * kernel_basis.transpose() == 0

    functionals = []
    functional_labels = []
    for required_label in required:
        functional = zero_vector(ZZ, len(allowed_indices))
        global_index = label_position[required_label]
        assert global_index in allowed_position
        functional[allowed_position[global_index]] = 1
        functionals.append(functional)
        functional_labels.append(required_label)
    for difference_label, positive_label, negative_label in difference_requirements:
        functional = zero_vector(ZZ, len(allowed_indices))
        positive_index = label_position[positive_label]
        negative_index = label_position[negative_label]
        functional[allowed_position[positive_index]] = 1
        functional[allowed_position[negative_index]] = -1
        functionals.append(functional)
        functional_labels.append(difference_label)

    # A finite union of proper hyperplanes cannot cover a rational vector
    # space.  Thus a transverse rational vector exists exactly when none of
    # the required functionals vanishes identically on this kernel.
    forced_zero = tuple(
        label
        for label, functional in zip(functional_labels, functionals)
        if kernel_basis * functional == 0
    )
    if forced_zero:
        return {
            "profile": left_choice + right_choice,
            "allowed_dimension": len(allowed_indices),
            "kernel_dimension": kernel_basis.nrows(),
            "forbidden": tuple(sorted(forbidden)),
            "required": tuple(functional_labels),
            "forced_zero": forced_zero,
            "witness": None,
        }

    witness_allowed = find_witness(kernel_basis, functionals)
    witness = zero_vector(ZZ, len(labels))
    for local_index, global_index in enumerate(allowed_indices):
        witness[global_index] = witness_allowed[local_index]
    assert relation_matrix.transpose() * witness == 0
    assert all(witness[label_position[label]] == 0 for label in forbidden)
    for label, functional in zip(functional_labels, functionals):
        assert functional * witness_allowed != 0, label

    return {
        "profile": left_choice + right_choice,
        "allowed_dimension": len(allowed_indices),
        "kernel_dimension": kernel_basis.nrows(),
        "forbidden": tuple(sorted(forbidden)),
        "required": tuple(functional_labels),
        "forced_zero": (),
        "witness": tuple(witness),
    }


parameters = parameter_records(surfaces)
parameter_labels = tuple(record["label"] for record in parameters)
parameter_matrix = matrix(
    ZZ, [record["ambient_functional"] for record in parameters]
)
assert parameter_matrix.dimensions() == (10, 26)
assert parameter_matrix.rank() == 10
assert smith_nonzero_diagonal(parameter_matrix) == (1,) * 10

local_labels = (
    tuple(f"N_{index}" for index in range(1, 27)) + parameter_labels
)
candidate_matrix = curve_map.stack(parameter_matrix)
assert candidate_matrix.dimensions() == (36, 26)
candidate_rank = candidate_matrix.rank()
assert candidate_rank == 21
candidate_smith = smith_nonzero_diagonal(candidate_matrix)
assert candidate_smith == (1,) * candidate_rank
assert candidate_matrix * vector(ZZ, ambient[3]) == 0

candidate_rows = tuple(tuple(row) for row in candidate_matrix.rows())
candidate_coloops = tuple(
    local_labels[index]
    for index in range(len(candidate_rows))
    if matrix(
        ZZ, candidate_rows[:index] + candidate_rows[index + 1:]
    ).rank() < candidate_rank
)
assert candidate_coloops == ("D_7,1:beta",)

profiles = tuple(
    profile_record(left, right, candidate_matrix, local_labels)
    for left, right in (("L", "L"), ("L", "P"), ("P", "L"), ("P", "P"))
)
assert tuple(record["forced_zero"] for record in profiles) == (
    ("N_9", "N_11", "D_7,1:beta", "D_7,2:beta"),
    ("N_9", "N_11", "D_7,1:beta", "D_7,2:beta", "D_6,1:s3"),
    ("N_9", "N_11", "D_7,1:beta", "D_7,2:beta", "D_6,2:s3"),
    ("N_9", "N_11", "D_7,1:beta", "D_7,2:beta"),
)

if not MIXED_CANDIDATE_QUIET:
    print("mixed_candidate: all exact assertions passed")
    print(
        f"candidate relation matrix: 36 x 26, rank {candidate_rank}, "
        f"kernel of transpose dimension {36 - candidate_rank}"
    )
    print(f"candidate coloops: {candidate_coloops}")
    print("local parameter rows:")
    for record in parameters:
        print(
            f"  {record['label']}: t={record['t']}, "
            f"lambda={record['lambda_cyclic']}, "
            f"Pic(E)^*={record['picard_functional']}"
        )
    print("candidate transverse profiles:")
    for record in profiles:
        print(
            f"  {record['profile']}: "
            f"allowed dimension={record['allowed_dimension']}, "
            f"kernel dimension={record['kernel_dimension']}"
        )
        if record["forced_zero"]:
            print(f"    forced zero={record['forced_zero']}")
        else:
            print(f"    witness={record['witness']}")
    print(
        "NOTE: the intrinsic kernel differs by branch-preserving local changes "
        "of coordinates; see Section 2.4 and Lemma 4.2 of the paper."
    )
