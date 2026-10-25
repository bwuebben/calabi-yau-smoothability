#!/usr/bin/env python3
"""Exact finite checks for the X9 component and reverse-transition results.

Run with SageMath:
    sage -python examples/verify_x9_global.py --output x9_global.json

Only the accompanying x9.json is read. The coefficient-space argument,
integral topology, projective blowup and analytic gluing are written proofs;
this program verifies their explicit polynomial, lattice and fan inputs.
"""
import argparse
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path

from sage.all import (
    QQ, ZZ, LaurentPolynomialRing, Polyhedron, PolynomialRing,
    gcd, matrix, prod, vector,
)


NEWTON_POINTS = [
    (1, 0, 0, 0), (0, 1, 0, 0), (-1, -1, 0, 0),
    (0, 0, 1, 0), (0, 0, 0, 1), (-4, -2, -1, 0),
    (-4, -2, 0, -1), (3, 1, 1, 1), (2, 1, 1, 1),
    (-1, 0, 0, 0), (-2, -1, 0, 0), (0, 0, 0, 0),
]
ESSENTIAL = {0, 1, 3, 4, 5, 6, 7}
MUTANT_POINTS = [
    (0, 0, 0, 1), (4, 1, 1, 1), (1, 0, 0, 0),
    (-1, 0, 0, 0), (0, 0, 0, 0), (-1, 0, 0, -1),
    (-2, 0, 0, -1), (-4, -1, 0, -1), (-4, 0, -1, -1),
]
SMOOTHING_RAYS = [
    (-1, -3, -1, 0), (-1, -3, 0, -1), (-1, -1, 0, 0),
    (-1, 0, 0, 0), (0, 1, 0, 0), (1, -1, 0, 0),
    (1, 3, 1, 1),
]
SMOOTHING_CONES = [
    (1, 4, 5, 6), (1, 3, 4, 6), (0, 4, 5, 6), (0, 3, 4, 6),
    (0, 1, 4, 5), (0, 1, 3, 4), (0, 1, 2, 5, 6),
    (1, 2, 3, 6), (0, 2, 3, 6), (0, 1, 2, 3),
]


def integer_rows(rows):
    return [[int(x) for x in row] for row in rows]


def fan_data(data):
    points = [tuple(p) for p in data["rays"]] + [(0, 0, 0, 0)]
    assert points == NEWTON_POINTS
    rays = [vector(ZZ, p) for p in data["triangulation"]["points"]]
    maximal = sorted(tuple(sorted(c)) for c in data["triangulation"]["tetrahedra"])
    assert len(rays) == 161 and len(maximal) == len(set(maximal)) == 813
    assert all(abs(matrix(ZZ, [rays[i] for i in c]).det()) == 1 for c in maximal)
    lattice_points = [p for p in itertools.product(
        range(-4, 4), range(-2, 2), range(-1, 2), range(-1, 2))
        if all(vector(ZZ, p) * r >= -1 for r in rays)]
    assert set(lattice_points) == set(points)
    faces = sorted({s for c in maximal for k in range(1, 5)
                    for s in itertools.combinations(c, k)}, key=lambda s: (len(s), s))
    census = dict(sorted(Counter(map(len, faces)).items()))
    assert census == {1: 161, 2: 974, 3: 1626, 4: 813}
    assert rays[8] == vector(ZZ, (1, -1, -1, -1))
    assert rays[11] == vector(ZZ, (0, 1, -1, -1))
    change = matrix(ZZ, [
        (1, 0, 0, 0), (1, 1, 0, 0),
        (2, 1, 1, 0), (2, 1, 0, 1),
    ])
    assert change.det() == 1
    assert change * rays[8] == vector(ZZ, (1, 0, 0, 0))
    assert change * rays[11] == vector(ZZ, (0, 1, 0, 0))
    monomials = [(vector(ZZ, p) + vector(ZZ, (2, 1, 0, 0))) * change.inverse()
                 for p in points]
    assert integer_rows(monomials) == [
        [2, 1, 0, 0], [0, 2, 0, 0], [1, 0, 0, 0],
        [0, 0, 1, 0], [0, 0, 0, 1], [0, 0, -1, 0],
        [0, 0, 0, -1], [1, 0, 1, 1], [0, 0, 1, 1],
        [0, 1, 0, 0], [0, 0, 0, 0], [1, 1, 0, 0],
    ]
    return points, rays, maximal, faces, census


def boundary_checks(points, rays, faces):
    records = {}
    for cone in faces:
        exponents = [tuple(int(vector(ZZ, n)*rays[j] + 1) for j in cone)
                     for n in points]
        assert all(min(e) >= 0 for e in exponents)
        orbit = [i for i, e in enumerate(exponents) if not any(e)]
        normal = [[i for i, e in enumerate(exponents)
                   if e == tuple(int(h == j) for h in range(len(cone)))]
                  for j in range(len(cone))]
        records[cone] = (orbit, normal)

    forced = []
    for i, cone in [(0, (19, 120)), (1, (136, 141)),
                    (3, (112, 134, 155)), (4, (121, 144, 149)),
                    (5, (127, 149, 150)), (6, (8, 132, 157))]:
        orbit, normal = records[cone]
        assert orbit == [i] and all(not e for e in normal)
        forced.append({"coefficient": i, "cone": cone,
                       "orbit_dimension": 4-len(cone),
                       "orbit_support": orbit, "normal_supports": normal})
    cone = (136, 154, 156)
    orbit, normal = records[cone]
    assert orbit == [7, 8] and all(not e for e in normal)
    both_zero = {"cone": cone, "orbit_support": orbit, "normal_supports": normal}

    fixed = []
    for cone in [(9, 80, 114, 136), (10, 85, 117, 136),
                 (114, 131, 136, 154), (117, 131, 136, 156)]:
        orbit, normal = records[cone]
        assert orbit == [7] and all(not e for e in normal)
        quadratic = {i: tuple(int(vector(ZZ, n)*rays[j]+1) for j in cone)
                     for i, n in enumerate(points)
                     if sum(vector(ZZ, n)*rays[j]+1 for j in cone) == 2}
        assert set(quadratic) == {0, 8}
        assert tuple(a+b for a, b in zip(quadratic[0], quadratic[8])) == (1, 1, 1, 1)
        fixed.append({"cone": cone, "orbit_support": orbit,
                      "normal_supports": normal, "quadratic_support": quadratic})
    binomials = []
    for cone, expected in [((9, 114, 154), [3, 8]),
                           ((10, 117, 156), [4, 8]),
                           ((9, 80, 136), [1, 8]),
                           ((10, 85, 136), [1, 8])]:
        orbit, normal = records[cone]
        assert orbit == [7]
        assert [e for e in normal if e] == [expected]
        difference = vector(ZZ, points[expected[0]]) - vector(ZZ, points[expected[1]])
        assert difference != 0 and gcd(difference) == 1
        assert all(difference*rays[j] == 0 for j in cone)
        binomials.append({"cone": cone, "orbit_support": orbit,
                          "normal_supports": normal,
                          "primitive_exponent_difference": list(map(int, difference))})

    counts = {"essential_affine_coloop": 0, "essential_singleton_normal": 0}
    chart_cones = []
    for cone in faces:
        if set(cone) <= {8, 11}:
            chart_cones.append(cone)
            continue
        orbit, normal = records[cone]
        rank = matrix(QQ, [[1]+list(points[i]) for i in orbit]).rank()
        if any(i in ESSENTIAL and matrix(QQ, [[1]+list(points[j]) for j in orbit if j != i]).rank() < rank
               for i in orbit):
            counts["essential_affine_coloop"] += 1
        else:
            assert any(len(e) == 1 and e[0] in ESSENTIAL for e in normal), cone
            counts["essential_singleton_normal"] += 1
    assert chart_cones == [(8,), (11,), (8, 11)]
    assert counts == {"essential_affine_coloop": 3555, "essential_singleton_normal": 16}
    return {"essential_nonzero_coefficients": sorted(ESSENTIAL),
            "forced_positive_dimensional_orbits": forced,
            "c7_c8_zero_witness": both_zero,
            "c7_zero_fixed_points": fixed,
            "c7_zero_c8_nonzero_binomial_orbits": binomials,
            "joint_chart_cones": chart_cones, "uniform_exclusion_counts": counts,
            "unexcluded_orbits_outside_joint_chart": []}


def compact_witness():
    P = PolynomialRing(QQ, names=("x", "y", "z", "w", "iv"), order="degrevlex")
    x, y, z, w, iv = P.gens()
    # zw is a unit on the joint chart, so multiplying F by zw preserves its singular scheme.
    poly = z*w*(y*y+(x*x+1)*y) + (z*w+1)*(z*w*(x+1)+z+w)
    inverse = iv*z*w-1
    singular = P.ideal([poly]+[poly.derivative(a) for a in (x, y, z, w)]+[inverse])
    expected = P.ideal([y, x*x+1, z+w-x-1, z*w+1, inverse])
    assert singular == expected
    assert singular.radical() == singular and singular.vector_space_dimension() == 4
    determinant = matrix(P, [[poly.derivative(a).derivative(b) for b in (x, y, z, w)]
                             for a in (x, y, z, w)]).det()
    assert singular + P.ideal(determinant) == P.ideal(1)
    coefficients = [x*x*y, y*y, x, z, w, iv*w, iv*z, x*z*w, z*w, y, P.one(), x*y]
    normal = [singular.reduce(f) for f in coefficients]
    monomials = sorted({m for f in normal for m in f.monomials()}, key=str)
    evaluation = matrix(QQ, [[f.monomial_coefficient(m) for m in monomials] for f in normal])
    assert evaluation.rank() == 3
    return {"coefficients": [1]*11+[0], "polynomial_after_multiplication_by_zw": str(poly),
            "singular_groebner_basis": [str(f) for f in singular.groebner_basis()],
            "radical": True, "degree": 4,
            "hessian_determinant_modulo_singular_ideal": str(singular.reduce(determinant)),
            "hessian_is_unit_modulo_singular_ideal": True,
            "coefficient_normal_forms": [str(f) for f in normal],
            "coefficient_evaluation_rank": int(evaluation.rank()),
            "compactness_argument": "The uniform boundary exclusion and this full joint-chart ideal exhaust the compact singular locus."}


def surface_checks(rays, maximal):
    ridx = 11
    m = vector(ZZ, (4, 2, 1, 1))
    quotient = matrix(ZZ, [(1, 0, 0, 0), (0, 1, 1, 0), (0, 1, 0, 1)])
    assert quotient*rays[ridx] == 0
    assert quotient.right_kernel().rank() == 1
    assert abs(quotient.matrix_from_columns((0, 2, 3)).det()) == 1
    assert vector(ZZ, (4, 1, 1))*quotient == m
    kernel_basis = matrix(ZZ, [(1, 0, -4), (0, 1, -1)])
    assert kernel_basis*vector(ZZ, (4, 1, 1)) == 0
    assert abs(kernel_basis.matrix_from_columns((0, 1)).det()) == 1
    star = [c for c in maximal if ridx in c]
    assert len(star) == 18
    fibre_rays = {}
    fibre_cones = set()
    for cone in star:
        values = [m*rays[i] for i in cone if i != ridx]
        assert not min(values) < 0 < max(values), cone
        zeros = [i for i in cone if i != ridx and m*rays[i] == 0]
        for i in zeros:
            q = quotient*rays[i]
            assert q[2] == -4*q[0]-q[1]
            fibre_rays[i] = (int(q[0]), int(q[1]))
        fibre_cones.update(itertools.combinations(sorted(zeros), 2))
    order = [5, 8, 158, 6, 131]
    assert [fibre_rays[i] for i in order] == [(-1, 1), (1, -2), (0, 1), (-1, 3), (-1, 2)]
    assert fibre_cones == {tuple(sorted((order[i], order[(i+1) % 5]))) for i in range(5)}
    assert all(matrix(ZZ, [fibre_rays[order[i]], fibre_rays[order[(i+1) % 5]]]).det() == 1
               for i in range(5))
    intersections = [-int(matrix(ZZ, [fibre_rays[order[i-1]], fibre_rays[order[(i+1) % 5]]]).det())
                     for i in range(5)]
    assert intersections == [0, 1, -1, -1, -2]
    return {"divisor_ray": ridx, "character": list(map(int, m)), "star_maximal_cones": star,
            "quotient_map_rows": integer_rows(quotient), "saturated_kernel_basis": integer_rows(kernel_basis),
            "fibre_ray_order": order, "fibre_rays": fibre_rays,
            "fibre_maximal_cones": sorted(fibre_cones), "self_intersections": intersections,
            "surface_K_squared": 12-len(order)}


def cayley_checks():
    rays = [vector(ZZ, v) for v in [(0, 0, 1, 0), (1, 1, 1, 0),
            (0, 0, 0, 1), (1, 0, 0, 1), (0, 1, 0, 1)]]
    M = matrix(ZZ, rays).transpose()
    relation = vector(ZZ, (-1, 1, 2, -1, -1))
    assert M*relation == 0 and M.right_kernel().rank() == 1
    cones = [(1, 2, 3, 4), (0, 1, 2, 4), (0, 1, 2, 3)]
    determinants = [int(M.matrix_from_columns(c).det()) for c in cones]
    assert determinants == [1, -1, 1]
    heights = vector(ZZ, (1, 0, 0, 1, 1))
    gaps = []
    cells = []
    for cone in cones:
        linear = matrix(QQ, [rays[i] for i in cone]).solve_right(vector(QQ, [heights[i] for i in cone]))
        omitted = [i for i in range(5) if i not in cone]
        gap = [heights[i]-linear*rays[i] for i in omitted]
        assert gap == [3]
        gaps.append(list(map(int, gap)))
        pairs = [rays[a]+rays[b] for a in cone if a < 2 for b in cone if b >= 2]
        cells.append(sorted({(int(v[0]), int(v[1])) for v in pairs}))
    assert cells == [[(1, 1), (1, 2), (2, 1)],
                     [(0, 0), (0, 1), (1, 1), (1, 2)],
                     [(0, 0), (1, 0), (1, 1), (2, 1)]]
    canonical = [sorted((1-x, 1-y) for x, y in cell) for cell in cells]
    assert canonical == [sorted(c) for c in [
        [(0, 0), (-1, 0), (0, -1)], [(0, -1), (1, 0), (1, 1), (0, 0)],
        [(1, 1), (0, 1), (-1, 0), (0, 0)]]]
    pentagon = Polyhedron(vertices=[(0, 0), (1, 0), (2, 1), (1, 2), (0, 1)])
    assert pentagon == Polyhedron(vertices=[(0, 0), (1, 1)]) + Polyhedron(vertices=[(0, 0), (1, 0), (0, 1)])
    assert sum(Polyhedron(vertices=c).volume() for c in cells) == pentagon.volume()
    R = PolynomialRing(QQ, names=("a", "b", "U", "V"))
    a, b, U, V = R.gens()
    assert R.ideal([(U-a*b*V).derivative(q) for q in R.gens()]).is_one()
    local = a*U-b*V
    assert R.ideal([local.derivative(q) for q in R.gens()]) == R.ideal(R.gens())
    assert matrix([[local.derivative(q).derivative(qq) for qq in R.gens()] for q in R.gens()]).det() == 1
    return {"rays": integer_rows(rays), "primitive_relation": list(map(int, relation)),
            "maximal_cones": cones, "determinants": determinants,
            "coherent_heights": list(map(int, heights)), "strict_convexity_gaps": gaps,
            "central_slice_cells": cells, "canonical_cells_T3_P4_P1": canonical,
            "local_critical_points": 2, "node_hessian_determinants": [1, 1]}


def mutant_checks():
    polytope = Polyhedron(vertices=MUTANT_POINTS, base_ring=QQ)
    assert len(polytope.vertices()) == 7
    assert sorted(map(tuple, polytope.integral_points())) == sorted(MUTANT_POINTS)
    transform = matrix(ZZ, [(0, 0, 0, 1), (1, 0, 0, -1), (0, 1, 0, 0), (0, 0, 1, 0)])
    assert transform.det() == -1
    assert sorted(tuple(transform*vector(v)) for v in polytope.vertices_list()) == sorted(SMOOTHING_RAYS)
    fan_polytope = Polyhedron(vertices=SMOOTHING_RAYS, base_ring=QQ)
    facets = []
    polar_vertices = []
    for inequality in fan_polytope.inequalities():
        assert inequality[0] > 0
        normal = vector(QQ, [x/inequality[0] for x in inequality[1:]])
        assert all(x in ZZ for x in normal)
        polar_vertices.append(tuple(map(int, normal)))
        facets.append(tuple(i for i, r in enumerate(SMOOTHING_RAYS) if 1+normal*vector(r) == 0))
    assert sorted(facets) == sorted(SMOOTHING_CONES)
    assert sum(fan_polytope.interior_contains(p) for p in fan_polytope.integral_points()) == 1
    polar = Polyhedron(vertices=polar_vertices, base_ring=QQ)
    assert len(polar.integral_points()) == 162

    R = LaurentPolynomialRing(QQ, names=("a", "b", "c", "v"))
    a, b, c, v = R.gens()
    G = 1+a**4*b*c
    H = a**(-1)+a**(-2)+a**(-4)*(b**(-1)+c**(-1))
    K = a+a**(-1)
    mutant = G*v+K+H/v
    assert set(map(tuple, mutant.exponents())) == set(MUTANT_POINTS)-{(0, 0, 0, 0)}
    fraction = R.fraction_field()
    t = fraction(G*v)
    assert t+fraction(K)+fraction(G*H)/t == fraction(mutant)
    jacobian = matrix(fraction, [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0],
        [a*G.derivative(a)/G, b*G.derivative(b)/G, c*G.derivative(c)/G, 1]])
    assert jacobian.det() == 1

    # Newton nondegeneracy of the exact all-one boundary, zero-origin member.
    S = PolynomialRing(QQ, names=("a", "b", "c", "v", "iv"), order="degrevlex")
    xs = S.gens()[:4]
    iv = S.gen(4)
    face_counts = {}
    for dim in range(1, 5):
        count = 0
        for face in polytope.faces(dim):
            face_polytope = face.as_polyhedron()
            support = [e for e in MUTANT_POINTS if any(e) and face_polytope.contains(e)]
            assert support
            minimum = [min(e[j] for e in support) for j in range(4)]
            f = sum(prod(xs[j]**(e[j]-minimum[j]) for j in range(4)) for e in support)
            critical = S.ideal([f]+[f.derivative(x) for x in xs]+[iv*prod(xs)-1])
            assert critical.is_one(), (dim, support)
            count += 1
        face_counts[dim] = count
    assert face_counts == {1: 18, 2: 21, 3: 10, 4: 1}
    # Every vertex polynomial is a single monomial with nonzero coefficient.
    assert all(any(v) for v in polytope.vertices_list())
    return {"newton_lattice_points": MUTANT_POINTS,
            "newton_vertices": integer_rows(sorted(polytope.vertices_list())),
            "unimodular_map_rows": integer_rows(transform), "unimodular_determinant": -1,
            "smoothing_fan_rays": SMOOTHING_RAYS, "smoothing_fan_maximal_cones": sorted(facets),
            "polar_vertices": sorted(polar_vertices), "polar_lattice_point_count": 162,
            "laurent_substitution_verified": True, "logarithmic_jacobian_determinant": 1,
            "all_one_mutant_nondegenerate": True, "tested_face_counts": face_counts,
            "vertex_faces": len(polytope.vertices())}


def verify(data_path):
    raw = data_path.read_bytes()
    points, rays, maximal, faces, census = fan_data(json.loads(raw))
    return {
        "schema_version": 1,
        "input_sha256": hashlib.sha256(raw).hexdigest(),
        "newton_lattice_points": points,
        "ambient_ray_count": len(rays), "unimodular_maximal_cone_count": len(maximal),
        "boundary_orbit_census": census,
        "boundary": boundary_checks(points, rays, faces),
        "compact_four_node_witness": compact_witness(),
        "compact_surface": surface_checks(rays, maximal),
        "original_side_cayley_model": cayley_checks(),
        "mirror_mutation": mutant_checks(),
        "result": "All exact finite checks passed.",
        "scope": "Uniform boundary exclusions, compact four-node witness and rank, toric surface, local Cayley model, and mutant Newton family.",
        "limits": "Connected-component classification, integral relation and torsion, projective small resolution, and analytic gluing additionally use the written proofs. This program does not prove a quantum or categorical mirror statement.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path(__file__).resolve().with_name("x9.json"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify(args.data)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(result["result"])
    print("Boundary orbits: 3574; uniform exclusions outside the joint chart: 3555 + 16.")
    print("Compact witness: four ordinary nodes, coefficient rank three; surface K^2 = 7.")
    print("Cayley central cells and nine-point mutant Newton family verified.")


if __name__ == "__main__":
    main()
