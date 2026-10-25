#!/usr/bin/env python3
"""Exact certificate for two components of the X9 four-node locus.

Run: sage -python examples/verify_x9_components.py --output /tmp/x9_components.json
Requires SageMath. Reads the accompanying X9 fan; no historical scripts.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

from sage.all import QQ, ZZ, NumberField, PolynomialRing, matrix, prod, vector


def verify(data_path):
    raw = data_path.read_bytes()
    data = json.loads(raw)
    points = [tuple(p) for p in data["rays"]] + [(0, 0, 0, 0)]
    expected = [
        (1, 0, 0, 0), (0, 1, 0, 0), (-1, -1, 0, 0),
        (0, 0, 1, 0), (0, 0, 0, 1), (-4, -2, -1, 0),
        (-4, -2, 0, -1), (3, 1, 1, 1), (2, 1, 1, 1),
        (-1, 0, 0, 0), (-2, -1, 0, 0), (0, 0, 0, 0),
    ]
    assert points == expected
    coefficients = [1, 1, 668, 324, 324, 1944, 1944, 108, 432, -21, 2496, 0]
    rays = [vector(ZZ, p) for p in data["triangulation"]["points"]]
    cones = [tuple(c) for c in data["triangulation"]["tetrahedra"]]
    lattice_points = [
        p for p in itertools.product(range(-4, 4), range(-2, 2), range(-1, 2), range(-1, 2))
        if all(vector(ZZ, p) * r >= -1 for r in rays)
    ]
    assert set(lattice_points) == set(points)
    assert len(rays) == 161 and len(cones) == 813
    assert all(abs(matrix(ZZ, [rays[i] for i in c]).det()) == 1 for c in cones)

    # Character basis x=t1, y=t1*t2, z=t1^2*t2*t3, w=t1^2*t2*t4.
    change = matrix(ZZ, [(1, 0, 0, 0), (1, 1, 0, 0), (2, 1, 1, 0), (2, 1, 0, 1)])
    assert change.det() == 1
    R = PolynomialRing(QQ, names=("x", "y", "z", "w"))
    x, y, z, w = R.gens()
    F = y**2 + (x**2 - 21)*y + 668*x + 2496 + 324*(z+w) + 1944*(1/z+1/w) + 108*(x+4)*z*w
    transformed = R.fraction_field().zero()
    for n, c in zip(points, coefficients):
        exponents = (vector(ZZ, n) + vector(ZZ, (2, 1, 0, 0))) * change.inverse()
        assert all(e in ZZ for e in exponents)
        transformed += c * prod(R.fraction_field()(a)**int(e) for a, e in zip(R.gens(), exponents))
    assert transformed == F
    assert z*F.derivative(z)-w*F.derivative(w) == 324*(z-w)*(1+6/(z*w))

    # The saturated singular scheme equals the two explicitly printed pairs.
    P = PolynomialRing(QQ, names=("x", "y", "z", "w", "iv"), order="degrevlex")
    xx, yy, zz, ww, iv = P.gens()
    G = P(z*w*F)
    inverse = iv*xx*yy*zz*ww-1
    J = P.ideal([G] + [G.derivative(a) for a in (xx, yy, zz, ww)] + [inverse])
    off = P.ideal([xx-5, yy+2, zz+ww-18, zz*ww+6, inverse])
    diagonal = P.ideal([zz-ww, ww**2+3, xx+4-3*ww, yy-16-12*ww, inverse])
    assert J == off.intersection(diagonal)
    assert off + diagonal == P.ideal(1)
    assert J.dimension() == 0 and J.vector_space_dimension() == 4
    assert J.radical() == J
    remainder = [J.reduce(a) for a in (P.one(), xx, zz, ww)]
    monomials = sorted(set(m for a in remainder for m in a.monomials()), key=str)
    evaluation = matrix(QQ, [[a.monomial_coefficient(m) for m in monomials] for a in remainder])
    assert evaluation.rank() == 4

    U = PolynomialRing(QQ, "t")
    t = U.gen()
    P9 = 4*t**9+12*t**7+3*t**6-32*t**5-26*t**4+48*t**3+60*t**2-72
    Q9 = 48*t**9-27*t**8+288*t**7+144*t**6-448*t**5-408*t**4+384*t**3+576*t**2-432
    P11, Q13 = (t**2+3)*P9, (t**2+3)**2*Q9
    xp = 18/z**3-3/z-4
    yp = -(xp**2-21)/2
    p_actual = U(F.derivative(x)(x=xp, y=yp, w=z).numerator().univariate_polynomial())
    q_actual = U(F(x=xp, y=yp, w=z).numerator().univariate_polynomial())
    assert p_actual.monic() == P11.monic() and q_actual.monic() == Q13.monic()
    assert P11.gcd(Q13).monic() == t**2+3
    resultants = [P9.resultant(Q9), P9.resultant(t**2+3)]
    assert resultants == [2**59*3**13*13**2*29**2*31*9311, 3**7*13*31]
    hess = matrix([[F.derivative(a).derivative(b) for b in R.gens()] for a in R.gens()])
    hessian_values = []
    for minimal, make_point, expected_det in [
        (t**2-87, lambda r: [5, -2, 9+r, 9-r], lambda r: 2**6*3**9*87),
        (t**2+3, lambda r: [-4+3*r, 16+12*r, r, r], lambda r: -2**6*3**9*(21+16*r)),
    ]:
        K = NumberField(minimal, "r")
        r = K.gen()
        point = make_point(r)
        assert all(a != 0 for a in point)
        assert all(f(*point) == 0 for f in [F] + [F.derivative(a) for a in R.gens()])
        determinant = matrix(K, [[h(*point) for h in row] for row in hess.rows()]).det()
        assert determinant == expected_det(r) and determinant != 0
        hessian_values.append(str(determinant))

    # Independent complete boundary-orbit check in smooth maximal-cone charts.
    face_chart = {}
    for c in cones:
        for k in range(1, 5):
            for face in itertools.combinations(c, k):
                face_chart.setdefault(tuple(sorted(face)), c)
    B = PolynomialRing(QQ, names=("b0", "b1", "b2", "b3", "inv"), order="degrevlex")
    b = B.gens()[:4]
    inv = B.gen(4)
    census = {k: {"total": 0, "independent": 0} for k in range(1, 5)}
    exceptions = []
    for face, containing in sorted(face_chart.items(), key=lambda item: (len(item[0]), item[0])):
        k = len(face)
        census[k]["total"] += 1
        support = [j for j, n in enumerate(points) if coefficients[j] and all(vector(ZZ, n)*rays[i] == -1 for i in face)]
        if matrix(QQ, [[1]+list(points[j]) for j in support]).rank() == len(support):
            census[k]["independent"] += 1
            continue
        chart = list(face)+[i for i in containing if i not in face]
        poly = B.zero()
        derivative_supports = [[] for _ in range(k)]
        for j, n in enumerate(points):
            ex = [int(vector(ZZ, n)*rays[i]+1) for i in chart]
            assert min(ex) >= 0
            poly += coefficients[j]*prod(b[i]**ex[i] for i in range(4))
            for normal in range(k):
                if ex[normal] == 1 and all(ex[i] == 0 for i in range(k) if i != normal) and coefficients[j]:
                    derivative_supports[normal].append(j)
        substitution = {b[i]: 0 for i in range(k)}
        f0 = poly.subs(substitution)
        first_jet = [f0]+[poly.derivative(b[i]).subs(substitution) for i in range(k)]+[f0.derivative(b[i]) for i in range(k, 4)]
        ideal = B.ideal(first_jet+[inv*prod(b[i] for i in range(k, 4))-1])
        assert list(ideal.groebner_basis()) == [1], face
        exceptions.append({"cone": list(face), "support": support, "normal_supports": derivative_supports, "chart": chart, "first_jet": [str(f) for f in first_jet], "groebner_basis": ["1"]})
    assert [(census[k]["total"], census[k]["independent"]) for k in range(1, 5)] == [(161, 157), (974, 971), (1626, 1626), (813, 813)]
    assert [e["cone"] for e in exceptions] == [[5], [6], [8], [11], [5, 11], [6, 11], [8, 11]]
    for exception in exceptions[:2]:
        other = [j for j in exception["support"] if j != 0]
        assert matrix(QQ, [[1]+list(points[j]) for j in other]).rank()+1 == matrix(QQ, [[1]+list(points[j]) for j in exception["support"]]).rank()
    assert exceptions[4]["normal_supports"][1] == [0]
    assert exceptions[5]["normal_supports"][1] == [0]
    assert exceptions[6]["normal_supports"][1] == [9]
    return {
        "input_sha256": hashlib.sha256(raw).hexdigest(),
        "coefficients": coefficients,
        "lattice_points": points,
        "all_maximal_cones_unimodular": True,
        "singular_scheme": {"equals_displayed_pair_intersection": True, "dimension": 0, "degree": 4, "radical": True},
        "diagonal_resultants": [str(v) for v in resultants],
        "hessian_determinants": hessian_values,
        "coefficient_evaluation_rank": int(evaluation.rank()),
        "coefficient_normal_forms": [str(f) for f in remainder],
        "boundary_census": census,
        "boundary_exceptions": exceptions,
        "conclusion": "Exactly four compact ordinary double points; smoothing-parameter rank four.",
        "limits": "The connected-component distinction additionally uses the proved Paper 6 smoothing-rank and local-system results. This calculation does not classify all components.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path(__file__).resolve().with_name("x9.json"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify(args.data)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(result["conclusion"])
    print("Boundary orbits:", sum(row["total"] for row in result["boundary_census"].values()))
    print("All exact checks passed.")
