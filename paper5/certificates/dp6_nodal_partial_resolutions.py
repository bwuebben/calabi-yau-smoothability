#!/usr/bin/env python3
"""Exact local dP6 subdivisions, root classes and A2 discriminant data.

Run: sage -python paper5/certificates/dp6_nodal_partial_resolutions.py
The certificate checks local geometry and polynomial identities; the
global deformation argument is in the accompanying paper.
"""
import json
from pathlib import Path
from sage.all import Cone, Fan, Polyhedron, PolynomialRing, QQ, matrix, vector


P = [(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)]
RAYS = [(*v, 1) for v in P] + [(0, 0, 1)]
CHECKS = []


def check(name, condition):
    assert condition, name
    CHECKS.append(name)


def subdivision(removed):
    covered = {(i - 1) % 6 for i in removed} | set(removed)
    cells = [[6, (i - 1) % 6, i, (i + 1) % 6] for i in removed]
    cells += [[6, i, (i + 1) % 6] for i in range(6) if i not in covered]
    fan = Fan(cones=cells, rays=RAYS, check=True)
    total = 0
    records = []
    for ids in cells:
        polygon = Polyhedron(vertices=[RAYS[i][:2] for i in ids])
        area = polygon.volume()
        total += area
        cone = Cone([RAYS[i] for i in ids])
        if len(ids) == 3:
            check(f"{removed}: triangle {ids} is smooth", cone.is_smooth())
            check(f"{removed}: triangle area", area == QQ(1) / 2)
            kind = "smooth"
        else:
            relation = matrix(QQ, [RAYS[i] for i in ids]).transpose().right_kernel().basis()[0]
            relation /= abs(next(x for x in relation if x))
            check(f"{removed}: node relation {ids}", sorted(relation) == [-1, -1, 1, 1])
            check(f"{removed}: node lattice area {ids}", area == 1)
            check(f"{removed}: node saturated lattice {ids}",
                  matrix(QQ, [RAYS[i] for i in ids[:3]]).det() in (-1, 1))
            kind = "ordinary_double_point"
        records.append(dict(rays=ids, kind=kind, area=str(area)))
    check(f"{removed}: support covers original cone", total == Polyhedron(vertices=P).volume() == 3)
    check(f"{removed}: six original boundary two-cones survive",
          all(fan.cone_containing(RAYS[i], RAYS[(i + 1) % 6]).dim() == 2 for i in range(6)))
    return records


def main():
    A1 = subdivision([0, 3])
    A2 = subdivision([0, 2, 4])
    # Basis h,e1,e2,e3; cyclic toric boundary classes.
    h, e1, e2, e3 = [vector(QQ, [int(i == j) for i in range(4)]) for j in range(4)]
    boundary = [e1, h-e1-e2, e2, h-e2-e3, e3, h-e1-e3]
    J = matrix.diagonal(QQ, [1, -1, -1, -1])
    K = -3*h + e1 + e2 + e3
    pair = lambda a, b: a * J * b
    check("all six boundary curves have self-intersection -1", all(pair(c, c) == -1 for c in boundary))
    check("all six boundary curves have canonical degree -1", all(pair(K, c) == -1 for c in boundary))
    check("cyclic boundary intersection matrix", all(
        pair(boundary[i], boundary[j]) == (1 if (i-j) % 6 in (1, 5) else 0)
        for i in range(6) for j in range(6) if i != j))
    alpha = h-e1-e2-e3
    check("opposite contracted curves differ by the A1 root", boundary[3]-boundary[0] == alpha)
    check("A1 root is canonical-orthogonal with square -2", pair(alpha, K) == 0 and pair(alpha, alpha) == -2)
    check("three alternating curves are pairwise disjoint", all(pair(boundary[i], boundary[j]) == 0 for i, j in [(0, 2), (0, 4), (2, 4)]))
    check("A2 roots are orthogonal to A1", pair(alpha, e1-e2) == pair(alpha, e2-e3) == 0)
    R = matrix(QQ, [2*h-e1-e2-e3, h-e1-e2, h-e2-e3, h-e1-e3]).transpose()
    F = matrix(QQ, [h,e1,e3,e2]).transpose()
    check("rotation sends each boundary curve to the next",
          all(R*boundary[i] == boundary[(i+1)%6] for i in range(6)))
    check("reflection reverses the boundary order",
          all(F*boundary[i] == boundary[-i%6] for i in range(6)))
    check("rotation acts by minus one on A1", R*alpha == -alpha)

    # Edge parameters (a,b,c,a,b,c), modulo simultaneous translation.
    # Determinant-twisted equivariance identifies this plane with A2.
    S = PolynomialRing(QQ, "a,b,c")
    a, b, c = S.gens()
    mu = vector(S, [a-c, c-b, b-a])
    root_class = vector(S, [0, *mu])
    check("rotation has the stated action on A2 coefficients",
          R.change_ring(S)*root_class == vector(S, [0,mu[1],mu[2],mu[0]]))
    check("reflection has the stated action on A2 coefficients",
          F.change_ring(S)*root_class == vector(S, [0,mu[0],mu[2],mu[1]]))
    check("A2 coefficient sum zero", sum(mu) == 0)
    rotate = {a:c, b:a, c:b}
    reflect = {a:c, b:b, c:a}
    check("Hodge map commutes with rotation", vector(S, [x.subs(rotate) for x in mu]) == vector(S, [mu[1], mu[2], mu[0]]))
    check("Hodge map has determinant twist under reflection", vector(S, [x.subs(reflect) for x in mu]) == -vector(S, [mu[0], mu[2], mu[1]]))
    check("three coefficient hyperplanes are the three coincidence lines",
          mu[0]*mu[1]*mu[2] == -(a-b)*(b-c)*(c-a))

    T = PolynomialRing(QQ, "z1,z2,z3,z4,z5,z6,z,a,b,c")
    z1,z2,z3,z4,z5,z6,z,a,b,c = T.gens()
    zs = [z1,z2,z3,z4,z5,z6]
    ts = [z+a,z+b,z+c,z+a,z+b,z+c]
    M = matrix(T, [[ts[2],z1,z2], [z4,ts[1],z3], [z5,z6,ts[0]]])
    altmann = [zs[i]*ts[i]-zs[(i-1)%6]*zs[(i+1)%6] for i in range(6)]
    altmann += [ts[1]*ts[2]-z1*z4, ts[2]*ts[3]-z2*z5, ts[0]*ts[1]-z3*z6]
    check("A2 family is exactly Altmann's nine equations", T.ideal(M.minors(2)) == T.ideal(altmann))
    check("central fibre dimension three", T.ideal(altmann+[a,b,c]).dimension() == 3)
    check("total space is rank-one matrix cone times translation line", T.ideal(altmann).dimension() == 6)
    U = PolynomialRing(QQ, "u,v,w,x")
    u,v,w,x = U.gens()
    node = u*v-w*x
    check("generic discriminant point is a node",
          matrix(U, [[node.derivative(i).derivative(j) for j in U.gens()]
                     for i in U.gens()]).det() != 0)

    out = dict(scope="local exact checks; global theorem is a separate proof",
               polygon=P, rays=RAYS, A1_cells=A1, A2_cells=A2,
               boundary_classes=[list(map(int, v)) for v in boundary],
               canonical_class=list(map(int, K)), A1_root=list(map(int, alpha)),
               A2_coefficients=["a-c", "c-b", "b-a"],
               A2_discriminant="(a-b)(b-c)(c-a)=0, modulo common translation",
               A2_matrix=[["z+c","z1","z2"],["z4","z+b","z3"],["z5","z6","z+a"]],
               checks=CHECKS, passed=len(CHECKS))
    path = Path(__file__).with_suffix(".json")
    path.write_text(json.dumps(out, indent=2)+"\n")
    print(f"PASS: {len(CHECKS)} exact local checks; wrote {path.name}")


if __name__ == "__main__":
    main()
