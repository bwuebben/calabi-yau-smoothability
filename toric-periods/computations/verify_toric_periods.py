#!/usr/bin/env python3
"""Exact Altmann moment-algebra checks. Requires SymPy; no random search/data files.

Run: python verify_toric_periods.py
All computations are over QQ or explicitly displayed rational-function fields.
The full moment ideal of an n-edge polygon is generated through degree n-1:
apply the characteristic polynomial of the t_i to the moments, using M_0=0.

Tests:
* a primitive hexagon with length-10 Altmann base, m^3 != 0 and m^4 = 0;
* an identity-derivative quadratic shear that fails on that actual base;
* source-base-dependent common cubic H and total-family unit identities;
* Altmann's Q8: exhaustive full ideal, its radical, and N^2 != 0, N^3 = 0.
"""
from itertools import combinations, combinations_with_replacement
from math import gcd
import sympy as s


def edge_vectors(vertices):
    n = len(vertices)
    return [(vertices[(i+1) % n][0]-p[0],
             vertices[(i+1) % n][1]-p[1]) for i, p in enumerate(vertices)]


def check_primitive_convex(vertices):
    edges = edge_vectors(vertices)
    assert all(gcd(abs(a), abs(b)) == 1 for a, b in edges)
    # All nonincident vertices lie strictly to the left of each oriented edge.
    n = len(vertices)
    for i, ((dx, dy), (px, py)) in enumerate(zip(edges, vertices)):
        for j, (qx, qy) in enumerate(vertices):
            if j not in (i, (i+1) % n):
                assert dx*(qy-py)-dy*(qx-px) > 0
    return edges


def moment(edges, parameters, k, component):
    return s.expand(sum(e[component]*t**k
                        for e, t in zip(edges, parameters)))


def monomials_of_degree(variables, degree):
    def compositions(n, d):
        if n == 1:
            yield (d,)
        else:
            for a in range(d+1):
                for tail in compositions(n-1, d-a):
                    yield (a,) + tail
    return [s.prod(x**a for x, a in zip(variables, powers))
            for powers in compositions(len(variables), degree)]


def standard_monomials(gb, variables, max_degree):
    leading = [p.LM(order=gb.order).exponents for p in gb.polys]
    degrees, result = [], []
    for d in range(max_degree+1):
        current = []
        for m in monomials_of_degree(variables, d):
            exps = s.Poly(m, *variables).monoms()[0]
            if not any(all(a >= b for a, b in zip(exps, lm))
                       for lm in leading):
                current.append(m)
        degrees.append(len(current))
        result.extend(current)
    return degrees, result


def hexagon_test():
    P = [(0,0), (1,0), (2,1), (2,2), (-2,3), (-1,1)]
    edges = check_primitive_convex(P)
    assert edges == [(1,0), (1,1), (0,1), (-4,1), (1,-2), (1,-1)]
    # Unit edges imply a nontrivial lattice Minkowski summand would give such a subset.
    for k in range(1, 6):
        for I in combinations(range(6), k):
            assert any(sum(edges[i][j] for i in I) != 0 for j in (0,1))
    x, y, z = variables = s.symbols('x y z')
    parameters = [x+5*y-3*z, -x-y+2*z, x, y, z, s.Integer(0)]
    assert all(moment(edges, parameters, 1, j) == 0 for j in (0,1))
    full = [moment(edges, parameters, k, j)
            for k in range(2,6) for j in (0,1)]
    gb = s.groebner(full, *variables, order='grlex', domain=s.QQ)
    gb23 = s.groebner(full[:4], *variables, order='grlex', domain=s.QQ)
    # Equality of the k<=3 and the exhaustive k<=5 ideals.
    assert all(gb23.reduce(f)[1] == 0 for f in full)
    expected_gb = [
        y*z**3, z**4,
        3*x*z**2+5*y*z**2-z**3,
        y**3+y*z**2-z**3,
        2*y**2*z-z**3,
        5*x**2-7*x*z-5*y**2+5*y*z-z**2,
        5*x*y-3*x*z+10*y**2-15*y*z+6*z**2,
    ]
    expected = s.groebner(expected_gb, *variables, order='grlex', domain=s.QQ)
    assert gb == expected
    hs, standard = standard_monomials(gb, variables, 4)
    assert hs == [1,3,4,2,0]
    assert set(standard) == {s.Integer(1), x,y,z, x*z,y*y,y*z,z*z,y*z*z,z**3}
    assert gb.reduce(z**3)[1] == z**3
    assert all(gb.reduce(m)[1] == 0 for m in monomials_of_degree(variables, 4))
    q = x*x+x*y-2*x*z+y*y-2*y*z+z*z
    assert gb.reduce(q)[1] == 0
    bad = gb.reduce(s.expand(q.xreplace({x:x+z*z})))[1]
    assert bad == -(7*y*z*z+4*z**3)/3 and bad != 0

    # For the following calculations degree>=4 can be discarded: m^4=0 was proved above.
    def quotient_nf(expr, g):
        p = s.Poly(s.expand(expr), *variables, domain=g.domain)
        low = sum(coef*s.prod(v**a for v,a in zip(variables,powers))
                  for powers,coef in p.terms() if sum(powers) <= 3)
        return g.reduce(s.expand(low))[1]

    def multiply(u, v, g):
        return quotient_nf(u*v, g)

    def power(u, n, g):
        out = s.Integer(1)
        for _ in range(n):
            out = multiply(out, u, g)
        return out

    def evaluate(poly, images, g):
        out = 0
        for powers, coef in s.Poly(poly, *variables, domain=s.QQ).terms():
            term = coef
            for u, n in zip(images, powers):
                term = multiply(term, power(u, n, g), g)
            out += term
        return quotient_nf(out, g)

    a, b = s.symbols('a b')
    ga = s.groebner(full[:4], *variables, order='grlex',
                    domain=s.QQ.frac_field(a,b))
    # Coefficients genuinely depend on the SOURCE base, as well as external a,b.
    def H(T):
        return (1+x-y)*T + (a+z)*T*T + b*T**3
    edge_images = [quotient_nf(H(t), ga) for t in parameters]
    free_images = [quotient_nf(H(t), ga) for t in variables]
    assert all(evaluate(f, free_images, ga) == 0 for f in full)
    assert all(quotient_nf(sum(e[j]*power(u,k,ga)
                               for e,u in zip(edges,edge_images)), ga) == 0
               for k in range(1,6) for j in (0,1))

    # For a common cubic H, divided differences have this form, with common coefficients.
    c0, c1, c2 = s.symbols('c0 c1 c2')
    gu = s.groebner(full[:4], *variables, order='grlex',
                    domain=s.QQ.frac_field(c0,c1,c2))
    units = [quotient_nf(c0+c1*t+c2*t*t, gu) for t in parameters]
    for j in (0,1):
        positive = negative = s.Integer(1)
        for e, U in zip(edges, units):
            if e[j] > 0:
                positive = multiply(positive, power(U,e[j],gu), gu)
            elif e[j] < 0:
                negative = multiply(negative, power(U,-e[j],gu), gu)
        assert quotient_nf(positive-negative, gu) == 0
    def log_unit_without_constant(T):
        return (c1/c0*T + (c2/c0-c1*c1/(2*c0*c0))*T*T
                + (-c1*c2/c0**2+c1**3/(3*c0**3))*T**3)
    assert all(quotient_nf(sum(e[j]*log_unit_without_constant(t)
                               for e,t in zip(edges,parameters)),gu) == 0
               for j in (0,1))
    return hs


def ideal_intersection(first, second, variables):
    aux = s.Dummy('aux')
    g = s.groebner([aux*f for f in first]+[(1-aux)*f for f in second],
                   aux,*variables,order='lex',domain=s.QQ)
    return [p.as_expr() for p in g.polys if not p.as_expr().has(aux)]


def octagon_test():
    P = [(0,0),(1,0),(2,1),(2,2),(1,3),(0,3),(-1,2),(-1,1)]
    edges = check_primitive_convex(P)
    assert edges == [(1,0),(1,1),(0,1),(-1,1),(-1,0),(-1,-1),(0,-1),(1,-1)]
    A,B,C,D,E = variables = s.symbols('A B C D E')
    parameters = [A+2*B+C-E, -A-B+D+E, A,B,C,D,E,s.Integer(0)]
    assert all(moment(edges,parameters,1,j)==0 for j in (0,1))
    initial = [moment(edges,parameters,k,j) for k in (2,3) for j in (0,1)]
    gb = s.groebner(initial,*variables,order='grlex',domain=s.QQ)
    # Exhaustive, not an arbitrary truncation: n=8 implies moments through 7 suffice.
    assert all(gb.reduce(moment(edges,parameters,k,j))[1]==0
               for k in range(1,8) for j in (0,1))
    components = [[A-E,B], [A-D,B,C-E], [A,B-E,C]]
    radical = components[0]
    for component in components[1:]:
        radical = ideal_intersection(radical,component,variables)
    gr = s.groebner(radical,*variables,order='grlex',domain=s.QQ)
    expected = [
        A*A-A*D-A*E-B*D+D*E,
        A*B,
        A*C-A*E-B*E-C*E+E*E,
        B*B-B*E,
        B*C,
    ]
    assert gr == s.groebner(expected,*variables,order='grlex',domain=s.QQ)
    assert all(gr.reduce(f)[1]==0 for f in initial)
    # The intersection is radical; its generators are nilpotent modulo the full ideal.
    assert all(gb.reduce(s.expand(r**3))[1]==0 for r in expected)
    # These two inclusions prove that this intersection IS sqrt(I).
    nil = B*(B-E)
    nil2 = gb.reduce(s.expand(nil**2))[1]
    assert nil2 == B*B*D*D-B*B*D*E-B*D*D*E+B*D*E*E
    assert nil2 != 0 and gb.reduce(s.expand(nil**3))[1] == 0
    # In fact the entire nilradical has exact nilpotency index 3.
    assert all(gb.reduce(s.expand(s.prod(triple)))[1]==0
               for triple in combinations_with_replacement(expected,3))

    a = s.symbols('a')
    ga = s.groebner(initial,*variables,order='grlex',domain=s.QQ.frac_field(a))
    images = {x:x+a*x*x for x in variables}
    assert all(ga.reduce(s.expand(f.xreplace(images)))[1]==0 for f in initial)
    c0,c1=s.symbols('c0 c1')
    gu=s.groebner(initial,*variables,order='grlex',domain=s.QQ.frac_field(c0,c1))
    units=[c0+c1*t for t in parameters]
    for j in (0,1):
        pos=s.prod(U**e[j] for e,U in zip(edges,units) if e[j]>0)
        neg=s.prod(U**(-e[j]) for e,U in zip(edges,units) if e[j]<0)
        assert gu.reduce(s.expand(pos-neg))[1]==0


def main():
    hs=hexagon_test()
    print('PASS hexagon: Hilbert (1,3,4,2), length 10; m^3 != 0, m^4 = 0; bad shear detected.')
    print('PASS hexagon: source-dependent common cubic H; both unit-product and log identities.')
    octagon_test()
    print('PASS Q8: full moments, exact radical, N^2 != 0 and N^3 = 0; common H and unit identities.')


if __name__ == '__main__':
    main()
