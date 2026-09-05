#!/usr/bin/env python3
"""Exact inputs to the counterexample's full deformation-germ calculation.

Run: sage -python paper5/certificates/counterexample_deformation_germ.py
The analytic normal form also uses the separately proved smoothness of the
four selected bases. These finite checks do not establish that theorem.
"""
from itertools import combinations
import json
from pathlib import Path
from sage.all import QQ, PolynomialRing, matrix, vector


def rows_json(A):
    return [[str(x) for x in row] for row in A]


def main():
    checks = []

    def check(name, condition):
        assert condition, name
        checks.append(name)

    data = json.loads(Path(__file__).with_name(
        'a2_branch_selection_counterexample.json').read_text())
    h, e1, e2, e3 = matrix.identity(QQ, 4).rows()
    intersection = matrix.diagonal(QQ, [1, -1, -1, -1])
    boundary = [e1, h-e1-e2, e2, h-e2-e3, e3, h-e1-e3]
    roots = [-h+e1+e2+e3, -e1+e3, e2-e3]
    rows = []
    for face in data['nodes']:
        row = vector(QQ, 25)
        for i, vertex in enumerate(face):
            row[vertex] = (-1)**i
        rows.append(row)
    for face in data['dp6_faces']:
        for root in roots:
            row = vector(QQ, 25)
            for vertex, curve in zip(face, boundary):
                row[vertex] = root*intersection*curve
            rows.append(row)
    M = matrix(QQ, rows)
    K = M.left_kernel().basis_matrix()
    check('full matrix size', M.dimensions() == (36, 25))
    check('full relation rank', M.rank() == 20)
    check('full kernel dimension', K.nrows() == 16)
    check('kernel multiplication', K*M == 0)
    node_kernel = M[:30, :].left_kernel()
    check('node-only kernel dimension', node_kernel.dimension() == 11)
    image = K.matrix_from_columns(range(30, 36)).row_space()
    expected_image = matrix(QQ, [[1, 0, 0, 0, 0, 0],
                                [0, 1, 0, 0, 0, 0],
                                [0, 0, 1, 0, 0, -1],
                                [0, 0, 0, 1, 0, 0],
                                [0, 0, 0, 0, 1, 0]])
    check('cone image has dimension five', image.dimension() == 5)
    check('only cone relation is b0+b1=0', image == expected_image.row_space())
    embedded_kernels = []
    for rec in data['profiles']:
        ids = list(range(30))
        for i, branch in enumerate(rec['profile']):
            ids += [30+3*i] if branch == 'L' else [31+3*i, 32+3*i]
        A = M.matrix_from_rows(ids)
        check(rec['profile']+' recovers saved matrix', A == matrix(QQ, rec['matrix']))
        KS = A.left_kernel().basis_matrix()
        embedding = matrix(QQ, len(ids), 36)
        for i, j in enumerate(ids):
            embedding[i, j] = 1
        embedded = KS*embedding
        embedded_kernels.extend(embedded.rows())
        check(rec['profile']+' embeds in full kernel', embedded*M == 0)
        check(rec['profile']+' kernel dimension', KS.nrows() == rec['kernel_dimension'])
    check('selected kernels span the full kernel',
          matrix(QQ, embedded_kernels).row_space() == K.row_space())

    R = PolynomialRing(QQ, names=['x', 'y', 'u', 'v', 'b'])
    x, y, u, v, b = R.gens()
    J = R.ideal(x*y, x*b, u*v, u*b)
    primes = {'LL': R.ideal(y, v, b), 'LP': R.ideal(y, b, u),
              'PL': R.ideal(x, v, b), 'PP': R.ideal(x, u)}
    total = R.ideal(1)
    for name, prime in primes.items():
        check(name+' contains four quadratic equations', all(f in prime for f in J.gens()))
        total = total.intersection(prime)
    check('quadratic ideal is intersection of four coordinate primes', total == J)
    check('four independent quadratic generators', len(J.groebner_basis()) == 4)
    check('largest component dimension three', J.dimension() == 3)
    for a, bname in combinations(primes, 2):
        check(a+' and '+bname+' are incomparable',
              any(f not in primes[bname] for f in primes[a].gens()) and
              any(f not in primes[a] for f in primes[bname].gens()))
    intersections = []
    for count in range(1, 5):
        for names in combinations(primes, count):
            summed = R.ideal(0)
            for name in names:
                summed += primes[name]
            dimension = int(32+summed.dimension())
            check('/'.join(names)+' intersection is a coordinate space',
                  all(f.total_degree() == 1 for f in summed.groebner_basis()))
            intersections.append(dict(profiles=list(names), dimension=dimension,
                                      ideal=[str(f) for f in summed.groebner_basis()]))
    check('component dimensions', [r['dimension'] for r in intersections[:4]] == [34, 34, 34, 35])
    check('triple and quadruple intersections have dimension 32',
          all(r['dimension'] == 32 for r in intersections if len(r['profiles']) >= 3))
    series = J.hilbert_series()
    t = series.parent().gen()
    check('Hilbert series', series == (1+2*t-t*t-t**3)/(1-t)**3)

    # Identities proving straightening for arbitrary analytic coefficient
    # functions c,p,q with c a unit. The inverse identities multiply by c.
    S = PolynomialRing(QQ, names=['x', 'y', 'u', 'v', 'b', 'c', 'p', 'q'])
    x, y, u, v, b, c, p, q = S.gens()
    H = c*b+p*y+q*v
    B = c*b+p*y
    check('straightening first product', x*B == c*x*b+p*x*y)
    check('straightening second product', u*B == u*H-q*u*v)
    check('inverse first product after multiplying by unit', c*x*b == x*B-p*x*y)
    check('inverse second product', u*H == u*B+q*u*v)
    check('LL profile straightening after multiplying by unit', c*b == B-p*y)
    check('PL profile straightening', H-B == q*v)

    target = Path(__file__).with_suffix('.json')
    target.write_text(json.dumps(dict(
        full_matrix=rows_json(M), full_kernel_basis=rows_json(K),
        full_relation_rank=int(M.rank()), full_kernel_dimension=K.nrows(),
        cone_coordinate_order=['lambda0', 'a0', 'b0', 'lambda1', 'a1', 'b1'],
        cone_image_basis=rows_json(image.basis_matrix()), node_only_kernel_dimension=11,
        locally_trivial_dimension=21, predicted_tangent_dimension=37,
        predicted_smooth_factor_dimension=32,
        predicted_analytic_ring='C{z1,...,z32,x,y,u,v,b}/(xy,xb,uv,ub)',
        component_intersections=intersections, transverse_hilbert_series=str(series),
        checks=checks, passed=len(checks),
        scope='exact linear and polynomial inputs; analytic normal form additionally uses the branch-smoothness proof'
    ), indent=2)+'\n')
    print(f'PASS: {len(checks)} exact deformation-germ input checks')
    print('full matrix rank 20; kernel 16; cone image 5; node-only kernel 11')
    print('predicted tangent dimension 37; smooth factor 32; four component dimensions 34,34,34,35')
    print('quadratic ideal equals intersection of four coordinate primes; straightening identities pass')


if __name__ == '__main__':
    main()
