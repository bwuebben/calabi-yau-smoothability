"""Reconstruct four compact examples from vertices and cyclic hexagon markings.

Run with Python 3.10 or later. All calculations use integer determinants and
rational elimination from the standard library. No charge matrices are inputs.
"""

from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, combinations_with_replacement, permutations, product
import json
from math import gcd
from pathlib import Path


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def det3(a, b, c):
    return (a[0] * (b[1]*c[2] - b[2]*c[1])
            - a[1] * (b[0]*c[2] - b[2]*c[0])
            + a[2] * (b[0]*c[1] - b[1]*c[0]))


def normal(a, b, c):
    """Primitive up to a positive factor, with sign from row order."""
    return tuple((-1)**j * det3(*(x[:j] + x[j+1:] for x in (a, b, c)))
                 for j in range(4))


def normals(vertices):
    """Find every supporting facet through four affinely independent vertices."""
    answer = set()
    for o, a, b, c in combinations(vertices, 4):
        v = normal(sub(a, o), sub(b, o), sub(c, o))
        d = gcd(*v)
        if not d:
            continue
        v = tuple(x//d for x in v)
        offset = dot(v, o)
        heights = [dot(v, p) for p in vertices]
        if max(heights) == offset:
            v, offset = tuple(-x for x in v), -offset
        elif min(heights) != offset:
            continue
        assert offset == -1  # origin interior; facet at lattice distance one
        answer.add(v)
    return sorted(answer)


def lattice_points(vertices, constraints):
    count = Counter()
    by_codim_one = set()
    for p in product(*(range(min(a), max(a) + 1) for a in zip(*vertices))):
        h = [dot(n, p) for n in constraints]
        if min(h) < -1:
            continue
        n_facets = h.count(-1)
        count[n_facets] += 1
        if n_facets >= 2:
            by_codim_one.add(p)
    assert count[0] == 1
    return sum(count.values()), count[1], by_codim_one


def rref(rows, width):
    a = [list(map(Q, row)) for row in rows]
    assert all(len(row) == width for row in a)
    pivots = []
    for col in range(width):
        k = len(pivots)
        j = next((i for i in range(k, len(a)) if a[i][col]), None)
        if j is None:
            continue
        a[k], a[j] = a[j], a[k]
        a[k] = [x/a[k][col] for x in a[k]]
        for i, row in enumerate(a):
            if i != k and row[col]:
                a[i] = [x - row[col]*y for x, y in zip(row, a[k])]
        pivots.append(col)
    return a[:len(pivots)], pivots


def kernel(rows, width):
    a, pivots = rref(rows, width)
    answer = []
    for j in range(width):
        if j in pivots:
            continue
        v = [Q(i == j) for i in range(width)]
        for row, p in zip(a, pivots):
            v[p] = -row[j]
        answer.append(v)
    return answer


def relations(rows):
    return rref(kernel(zip(*rows), len(rows)), len(rows))[0]


def format_rows(rows):
    return [tuple(int(x) if x.denominator == 1 else str(x) for x in row)
            for row in rows]


def local_classes(number):
    data = json.loads(Path(__file__).with_name(f'record_{number}.json').read_text())
    v = list(map(tuple, data['vertices']))
    marks = data['hexagons']
    dual = normals(v)
    # Distance-one support in the opposite direction also detects redundant
    # listed points: every input point must actually occur as a polar facet.
    assert set(normals(dual)) == set(v)
    facets = [{i for i,p in enumerate(v) if dot(p, u) == -1} for u in dual]
    faces = sorted({tuple(sorted(a & b)) for a,b in combinations(facets, 2)
                    if len(a & b) >= 3})
    edges = {tuple(sorted(set(f) & a)) for f in faces for a in facets
             if len(set(f) & a) == 2}
    assert all(gcd(*sub(v[i], v[j])) == 1 for i,j in edges)
    assert {tuple(sorted(f)) for f in marks} == {f for f in faces if len(f) == 6}
    centers = []
    for f in marks:
        sums = [sum(v[i][j] for i in f) for j in range(4)]
        assert all(x % 6 == 0 for x in sums)
        centers.append(tuple(x//6 for x in sums))
    rays = v + centers
    assert len(set(rays)) == len(rays)
    prim, prim_interior, visible = lattice_points(v, dual)
    polar, polar_interior, _ = lattice_points(dual, v)
    assert visible == set(rays)
    nodes, roots, walls = [], [], []
    def row(coeffs):
        r = [coeffs.get(i, 0) for i in range(len(rays))]
        assert not any(dot(r, c) for c in zip(*rays))
        return r
    def is_unit_triangle(i, j, k):
        return gcd(*normal(rays[i], rays[j], rays[k])) == 1
    for f in faces:
        if len(f) == 3:
            assert is_unit_triangle(*f)
            continue
        adjacent = [u for u,facet in zip(dual, facets) if set(f) <= facet]
        assert len(adjacent) == 2 and gcd(*sub(*adjacent)) == 1
        if len(f) == 4:
            i,*others = f
            j, = [j for j in others if (i,j) not in edges]
            k,l = sorted(set(others) - {j})
            assert is_unit_triangle(i,j,k)
            nodes.append(row({i:1,j:1,k:-1,l:-1}))
        else:
            assert len(f) == 6
    for m, f in enumerate(marks):
        B = []
        center = len(v) + m
        for j, i in enumerate(f):
            n, p = f[j-1], f[(j+1) % 6]
            assert tuple(sorted((i,p))) in edges
            assert is_unit_triangle(center,i,p)
            B.append(row({n:1,p:1,i:-1,center:-1}))
        roots.append([sub(B[0],B[3]), sub(B[4],B[0]), sub(B[2],B[4])])
        walls.append(B)
    h11, h21 = prim - prim_interior - 5, polar - polar_interior - 5
    print(f'record {number}: f-vector = {(len(v),len(edges),len(faces),len(dual))}')
    print('  polygons',dict(sorted(Counter(map(len,faces)).items())),
          'lattice points/interior-to-facet',((prim,prim_interior),(polar,polar_interior)),
          'Hodge numbers',(h11,h21))
    return nodes, roots, walls, h21


def restrictions(k, n_nodes, n_hex):
    """Rows in Sym^2(K*) ordered by i<=j, with ordinary monomial coefficients."""
    degree = list(combinations_with_replacement(range(len(k)), 2))
    quadrics = []
    for m in range(n_hex):
        i = n_nodes + 3*m
        for j in (i+1, i+2):
            quadrics.append([k[a][i]*k[b][j] + int(a!=b)*k[b][i]*k[a][j]
                             for a,b in degree])
    return rref(quadrics, len(degree))[0]


def chosen(nodes, roots, keep):
    return nodes + [row for m in keep for row in roots[m]]


def forced_zeros(k, width):
    return [i for i in range(width) if all(v[i] == 0 for v in k)]


def control_248():
    n, roots, _, h = local_classes(248)
    k = relations(chosen(n,roots,(0,1)))
    assert h == 25 and len(n) == 12 and not relations(n) and len(k) == 1
    witness = [x/k[0][12] for x in k[0]]
    assert witness[12:] == [1,1,1,-1,-1,-1]
    print('  unique relation',format_rows([witness])[0])
    assert restrictions(k,len(n),2) == [[1]]
    for profile in product('LP', repeat=2):
        rows = n + [r for m,p in enumerate(profile)
                    for r in (roots[m][:1] if p=='L' else roots[m][1:])]
        assert not relations(rows)
    print('  every selected kernel = 0; transverse Hilbert function (1,1)')


def control_528():
    n, roots, walls, h = local_classes(528)
    k = relations(chosen(n,roots,(0,1)))
    assert len(n) == 30 and h == 18 and len(relations(n)) == 10 and len(k) == 14
    projected = [v[len(n):] for v in k]
    eq = rref(kernel(projected,6),6)[0]
    assert eq == [[0,1,0,0,0,0],[0,0,1,1,0,0]]
    print('  cone equations',format_rows(eq))
    expected = {'LL':(29,[17,32,33]), 'LP':(31,[]), 'PL':(29,[17,30]),
                'PP':(30,[30,31,32]), '00':(28,[17])}
    for profile,(dim,zeros) in expected.items():
        rows = n + [walls[m][i] for m,p in enumerate(profile)
                    for i in {'L':(0,3),'P':(0,2,4),'0':()}[p]]
        rel = relations(rows)
        assert (h + len(rel),forced_zeros(rel,len(rows))) == (dim,zeros)
        print(' ',profile,dim,'forced incidence nodes',zeros)
    print('  maximal component dimensions 31,29; intersection dimension 28')


def cluster(roots, indices):
    """Three pure-profile relations sum to the mixed relation on the triple."""
    k = relations(chosen([], roots, indices[:3]))
    assert len(k) == 1
    target = k[0] + [0,0,0]
    vectors = []
    for line in range(3):
        allowed = [3*m+j for m in range(4)
                   for j in ((0,) if m==line else (1,2))]
        full = chosen([],roots,indices)
        selected = relations([full[j] for j in allowed])
        assert len(selected) == 1
        v = [0]*12
        for i,a in zip(allowed,selected[0]):
            v[i] = a
        vectors.append(v)
    # Solve for coefficients, rejecting a merely set-theoretic decomposition.
    relation = relations(vectors+[target])
    assert len(relation) == 1 and relation[0][-1]
    coeff = [-a/relation[0][-1] for a in relation[0][:-1]]
    summands = [[a*x for x in v] for a,v in zip(coeff,vectors)]
    assert list(map(sum,zip(*summands))) == target
    print('  cluster',tuple(indices),'triple tangent',format_rows([target])[0])
    print('    LPPP, PLPP, PPLP summands',format_rows(summands))
    return target


def symmetric_model():
    # These functionals give a unimodular map to the sum-zero lattice.
    original = json.loads(Path(__file__).with_name('record_575.json').read_text())
    axes = [(0,1,-1,-1),(1,0,1,1),(-1,-1,0,-1),(-1,-1,-1,0),(1,1,1,1)]
    assert all(sum(c)==0 for c in zip(*axes)) and abs(dot(normal(*axes[:3]),axes[3])) == 1
    census_vertices = [tuple(dot(p,a) for a in axes) for p in original['vertices']]
    vertices = sorted(set(permutations((1,1,-1,-1,0))))
    assert sorted(census_vertices) == vertices
    lookup = {p:i for i,p in enumerate(vertices)}
    cyclic = [(-1,0,1),(0,-1,1),(1,-1,0),(1,0,-1),(0,1,-1),(-1,1,0)]
    roots = {}
    for i,j in permutations(range(5),2):
        free = [k for k in range(5) if k not in (i,j)]
        mark = []
        for a in cyclic:
            p = [0]*5
            p[i],p[j] = 1,-1
            for k,val in zip(free,a):
                p[k] = val
            mark.append(lookup[tuple(p)])
        walls = [[int(q in (mark[t-1],mark[(t+1)%6])) - int(q==mark[t])
                  for q in range(len(vertices))] for t in range(6)]
        roots[(i+1,j+1)] = [sub(walls[0],walls[3]),sub(walls[4],walls[0]),sub(walls[2],walls[4])]
    named = []
    for f in original['hexagons']:
        plus, = [i+1 for i in range(5) if all(census_vertices[p][i]==1 for p in f)]
        minus, = [i+1 for i in range(5) if all(census_vertices[p][i]==-1 for p in f)]
        named.append((plus,minus))
    a = [(5,4),(5,3),(5,1)]
    b = [(5,4),(2,4),(1,4)]
    c = [(2,4),(2,3),(2,1)]
    assert a == [named[i] for i in (0,1,2)]
    assert b == [named[i] for i in (0,4,5)]
    assert c == [named[i] for i in (4,7,10)]
    print('  symmetric cube slice, triples in cyclic markings of the manuscript')
    expected = [(1,0,-1,-1,0,1,-1,-1,-1), (1,1,0,1,0,1,-1,0,-1),
                (1,-1,0,-1,1,0,1,-1,0)]
    for g,extra,vector in zip((a,b,c),((5,2),(3,4),(2,5)),expected):
        r = relations(chosen([],roots,g))
        assert format_rows(r) == [vector]
        print('   ',g,vector)
        cluster(roots,g+[extra])
    for keep,rank in [(a+b[1:],13),(a+c,16)]:
        rows = chosen([],roots,keep)
        assert len(rows)-len(relations(rows)) == rank
        print('    union',keep,'curve matrix rank',rank)


def new_examples():
    n, roots, _, h = local_classes(575)
    assert not n and len(roots) == 20
    for f in [(0,1,2,3),(0,4,5,6),(4,7,10,11)]:
        cluster(roots,f)
    for keep, expected in [((0,1,2,4,5), [[1,0,0],[0,1,0],[0,0,1]]),
                            ((0,1,2,4,7,10), [[1,0,0],[0,0,1]]),
                            ((0,1,2,5,6), [[1,0,1],[0,1,0]])]:
        rows = chosen(n,roots,keep)
        k = relations(rows)
        assert len(k) == 2
        quad = restrictions(k,0,len(keep))
        assert quad == expected
        for p in product('LP',repeat=len(keep)):
            selected = [r for m,s in zip(keep,p)
                        for r in (roots[m][:1] if s=='L' else roots[m][1:])]
            assert not relations(selected)
        print('  retained',keep,'relation basis',format_rows(k))
        print('    restricted quadrics (x^2,xy,y^2)',format_rows(quad))
        print('    selected dimensions all zero; unrestricted tangent',h+len(k))
    symmetric_model()
    print('  overlapping completed clusters: transverse Hilbert function (1,2)')
    print('  disjoint completed clusters: transverse Hilbert function (1,2,1)')
    print('  retained (0,1,2,5,6): predicted transverse Hilbert function (1,2,1)')


if __name__ == '__main__':
    control_248()
    control_528()
    new_examples()
