"""Exact fan and moment calculations for the compact cubic example."""

from math import gcd

from sympy import Matrix, QQ, expand, groebner, symbols


polygon = [(0, -1), (1, 0), (0, 1), (-1, 1), (-2, 0), (-1, -1)]


def cyclic_pairs(vertices):
    return zip(vertices, vertices[1:] + vertices[:1])


def determinant(a, b, c):
    return Matrix([a, b, c]).det()


def check_fan():
    upper = [(*point, 1) for point in polygon]
    apex = (0, 0, -1)
    midpoint = (-1, 0, 0)
    centre, second = (0, 0, 1), (-1, 0, 1)
    side_normals = [(2, -2, -1), (2, 2, -1), (0, 2, -1),
                    (-1, 1, -1), (-1, -1, -1), (0, -2, -1)]
    determinants = []
    for (a, b), normal in zip(cyclic_pairs(upper), side_normals):
        on_facet = []
        for vertex in [*upper, apex]:
            value = sum(x*y for x, y in zip(normal, vertex))
            assert value <= 1
            if value == 1:
                on_facet.append(vertex)
        assert set(on_facet) == {a, b, apex}
        determinants.append(abs(determinant(a, b, apex)))
        # The two subdivisions bisect the edges opposite the same boundary ray.
        if determinants[-1] == 2:
            other = next(v for v in (a, b) if v != (-2, 0, 1))
            for end, middle in ((apex, midpoint), (centre, second)):
                assert abs(determinant(other, end, middle)) == 1
                assert abs(determinant(other, middle, (-2, 0, 1))) == 1
        else:
            assert abs(determinant(a, b, centre)) == 1
    assert determinants == [1, 1, 1, 2, 2, 1]
    assert determinant(upper[1], upper[2], upper[3]) == 1


def check_moments():
    edges = [(c-a, d-b) for (a, b), (c, d) in cyclic_pairs(polygon)]
    assert all(gcd(dx, dy) == 1 for dx, dy in edges)
    u, v, s = symbols("u v s")
    heights = (u + 3*s, v + s, 4*s, u + s, v + 3*s, 0)

    def moment(k):
        return [expand(sum(d[j]*t**k for d, t in zip(edges, heights)))
                for j in (0, 1)]

    assert moment(1) == [0, 0]
    quadrics = groebner(moment(2), u, v, s, domain=QQ)
    assert list(quadrics) == [u*s, v*s]
    assert [quadrics.reduce(f)[1] for f in moment(3)] == [-12*s**3, 0]
    # For six heights the characteristic-polynomial recurrence needs k <= 5.
    equations = [f for k in range(2, len(heights)) for f in moment(k)]
    basis = groebner(equations, u, v, s, domain=QQ)
    assert list(basis) == [u*s, v*s, s**3]
    assert basis.reduce(s**2)[1] == s**2


if __name__ == "__main__":
    check_fan()
    check_moments()
    print("Unimodular subdivisions verified. The moment ideal is (us, vs, s³) and s² survives.")
