"""Verify the generators and products used in the octagon example in Section 5."""

from itertools import combinations_with_replacement

from sympy import QQ, expand, groebner, prod, symbols


def main():
    a, b, c, d, e = variables = symbols("a b c d e")
    heights = (a + 2*b + c - e, -a - b + d + e, a, b, c, d, e, 0)
    # The two triples on each side of the x and y moment equations.
    edge_rows = (((0, 1, 7), (3, 4, 5)), ((1, 2, 3), (5, 6, 7)))
    moments = []
    for positive, negative in edge_rows:
        for degree in (2, 3):
            moments.append(sum(heights[i]**degree for i in positive)
                           - sum(heights[i]**degree for i in negative))
    basis = groebner(moments, *variables, order="grlex", domain=QQ)

    def normal_form(polynomial):
        return basis.reduce(expand(polynomial))[1]

    u, v, w = a*b, b*(b-e), b*c
    f, h = (a-e)*(a-d), (a-e)*(c-e)
    assert normal_form(f - d*b + u + v) == 0
    assert normal_form(h - e*b + 2*u + v + 2*w) == 0
    for relation in (u*(u+v), v*(u+v), u*w, v*w, v**3, w**3):
        assert normal_form(relation) == 0
    for triple in combinations_with_replacement((u, v, w), 3):
        assert normal_form(prod(triple)) == 0
    square = normal_form(v**2)
    assert square == expand(b*d*(b-e)*(d-e)) and square != 0
    print("Octagon: the three generators have zero triple products; v^2 =", square)


if __name__ == "__main__":
    main()
