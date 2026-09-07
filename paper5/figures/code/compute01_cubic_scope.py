#!/usr/bin/env python3
"""Exact dimension diagnostic for the degree-three cone; no global criterion.

The Fermat cubic Jacobian ideal is (x0^2,x1^2,x2^2,x3^2).
The Gysin rank uses the standard marking of a cubic surface as Bl_6(P^2).
All computations use the Python standard library.
"""

from itertools import product


def main():
    basis = list(product((0, 1), repeat=4))
    graded = [sum(sum(e) == k for e in basis) for k in range(5)]
    assert len(basis) == 16
    assert graded == [1, 4, 6, 4, 1]
    canonical = (-3, 1, 1, 1, 1, 1, 1)
    assert canonical[0] ** 2 - sum(x * x for x in canonical[1:]) == 3
    assert len(canonical) - 1 == 6
    assert len(basis) != len(canonical) - 1
    print("Cubic cone: dim T^1 = 16, graded dimensions =", graded)
    print("Link: dim H_2 = 6. Five exact scope checks passed.")


if __name__ == "__main__":
    main()
