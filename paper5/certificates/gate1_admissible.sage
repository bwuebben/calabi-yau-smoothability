"""Admissibility refuter for the local dP7 crosscut in Theorem A.

The retired crosscut search asks only whether the crosscut is a Minkowski sum
of polyhedra whose cones are regular.  That is not enough for an
Ilten--Vollmert deformation: Definition 2.2 also requires that, for every
character, at most one summand face contain no lattice point.

For a requested degree this script:

* computes the extreme rays of the Minkowski-summand cone;
* keeps the rays whose summand cones are regular;
* enumerates vertex representations of the all-ones dilation by those rays;
* rescales the summands by the representation coefficients and rechecks
  regularity (regularity is not invariant under rational scaling);
* adds the point translation needed to recover the actual crosscut; and
* tests Definition 2.2 on every vertex chamber of the common normal fan.

The translation convention is canonical: all positive-dimensional summands
start at the lattice point zero and the residual translation is a point
summand.  Thus a PASS is a genuine admissible decomposition.  A failure is a
refuter result for this convention, not by itself a proof that no redistribution
of the translations can work.

Run from paper5/:
    sage certificates/gate1_admissible.sage
"""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper4", "certificates"))
from batyrev_global import dot
from examples import V_19

V = [tuple(v) for v in V_19]
Vv = [vector(ZZ, v) for v in V]
PENT = [14, 15, 16, 17, 18]
SP = Polyhedron(rays=[list(Vv[i]) for i in PENT], base_ring=QQ)
LAT = matrix(ZZ, [list(Vv[i]) for i in PENT]).row_space().intersection(ZZ**4)
B = matrix(ZZ, LAT.basis())


def co(x):
    """Coordinates in the saturated rank-three lattice of the germ."""
    return vector(QQ, B.transpose().solve_right(vector(QQ, x)))


def summand_poly(edges, dilations):
    """Polygon obtained by walking the cyclic edge dilations."""
    points = [vector(QQ, [0, 0, 0])]
    current = vector(QQ, [0, 0, 0])
    for edge, dilation in zip(edges, dilations):
        current += dilation * edge
        points.append(current)
    return Polyhedron(vertices=[list(point) for point in points], base_ring=QQ)


def cone_smooth_at_height1(poly):
    """Whether the cone over a rational polygon at height one is regular."""
    generators = []
    for vertex in poly.vertices_list():
        generator = vector(QQ, list(vertex) + [1])
        generator *= lcm([entry.denominator() for entry in generator])
        content = gcd([ZZ(entry) for entry in generator])
        generators.append(vector(ZZ, [ZZ(entry) / content for entry in generator]))
    generators = list({tuple(generator) for generator in generators})
    matrix_generators = matrix(ZZ, generators)
    rank = matrix_generators.rank()
    return len(generators) == rank and gcd(matrix_generators.minors(rank)) == 1


def lattice_face(poly, u):
    """The u-minimal face and whether it contains a lattice point."""
    vertices = [vector(QQ, v) for v in poly.vertices_list()]
    values = [u.dot_product(v) for v in vertices]
    minimum = min(values)
    face = Polyhedron(vertices=[list(v) for v, a in zip(vertices, values)
                                if a == minimum], base_ring=QQ)
    return face, bool(face.integral_points())


def minkowski_sum(polys):
    out = Polyhedron(vertices=[[0, 0, 0]], base_ring=QQ)
    for poly in polys:
        out = out + poly
    return out


def vertex_test_directions(poly):
    """Integral characters in the interiors of all vertex normal cones."""
    inequalities = [vector(QQ, q) for q in poly.inequalities_list()]
    out = []
    for v in [vector(QQ, x) for x in poly.vertices_list()]:
        active = [q for q in inequalities
                  if q[0] + vector(QQ, q[1:]).dot_product(v) == 0]
        u = sum((vector(QQ, q[1:]) for q in active), vector(QQ, [0, 0, 0]))
        den = lcm([a.denominator() for a in u])
        out.append(vector(ZZ, [ZZ(den * a) for a in u]))
    return out


def data_for_degree(R0):
    R = vector(ZZ, R0)
    pairings = [dot(R, V[i]) for i in PENT]
    if min(pairings) <= 0:
        raise ValueError("this diagnostic expects a strictly positive crosscut")

    crosscut = SP & Polyhedron(eqns=[[-1] + list(R)], base_ring=QQ)
    poly = Polyhedron(vertices=[list(co(v)) for v in crosscut.vertices_list()],
                      base_ring=QQ)
    vertices = [vector(QQ, v) for v in poly.vertices_list()]
    adjacency = {}
    for edge in poly.faces(1):
        a, b = [tuple(vector(QQ, x.vector())) for x in edge.vertices()]
        adjacency.setdefault(a, []).append(b)
        adjacency.setdefault(b, []).append(a)
    start = tuple(vertices[0])
    cycle = [start]
    previous, current = None, start
    while len(cycle) < len(vertices):
        nxt = [x for x in adjacency[current] if x != previous][0]
        cycle.append(nxt)
        previous, current = current, nxt
    edges = [vector(QQ, cycle[(i + 1) % len(cycle)]) - vector(QQ, cycle[i])
             for i in range(len(cycle))]

    rows = [[edges[j][i] for j in range(len(edges))] for i in range(3)]
    kernel = matrix(QQ, rows).right_kernel().basis()
    cone = (Polyhedron(rays=[list(b) for b in kernel] +
                              [list(-b) for b in kernel], base_ring=QQ) &
            Polyhedron(ieqs=[[0] + [1 if k == j else 0
                                    for k in range(len(edges))]
                              for j in range(len(edges))], base_ring=QQ))
    rays = [vector(QQ, r) for r in cone.rays()]
    good = [(r, summand_poly(edges, r)) for r in rays
            if cone_smooth_at_height1(summand_poly(edges, r))]
    return pairings, poly, cycle, edges, good


def check_degree(R0, verbose=True):
    pairings, poly, cycle, edges, good = data_for_degree(R0)
    if verbose:
        print(f"\nR = {tuple(R0)}, pentagon pairings = {tuple(pairings)}")
        print(f"  regular extreme summands: {len(good)}")
    if not good:
        return False

    G = matrix(QQ, [list(r) for r, _ in good]).transpose()
    target = vector(QQ, [1] * len(edges))
    coefficient_poly = Polyhedron(
        eqns=[[-target[i]] + list(G.row(i)) for i in range(G.nrows())],
        ieqs=[[0] + [1 if i == j else 0 for i in range(G.ncols())]
              for j in range(G.ncols())], base_ring=QQ)
    if coefficient_poly.is_empty():
        if verbose:
            print("  the regular rays do not span the all-ones dilation")
        return False

    directions = vertex_test_directions(poly)
    passed = False
    for number, coeffs0 in enumerate(coefficient_poly.vertices_list(), start=1):
        coeffs = vector(QQ, coeffs0)
        summands = []
        regular = True
        for coefficient, (_, summand) in zip(coeffs, good):
            if coefficient == 0:
                continue
            scaled = coefficient * summand
            regular = regular and cone_smooth_at_height1(scaled)
            summands.append(scaled)
        edge_sum = minkowski_sum(summands)
        translation = vector(QQ, cycle[0])
        point = Polyhedron(vertices=[list(translation)], base_ring=QQ)
        summands_with_point = summands + [point]
        exact = minkowski_sum(summands_with_point) == poly
        violations = 0
        for u in directions:
            missing = sum(not lattice_face(summand, u)[1]
                          for summand in summands_with_point)
            violations += missing > 1
        admissible = regular and exact and violations == 0
        if verbose:
            print(f"  representation {number}: coefficients={tuple(coeffs)}, "
                  f"regular={regular}, exact={exact}, "
                  f"admissibility violations={violations}")
        passed = passed or admissible
    return passed


if __name__ == "__main__":
    tests = [
        (-3, -8, 5, 0),       # pairings (3,3,3,3,3): false positive without D2
        (-2, -3, 2, 1),       # pairings (1,2,1,2,3): published off-line claim
        (-3, -3, -2, -1),     # pairings (5,4,4,2,2): below the old area bound
    ]
    results = [check_degree(R) for R in tests]
    print("\nSummary:", list(zip(tests, results)))
