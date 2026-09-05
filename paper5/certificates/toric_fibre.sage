"""Exact toric description of the general ambient fibre for X_9.

The positive slice is the whole tail fan. Splitting the negative slice
therefore leaves only two nontrivial slices. The original T^3-action extends
to T^4 by Suess, Canonical divisors on T-varieties, Corollary 1.9:
https://arxiv.org/pdf/0811.0626

This certificate reconstructs the actual displaced coefficients, including
their markings, and checks them against an explicit complete four-dimensional
fan. An empty coefficient is an omitted point in one chart's locus; it does
not remove a cell supplied by another chart from the global subdivision.

The fan also verifies that -K is ample Cartier, h^0(-K)=162, and the singular
locus is one fixed point. This checks the finite data of the stated family;
the separate slice_admissible.sage verifies its admissibility axioms.

Run from paper5/:
    sage certificates/toric_fibre.sage
or:
    sage -python certificates/toric_fibre.sage

For a proposed copy outside the repository:
    sage -python toric_fibre.sage --project-root /path/to/cy_smoothing
"""

import argparse
from collections import Counter
import itertools
import json
from pathlib import Path
import sys

from sage.all import Cone, Fan, Polyhedron, QQ, ToricLattice, ZZ
from sage.all import binomial, gcd, lcm, matrix, vector


DECOMPOSED_FACET = 11
EDGE_ORDER = (
    (2, 5), (2, 6), (2, 7), (3, 6), (3, 7), (3, 8),
    (4, 5), (4, 7), (4, 8), (5, 6), (7, 8),
)
DILATIONS = (0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0)

# Coordinates are (t, x, z, w), with the last three coordinates in ker(R).
# These printed data are checked against the reconstructed fan, not used to
# construct its coefficients or its cones.
EXPECTED_RAYS = (
    (-1, -3, -1, 0), (-1, -3, 0, -1), (-1, -1, 0, 0),
    (-1, 0, 0, 0), (0, 1, 0, 0), (1, -1, 0, 0), (1, 3, 1, 1),
)
EXPECTED_MAXIMAL_CONES = {
    (1, 4, 5, 6), (1, 3, 4, 6), (0, 4, 5, 6), (0, 3, 4, 6),
    (0, 1, 4, 5), (0, 1, 3, 4), (0, 1, 2, 5, 6),
    (1, 2, 3, 6), (0, 2, 3, 6), (0, 1, 2, 3),
}


def primitive(ray):
    """Canonical primitive generator of a rational ray."""
    ray = vector(QQ, ray)
    integral = vector(ZZ, lcm([x.denominator() for x in ray]) * ray)
    return tuple(x / gcd(integral) for x in integral)


def poly_key(poly):
    """Exact representation, independent of the order of vertices and rays."""
    if poly.is_empty():
        return ("empty",)
    if poly.lines():
        raise AssertionError("A coefficient polyhedron contains a line")
    return (
        tuple(sorted(tuple(v) for v in poly.vertices())),
        tuple(sorted(primitive(r) for r in poly.rays())),
    )


def face_closure(polys):
    """The nonempty polyhedral complex determined by its maximal cells."""
    result = set()
    for poly in polys:
        if not poly.is_empty():
            for dim in range(poly.dim() + 1):
                result.update(poly_key(f.as_polyhedron()) for f in poly.faces(dim))
    return result


def marked_maxima(records):
    """Marks belong to full-dimensional cells, independent of chart repeats."""
    result = {}
    for poly, marked in records:
        if poly.dim() == 3:
            key = poly_key(poly)
            result[key] = result.get(key, False) or marked
    return result


def cone_indices(cone, rays):
    index = {ray: i for i, ray in enumerate(rays)}
    return tuple(sorted(index[primitive(r)] for r in cone.rays()))


def upgraded_slice(cone, level):
    """Slice a four-dimensional fan cone by its first coordinate."""
    poly = Polyhedron(vertices=[[0, 0, 0, 0]],
                      rays=[list(r) for r in cone.rays()], base_ring=QQ)
    section = poly & Polyhedron(eqns=[[-level, 1, 0, 0, 0]], base_ring=QQ)
    if section.is_empty():
        return Polyhedron(ambient_dim=3, base_ring=QQ)
    return Polyhedron(vertices=[list(v)[1:] for v in section.vertices()],
                      rays=[list(r)[1:] for r in section.rays()], base_ring=QQ)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path,
                        default=Path(__file__).resolve().parents[2])
    root = parser.parse_args().project_root.resolve()
    candidate = root / "paper5" / "certificates" / "v09_candidate.json"
    if not candidate.is_file() or not (root / "src" / "batyrev_global.py").is_file():
        parser.error("--project-root must name the cy_smoothing directory")
    sys.path.insert(0, str(root / "src"))
    from batyrev_global import facets, two_faces

    checks = []

    def check(label, condition):
        if not condition:
            raise AssertionError(label)
        checks.append(label)
        print("  [ok] " + label)

    with candidate.open() as stream:
        data = json.load(stream)
    # batyrev_global uses fractions.Fraction and therefore receives Python ints.
    vertices_int = [tuple(int(x) for x in v) for v in data["V"]]
    vertices = [vector(ZZ, v) for v in vertices_int]
    facet_data = facets(vertices_int)
    facet_vertices = facet_data[DECOMPOSED_FACET][2]
    degree = vector(ZZ, [-x for x in facet_data[DECOMPOSED_FACET][0]])
    splitting = vertices[1]
    levels = [degree.dot_product(v) for v in vertices]
    check("R=(0,1,-1,-1), with levels (0,1,-1,...,-1)",
          tuple(degree) == (0, 1, -1, -1)
          and levels == [0, 1, -1, -1, -1, -1, -1, -1, -1])

    # e_x, e_y+e_z, e_y+e_w are a saturated basis of ker(R).
    def kernel_coords(value):
        value = vector(QQ, value)
        if degree.dot_product(value) != 0:
            raise AssertionError("A vector is outside ker(R)")
        return vector(QQ, (value[0], value[2], value[3]))

    def project(value):
        value = vector(QQ, value)
        return kernel_coords(value - degree.dot_product(value) * splitting)

    whole = frozenset(range(len(vertices)))
    seen, pending = {whole}, [whole]
    while pending:
        face = pending.pop()
        for _, _, facet in facet_data:
            intersection = face & facet
            if intersection and intersection not in seen:
                seen.add(intersection)
                pending.append(intersection)
    faces = sorted(seen - {whole}, key=lambda face: (len(face), sorted(face)))
    check("72 nonzero face cones; {v0,v1} is an edge",
          len(faces) == 72 and frozenset((0, 1)) in faces)

    # Derive the facet edge graph independently of the printed dilation order.
    boundary_faces = [face for face, adjacent in two_faces(vertices_int, facet_data)
                      if DECOMPOSED_FACET in adjacent]
    edges = sorted({tuple(sorted(a & b))
                    for a, b in itertools.combinations(boundary_faces, 2)
                    if len(a & b) == 2})
    check("the eleven facet edges agree with the printed dilation order",
          tuple(edges) == EDGE_ORDER)
    adjacency = {i: [] for i in facet_vertices}
    for (a, b), dilation in zip(edges, DILATIONS):
        adjacency[a].append((b, dilation))
        adjacency[b].append((a, dilation))
    base = min(facet_vertices)
    phi0, pending = {base: vector(ZZ, 4)}, [base]
    while pending:
        a = pending.pop()
        for b, dilation in adjacency[a]:
            value = phi0[a] + dilation * (vertices[b] - vertices[a])
            if b in phi0:
                if value != phi0[b]:
                    raise AssertionError("Inconsistent edge displacement")
            else:
                phi0[b] = value
                pending.append(b)
    check("displacements are defined consistently on all seven facet vertices",
          set(phi0) == set(facet_vertices))
    phi1 = {i: vertices[i] - vertices[base] - phi0[i] for i in phi0}
    a_vertex = {i: kernel_coords(phi0[i] + vertices[base] + splitting) for i in phi0}
    b_vertex = {i: kernel_coords(phi1[i]) for i in phi0}

    level_planes = {level: Polyhedron(eqns=[[-level] + list(degree)], base_ring=QQ)
                    for level in (-1, 0, 1)}
    cells = {}
    for face in faces:
        cone = Polyhedron(rays=[list(vertices[i]) for i in sorted(face)], base_ring=QQ)
        zero = cone & level_planes[0]
        tail = Polyhedron(vertices=[[0, 0, 0]],
                          rays=[list(kernel_coords(r)) for r in zero.rays()], base_ring=QQ)

        def original_section(level):
            section = cone & level_planes[level]
            if section.is_empty():
                return Polyhedron(ambient_dim=3, base_ring=QQ)
            return Polyhedron(vertices=[list(project(v)) for v in section.vertices()],
                              rays=[list(kernel_coords(r)) for r in section.rays()],
                              base_ring=QQ)

        positive, negative = original_section(1), original_section(-1)
        on_facet = sorted(face & facet_vertices)
        if negative.is_empty():
            a_poly = b_poly = Polyhedron(ambient_dim=3, base_ring=QQ)
        else:
            if not on_facet:
                raise AssertionError("A nonempty negative cell misses the facet")
            a_poly = Polyhedron(vertices=[list(a_vertex[i]) for i in on_facet],
                                rays=tail.rays_list(), base_ring=QQ)
            b_poly = Polyhedron(vertices=[list(b_vertex[i]) for i in on_facet],
                                rays=tail.rays_list(), base_ring=QQ)
        cells[face] = dict(tail=tail, positive=positive, negative=negative,
                           A=a_poly, B=b_poly,
                           marked=not positive.is_empty() and not negative.is_empty())

    positive_cells = [c["positive"] for c in cells.values() if not c["positive"].is_empty()]
    negative_cells = [c for c in cells.values() if not c["negative"].is_empty()]
    tails = {poly_key(c["tail"]): c["tail"] for c in cells.values()}
    check("every nonempty positive coefficient equals its tail",
          all(c["positive"] == c["tail"] for c in cells.values()
              if not c["positive"].is_empty()))
    check("the positive slice equals all 29 distinct tail cones",
          len(positive_cells) == len(tails) == 29
          and set(map(poly_key, positive_cells)) == set(tails))
    check("A+B equals the original negative coefficient on all 69 nonempty cells",
          len(negative_cells) == 69
          and all(c["A"] + c["B"] == c["negative"] for c in negative_cells))

    tail_fan = Fan([Cone(p.rays_list(), lattice=ToricLattice(3)) for p in tails.values()],
                   discard_faces=True, check=True)
    check("the tail fan is complete, with seven rays and eight maximal cones",
          tail_fan.is_complete() and len(tail_fan.rays()) == 7
          and len(tail_fan.generating_cones()) == 8)
    tail_rays = tuple(sorted(primitive(r) for r in tail_fan.rays()))
    print("TAIL RAYS:", tail_rays)
    print("TAIL MAXIMAL CONES:", sorted(cone_indices(c, tail_rays)
                                        for c in tail_fan.generating_cones()))

    rows = [(i, cells[face]) for i, (_, _, face) in enumerate(facet_data)]
    complete = [i for i, c in rows if c["marked"]]
    affine = [i for i, c in rows if not c["marked"]]
    check("eight complete-locus facets and four affine-locus facets",
          complete == [0, 1, 2, 3, 4, 8, 9, 10] and affine == [5, 6, 7, 11])

    lattice4 = ToricLattice(4)

    def lifted(poly, sign):
        return [[sign] + list(v) for v in poly.vertices()]

    rebuilt, candidates = [], []
    for i, c in rows:
        tail_rays4 = [[0] + list(r) for r in c["tail"].rays()]
        old = Cone(lifted(c["positive"], 1) + lifted(c["negative"], -1) + tail_rays4,
                   lattice=lattice4)
        target = Cone([[degree.dot_product(vertices[j])] + list(project(vertices[j]))
                       for j in facet_data[i][2]], lattice=lattice4)
        rebuilt.append(set(map(primitive, old.rays())) == set(map(primitive, target.rays())))

        a_rays, b_rays = lifted(c["A"], 1), lifted(c["B"], -1)
        if c["marked"]:
            candidates.append(Cone(a_rays + b_rays + tail_rays4, lattice=lattice4))
        else:
            # Affine-locus pieces remain separate. Combining them would impose
            # a complete-locus construction on a chart where it does not apply.
            candidates.extend([Cone(a_rays + tail_rays4, lattice=lattice4),
                               Cone(b_rays + tail_rays4, lattice=lattice4)])
    check("the original two slices reconstruct all twelve special-fibre cones",
          len(rebuilt) == 12 and all(rebuilt))
    check("all reconstructed general-fibre cones are strongly convex",
          all(c.is_strictly_convex() for c in candidates))
    general_fan = Fan(candidates, discard_faces=True, check=True)
    check("the reconstructed general fan is complete", general_fan.is_complete())
    rays = tuple(sorted(primitive(r) for r in general_fan.rays()))
    maximal_cones = {cone_indices(c, rays) for c in general_fan.generating_cones()}
    check("the reconstructed seven rays and ten maximal cones match the printed data",
          rays == EXPECTED_RAYS and maximal_cones == EXPECTED_MAXIMAL_CONES)
    print("GENERAL RAYS:")
    for i, ray in enumerate(rays):
        print("  r{} = {}".format(i, ray))
    print("GENERAL MAXIMAL CONES:", sorted(maximal_cones))

    # Equality of maximal marked subdivisions identifies the T-variety.
    # Check all lower-dimensional input coefficients as faces as well.
    expected_counts = {"A": (8, 8), "B": (10, 8), "tail": (8, 8)}
    for name, level in (("A", 1), ("B", -1), ("tail", 0)):
        source = [(c[name], c["marked"]) for _, c in rows]
        downgraded = []
        for cone in general_fan.generating_cones():
            marked = (any(r[0] > 0 for r in cone.rays())
                      and any(r[0] < 0 for r in cone.rays()))
            downgraded.append((upgraded_slice(cone, level), marked))
        source_max, actual_max = marked_maxima(source), marked_maxima(downgraded)
        source_complex = face_closure([p for p, _ in source if p.dim() == 3])
        actual_complex = face_closure([p for p, _ in downgraded if p.dim() == 3])
        counts = (len(actual_max), sum(actual_max.values()))
        check("{} slice: all cells and marks agree ({} maximal, {} marked)".format(name, *counts),
              source_max == actual_max and counts == expected_counts[name]
              and source_complex == actual_complex
              and all(poly_key(c[name]) in source_complex for c in cells.values()
                      if not c[name].is_empty()))

    singular = []
    for cone in general_fan.generating_cones():
        indices = cone_indices(cone, rays)
        if not cone.is_smooth():
            singular.append(indices)
        for dim in range(cone.dim()):
            if any(not face.is_smooth() for face in cone.faces(dim)):
                raise AssertionError("Nonsmooth proper face in cone {}".format(indices))
    check("all proper cone faces are smooth, with exactly one singular maximal cone",
          singular == [(0, 1, 2, 5, 6)])
    print("SINGULAR MAXIMAL CONE:", singular[0])

    anticanonical = Polyhedron(ieqs=[[1] + list(r) for r in rays], base_ring=QQ)
    check("the anticanonical polytope is full-dimensional and bounded",
          anticanonical.dim() == 4 and anticanonical.is_compact())
    cartier = {}
    for cone in general_fan.generating_cones():
        indices = cone_indices(cone, rays)
        character = matrix(QQ, cone.rays()).solve_right(vector(QQ, [-1] * cone.nrays()))
        if not all(x in ZZ for x in character):
            raise AssertionError("Nonintegral Cartier character in cone {}".format(indices))
        for i, ray in enumerate(rays):
            pairing = character.dot_product(vector(QQ, ray))
            if (i in indices and pairing != -1) or (i not in indices and pairing <= -1):
                raise AssertionError("Anticanonical strict convexity fails")
        cartier[indices] = tuple(character)
    check("integral strictly convex Cartier data identify the anticanonical vertices",
          set(cartier.values()) == {tuple(v) for v in anticanonical.vertices()})
    print("ANTICANONICAL CARTIER CHARACTERS:")
    for indices in sorted(cartier):
        print("  {}: {}".format(indices, cartier[indices]))

    points = anticanonical.integral_points()
    grouped = Counter(int(point[1]) for point in points)
    hand_count = {
        -1: sum(binomial(5 - s + 3, 2) for s in range(-2, 2)),
        0: sum(binomial(5 - s, 2) for s in range(-1, 2)),
        1: binomial(2, 2),
    }
    check("the anticanonical lattice count is 130+31+1=162",
          grouped == hand_count == {-1: 130, 0: 31, 1: 1} and len(points) == 162)
    print("LATTICE COUNT BY a:", dict(sorted(grouped.items())))
    print("h^0(-K) =", len(points))
    print("\n{} exact checks passed: the general fibre is toric Fano, with one singular fixed point.".format(len(checks)))


if __name__ == "__main__":
    main()
