#!/usr/bin/env sage-python
"""Exact projective partial resolution of the standard Delta_19 hypersurface.

Run: sage -python paper5/certificates/d19_quadric_partial_resolution.py

The pentagon is subdivided into a unimodular triangle and the reflexive
parallelogram defining the anticanonical cone over P1 x P1. All fourteen
nodal cones are preserved. This certifies geometry and relation matrices;
the all-orders deformation argument is in the accompanying paper.
No recursion data or hypersurface coefficient witness is used.
"""

from itertools import combinations, product
import json
from pathlib import Path
import sys

from sage.all import (Cone, Fan, Polyhedron, QQ, ZZ, ceil, gcd, lcm,
                      matrix, singular, vector)
from sage.schemes.toric.ideal import ToricIdeal

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT / "src"))
from face_data import V_19
from batyrev_global import (dual_edge_length, facets, face_lattice_polygon,
                            two_faces)


def ints(v):
    return [int(x) for x in v]


def primitive(v):
    w = vector(ZZ, v * lcm(x.denominator() for x in v))
    return vector(ZZ, w / gcd(list(w)))


def main():
    vs = [vector(QQ, v) for v in V_19]
    fs = facets(V_19)
    face_data = two_faces(V_19, fs)
    node_faces, node_rows = [], []
    for face, pair in face_data:
        if len(face) != 4:
            continue
        ids = sorted(face)
        b = primitive(matrix(QQ, [vs[i] for i in ids]).transpose()
                      .right_kernel().basis()[0])
        assert sorted(b) == [-1, -1, 1, 1]
        assert dual_edge_length(fs[pair[0]][0], fs[pair[1]][0]) == 1
        row = vector(QQ, 19)
        for i, c in zip(ids, b):
            row[i] = c
        node_faces.append(ids)
        node_rows.append(row)
    assert len(node_faces) == 14

    pentagon = set(range(14, 19))
    quad = [15, 16, 17, 18]
    triangle = [14, 15, 16]
    root = vector(QQ, [0] * 15 + [1, -1, 1, -1])
    assert sum((root[i] * vs[i] for i in range(19)), vector(QQ, 4)) == 0
    nodes = matrix(QQ, node_rows)
    mixed = nodes.stack(matrix(QQ, [root]))
    assert nodes.rank() == mixed.rank() == 12
    theta = vector(QQ, [-2, 2, -1, 1, 1, -3, -2, 4,
                        -2, 2, 3, 1, -1, -3, -1])
    assert all(theta) and theta * mixed == 0
    assert mixed.left_kernel().dimension() == 3

    # The original two-dimensional cone block adds one transverse constraint.
    original = nodes
    for b in matrix(QQ, [vs[i] for i in sorted(pentagon)]).transpose().right_kernel().basis():
        row = vector(QQ, 19)
        for i, c in zip(sorted(pentagon), b):
            row[i] = c
        original = original.stack(matrix(QQ, [row]))
    assert original.nrows() == 16 and original.rank() == 13
    assert original.left_kernel().dimension() == 3

    h = vector(QQ, [0, -1, 0, 0, 1, 1, 0, 0, 0, -1,
                    -1, 0, -1, -2, 0, -1, -1, -1, -1])
    assert mixed * h == 0
    ell = matrix(QQ, [vs[i] for i in quad]).solve_right(vector(QQ, [h[i] for i in quad]))
    assert h[14] - ell * vs[14] == 1

    # Every lower cell is enumerated exactly. A common strictly convex
    # support function certifies projectivity of the complete new fan.
    cells, ratio = {}, QQ(0)
    for fi, (n, c, face) in enumerate(fs):
        ids = sorted(face)
        normal = vector(QQ, n) / c
        assert all(normal * vs[i] == 1 for i in ids)
        for four in combinations(ids, 4):
            a = matrix(QQ, [vs[i] for i in four])
            if a.det() == 0:
                continue
            m = a.solve_right(vector(QQ, [h[i] for i in four]))
            if any(m * vs[i] > h[i] for i in ids):
                continue
            cell = tuple(i for i in ids if m * vs[i] == h[i])
            if cell in cells:
                continue
            cells[cell] = dict(parent=fi, linear=m)
            for i in range(19):
                if i not in ids:
                    ratio = max(ratio, (m * vs[i] - h[i]) / (1 - normal * vs[i]))
    shift = ceil(ratio) + 1
    assert shift == 3
    H = h + vector(QQ, [shift] * 19)
    supports = []
    for cell, data in sorted(cells.items()):
        n, c, _ = fs[data["parent"]]
        support = data["linear"] + shift * vector(QQ, n) / c
        assert all(support * vs[i] == H[i] for i in cell)
        gaps = [H[i] - support * vs[i] for i in range(19) if i not in cell]
        assert min(gaps) > 0
        supports.append(dict(rays=list(cell), parent=data["parent"],
                             support=[str(x) for x in support],
                             minimum_strict_gap=str(min(gaps))))
    fan = Fan(cones=list(cells), rays=V_19, check=True)
    assert fan.is_complete() and len(fs) == 27 and len(cells) == 39
    faces3 = {frozenset(c.ambient_ray_indices()) for c in fan(3)}
    assert {f for f in faces3 if f <= pentagon} == {frozenset(triangle), frozenset(quad)}
    assert all(frozenset(f) in faces3 for f in node_faces)
    assert Cone([V_19[i] for i in triangle]).is_smooth()
    for face, _ in face_data:
        if len(face) == 3:
            assert face in faces3 and Cone([V_19[i] for i in face]).is_smooth()

    # Check the parallelogram in the saturated face lattice, not its cover.
    pair = next(pair for face, pair in face_data if set(face) == pentagon)
    q2, edges, lengths = face_lattice_polygon(V_19, quad, fs[pair[0]][0], fs[pair[1]][0])
    assert lengths == [1, 1, 1, 1]
    assert vector(edges[0]) == -vector(edges[2])
    assert vector(edges[1]) == -vector(edges[3])
    assert abs(matrix(ZZ, edges[:2]).det()) == 2
    polygon = Polyhedron(vertices=q2)
    centre = sum((vector(QQ, p) for p in q2), vector(QQ, 2)) / 4
    assert all(x in ZZ for x in centre)
    assert len(polygon.integral_points()) == 5
    centred = [vector(ZZ, vector(QQ, p) - centre) for p in q2]
    assert all(abs(matrix(ZZ, [centred[i], centred[(i + 1) % 4]]).det()) == 1
               for i in range(4))

    # Independent local cotangent-cohomology count. Smooth versal base does
    # NOT mean vanishing T2: the six-dimensional space is an obstruction space.
    hb = Cone([tuple(p) + (1,) for p in q2]).dual().Hilbert_basis()
    ideal = ToricIdeal(matrix(ZZ, list(hb)).transpose())
    assert ideal.ring().ngens() == 9 and len(ideal.gens()) == 20
    singular.lib("sing.lib")
    si = singular(ideal)
    t1 = singular.eval(f"module d19quad_t1=T_1({si.name()});vdim(d19quad_t1);")
    t2 = singular.eval(f"module d19quad_t2=T_2({si.name()});vdim(d19quad_t2);")
    assert int(t1.strip().splitlines()[-1]) == 1
    assert int(t2.strip().splitlines()[-1]) == 6

    # Batyrev h21. Interior points of polar facets / codimension-two faces
    # are recognised by the set of active original vertex inequalities.
    polar = [tuple(-x for x in n) for n, _, _ in fs]
    boxes = [range(min(p[i] for p in polar), max(p[i] for p in polar) + 1) for i in range(4)]
    polar_points = [p for p in product(*boxes) if all(vector(p) * v >= -1 for v in vs)]
    facet_interior, codim2_interior = [], []
    correction = 0
    for p in polar_points:
        active = [v for v in vs if vector(p) * v == -1]
        rank = matrix(QQ, active).rank() if active else 0
        if rank == 1:
            facet_interior.append(p)
        if rank == 2:
            codim2_interior.append(p)
            assert len(active) == 2
            correction += gcd(int(x) for x in active[0] - active[1]) - 1
    assert len(polar_points) == 32 and not facet_interior
    assert len(codim2_interior) == 4 and correction == 0
    h21 = len(polar_points) - 5 - len(facet_interior) + correction
    assert h21 == 27

    result = dict(
        vertices=[list(v) for v in V_19], node_faces=node_faces,
        node_rows=[ints(r) for r in node_rows], quad=quad, triangle=triangle,
        quad_root_row=ints(root), all_nonzero_relation=ints(theta),
        node_rank=12, modified_mixed_rank=12, modified_kernel_dimension=3,
        original_full_rank=13, original_kernel_dimension=3,
        height=ints(h), ear_gap=1, projective_height=ints(H),
        original_maximal_cones=27, refined_maximal_cones=39,
        maximal_cones_and_supports=supports,
        complete=True, all_fourteen_nodes_preserved=True,
        quad_saturated_polygon=[list(p) for p in q2], quad_area=2,
        quad_interior_points=1, quad_hilbert_basis=[ints(v) for v in hb],
        quad_ideal_generators=20, quad_T1_dimension=1, quad_T2_dimension=6,
        polar_lattice_points=32, polar_facet_interior_points=0,
        polar_codimension_two_interior_points=4, polar_correction_term=0,
        crepant_resolution_h21=int(h21),
        scope="Exact toric geometry and linear algebra; deformation proof is separate.")
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS: complete projective fan, 27 -> 39 maximal cones; all 14 nodes preserved")
    print("PASS: pentagon -> smooth triangle + anticanonical P1 x P1 cone")
    print("PASS: mixed rank 12, kernel dimension 3, relation nonzero at all 15 germs")
    print("PASS: local T1 = 1, T2 = 6; h21 of crepant resolution = 27")
    print(output)


if __name__ == "__main__":
    main()
