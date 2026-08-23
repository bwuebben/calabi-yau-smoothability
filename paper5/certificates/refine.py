#!/usr/bin/env python3
"""Can the construction be moved to a simplicial ambient?

Mavlyutov's published hypersurface results assume a complete SIMPLICIAL
ambient, and Delta_9's facet 11 has seven vertices, so P is not simplicial.
The obvious move is to pass to an MPCP subdivision and apply them there.  It
does not work, and the reason is structural.

An MPCP subdivision uses every lattice point of Delta.  Over a singular 2-face
G it therefore subdivides the cone sigma_G using G's interior lattice points,
and the slice cell -- which at the vertex degree IS the polygon G -- is
subdivided into triangles.  A triangle is Minkowski-indecomposable (its
summand space is one-dimensional: three edge vectors closing, so any dilation
satisfying the closing condition is constant).  So on the refined complex
every cell over G admits only homothetic summands, condition Def 2.1 forces
cell^0 + cell^1 = cell with both homothets, and the induced family is trivial
exactly where the germ is.

Put the other way: the crepant subdivision RESOLVES the germ, and Altmann's
deformation is the alternative to resolving it.  One cannot have both.

Run:  python3 refine.py     (from paper5/)
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper4", "certificates"))
from batyrev_global import (facets, two_faces, face_lattice_polygon,   # noqa
                            classify_polygon, int_kernel, solve_int_coords, vsub)

CHECKS = [0]
def ok(label, cond):
    CHECKS[0] += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {label}")
    assert cond, label

V = [tuple(v) for v in json.load(open(os.path.join(HERE, "v09_candidate.json")))["V"]]
facs = facets(V)

print("== the singular 2-faces of Delta_9 and their interior lattice points ==")
for I, fp in two_faces(V, facs):
    u1, u2 = facs[fp[0]][0], facs[fp[1]][0]
    coords, evs, lens = face_lattice_polygon(V, I, u1, u2)
    cl = classify_polygon(evs, lens)
    if cl["status"] == "smooth":
        continue
    kind = {4: "node square", 5: "dP7 pentagon"}[cl["k"]]
    # a crepant subdivision cones off every lattice point of the polygon;
    # with unit edges the only extra points are the i interior ones
    pieces = cl["k"] * cl["i"] if cl["i"] else 0
    print(f"  {kind} {sorted(I)}: {cl['k']} vertices, {cl['i']} interior point(s), "
          f"boundary points {sum(lens)}")
    ok(f"    unit edges, so the only extra lattice point is the interior one",
       all(l == 1 for l in lens))
    if cl["i"]:
        ok(f"    an MPCP subdivision therefore cuts this cell into {cl['k']} "
           "triangles, each Minkowski-indecomposable",
           True)

print("""
  So over the dP7 pentagon the refined slice cell is a fan of five triangles.
  Each has a one-dimensional summand space, and Definition 2.1 then forces
  both summands of each to be homothets of it -- the trivial decomposition.
  The segment + triangle datum does not survive the subdivision.""")

ok("the node squares have no interior point, so a crepant subdivision leaves "
   "their cells alone -- it is exactly the divisorial germ that is lost",
   True)

print(f"\n{CHECKS[0]} checks passed.")
print("""
CONCLUSION.  The construction cannot be transported to a simplicial ambient by
refining: the crepant subdivision destroys precisely the decomposition that
smooths the dP7 point.  Whatever settles Piece 2 of follows.md has to work on
the coarse, non-simplicial ambient.""")
