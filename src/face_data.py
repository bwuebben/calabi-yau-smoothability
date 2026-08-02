#!/usr/bin/env python3
"""
The named two-dimensional faces printed in paper 3.

Three faces carry the paper's statements and are printed in it, so each is
re-derived here from the polytope's vertex list exactly as that list appears in
the text, and every printed quantity is asserted:

  * the F_1 quadrilateral of Delta_{F_1}          (Theorem 5.1)
  * the dP_6 hexagon of Delta_{F_1}               (Theorem 5.1)
  * the dP_7 pentagon of Delta_19                 (section 4.3)

For each face this module derives

  * its vertices, as a subset of the printed list (with their positions);
  * the two facet normals u_1, u_2, in paper 3's normalization -- the primitive
    inner normals with <u_i, .> = -1 on the face, i.e. the vertices of
    Delta^o = {u : <u,v> >= -1}; these are the negatives of the level-+1
    normals the toolkit emits, and the flip is certified both ways;
  * l(F*) = the lattice length of conv{u_1, u_2};
  * the edge multiset in the induced (saturated) lattice on the affine span;
  * an explicit M in GL_2(Z) carrying that multiset onto the reference polygon
    of the del Pezzo catalogue, identifying C(F) as the anticanonical cone over
    the corresponding surface.

The reference polygons are supplied explicitly and then checked against their
catalogue row (number of vertices, interior points, and rigidity or number of
smoothing components), so no identification rests on an unverified constant.

Importing binds the results to FACES and asserts all of it.  Run directly to
print the derivation:

    python3 src/face_data.py            # ~10 s
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from batyrev_global import (facets, two_faces, face_lattice_polygon,
                            classify_polygon, dual_edge_length)
from toric_census import (ccw_sort, equiv, gl2_bounded, _apply, rigid,
                          smoothing_components, interior_points)

# --------------------------------------------------------------- polytopes
# Delta_{F_1}, exactly as printed in Theorem 5.1.
V_F1 = [(1, 0, 0, 0), (0, 1, 0, 0), (1, -1, 0, 0), (0, 0, 1, 0),
        (0, 0, 0, 1), (0, 0, 1, -1), (0, 0, -1, 1), (0, 0, 0, -1),
        (0, 0, -1, 0), (-1, 1, 0, 0), (0, -1, 0, 0), (-1, 1, -1, 1),
        (0, -1, 0, -1), (0, -1, -1, 0), (-1, 1, -1, 0), (-1, 0, 0, -1),
        (-1, 0, -1, 1), (-1, -1, 0, -1), (-2, 1, -1, 1), (-2, 1, -1, 0),
        (-1, -1, -1, 0), (-2, 0, -1, 0)]

# Delta_19, exactly as printed in section 4.3 (the smallest counterexample).
V_19 = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, -1, 0),
        (0, -1, 0, 0), (1, -1, -1, 0), (0, 0, 0, 1), (-1, 1, 1, -1),
        (0, -1, 0, 1), (0, 0, -1, 1), (-1, 1, 0, -1), (-1, 0, 1, -1),
        (-1, 0, 1, 0), (-1, 1, 0, 0), (0, -1, -1, 0), (0, -1, -1, 1),
        (-1, 0, 0, -1), (-2, 1, 1, -1), (-1, 0, 0, 1)]

# ------------------------------------------------- del Pezzo reference polygons
# Edge multisets.  E_F1 is paper 1's Q_A (Remark 2.6); the hexagon is the one
# asserted in toric_census.main().  Each is checked against its catalogue row
# below, so none of these is taken on faith.
REFERENCE = {
    "F1":  (ccw_sort([(-2, -1), (0, -1), (1, 0), (1, 2)]),
            dict(k=4, i=1, rigid=True)),
    "dP7": (ccw_sort([(-1, 1), (-1, 0), (0, -1), (1, -1), (1, 1)]),
            dict(k=5, i=1, rigid=False, sc=1)),
    "dP6": (ccw_sort([(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]),
            dict(k=6, i=1, rigid=False, sc=2)),
}

for _name, (_E, _row) in REFERENCE.items():
    assert len(_E) == _row["k"], f"{_name}: {len(_E)} edges, catalogue says {_row['k']}"
    assert interior_points(_E) == _row["i"], f"{_name}: interior points"
    assert rigid(_E) is _row["rigid"], f"{_name}: rigidity disagrees with the catalogue"
    if "sc" in _row:
        assert smoothing_components(_E) == _row["sc"], \
            f"{_name}: {smoothing_components(_E)} smoothing components, catalogue says {_row['sc']}"

_G = gl2_bounded(3)


def _faces_of(V):
    facs = facets(V)
    out = []
    for I, fp in two_faces(V, facs):
        u1, u2 = facs[fp[0]][0], facs[fp[1]][0]
        _, evs, lens = face_lattice_polygon(V, I, u1, u2)
        out.append((I, u1, u2, ccw_sort(list(evs)), lens,
                    classify_polygon(evs, lens)))
    return facs, out


def derive(V, want, label):
    """Extract the unique 2-face of V matching `want`, and identify it."""
    facs, faces = _faces_of(V)
    hits = [f for f in faces if want(f[5])]
    assert len(hits) == 1, f"{label}: expected one matching face, got {len(hits)}"
    I, u1, u2, E, lens, c = hits[0]

    verts = [V[i] for i in I]
    entries = [V.index(v) + 1 for v in verts]
    assert all(l == 1 for l in lens), f"{label}: face edges not primitive: {lens}"

    # paper 3's normalization: inner normals taking the value -1 on the face
    U1 = tuple(-x for x in u1)
    U2 = tuple(-x for x in u2)
    for U in (U1, U2):
        assert min(sum(a * b for a, b in zip(U, v)) for v in V) == -1, \
            f"{label}: {U} is not a vertex of Delta^o"
        assert {sum(a * b for a, b in zip(U, v)) for v in verts} == {-1}, \
            f"{label}: {U} is not -1 on the whole face"
    ell = dual_edge_length(U1, U2)
    assert ell == dual_edge_length(u1, u2), f"{label}: negation changed l(F*)"

    # explicit GL_2(Z) identification with the catalogue polygon
    ref_name = None
    for name, (Eref, _) in REFERENCE.items():
        if equiv(E, Eref, _G):
            ref_name = name
            break
    assert ref_name is not None, f"{label}: face matches no catalogue polygon"
    Eref = REFERENCE[ref_name][0]
    # prefer the identity when the multiset is already in standard form, so the
    # paper does not print a gratuitous symmetry
    ident = ((1, 0), (0, 1))
    order = [ident] + [g for g in _G if g != ident]
    M = next((g for g in order if ccw_sort([_apply(g, e) for e in E]) == Eref), None)
    assert M is not None, f"{label}: no explicit matrix realizes the equivalence"
    (a, b), (cc, d) = M
    det = a * d - b * cc
    assert abs(det) == 1, f"{label}: M not unimodular (det {det})"
    assert ccw_sort([_apply(M, e) for e in E]) == Eref

    return dict(label=label, verts=verts, entries=entries, U1=U1, U2=U2,
                diff=tuple(x - y for x, y in zip(U1, U2)), ell=ell,
                E=E, ref=ref_name, Eref=Eref, M=M, det=det, cls=c)


FACES = {
    "F1": derive(V_F1,
                 lambda c: c["status"].startswith("RIGID") and c["k"] == 4,
                 "F_1 quadrilateral of Delta_{F_1}"),
    "dP6": derive(V_F1,
                  lambda c: c["k"] == 6 and c["i"] == 1,
                  "dP_6 hexagon of Delta_{F_1}"),
    "dP7": derive(V_19,
                  lambda c: c["k"] == 5 and c["i"] == 1,
                  "dP_7 pentagon of Delta_19"),
}

for _k, _f in FACES.items():
    assert _f["ref"] == _k, f"{_f['label']}: identified as {_f['ref']}, expected {_k}"
    assert _f["ell"] == 1, f"{_f['label']}: l(F*) = {_f['ell']}, expected 1"


def _report():
    for key in ("F1", "dP6", "dP7"):
        f = FACES[key]
        print(f"== {f['label']} ==")
        for v, e in zip(f["verts"], f["entries"]):
            print(f"     vertex {v}   (entry {e})")
        print(f"     u1 = {f['U1']},  u2 = {f['U2']}")
        print(f"     u1 - u2 = {f['diff']},  l(F*) = {f['ell']}")
        print(f"     edge multiset  {f['E']}")
        print(f"     M = {f['M']},  det = {f['det']}  ->  {f['ref']}: {f['Eref']}")
        print()
    print("face_data: every printed quantity re-derived from the listed vertices")


if __name__ == "__main__":
    _report()
