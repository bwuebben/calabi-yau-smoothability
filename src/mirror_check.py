#!/usr/bin/env python3
"""
Mirror experiment (§6.3): polar-dualize each of the five worked examples and
re-run the face classifier + Hodge numbers on the Batyrev MIRROR.

For a reflexive Delta in Z^4, Batyrev--Borisov mirror symmetry is polar duality
Delta <-> Delta*, and the mirror CY threefold X* is the generic anticanonical
hypersurface in P_{Delta*}.  The vertices of Delta* are exactly the primitive
facet normals of Delta (all at level 1 by reflexivity), which `facets()`
already returns.

For each example we ask the three questions of the plan:
  (1) Is the mirror X* itself smoothable?  -> classify the 2-faces of Delta*
      (R/D/S), and note whether Delta* even has unit edges (else X* has
      non-isolated A_n curve singularities and is 'more singular' outright).
  (2) Where does the non-smoothable 2-face F of Delta go?  Under duality a
      2-face F of Delta is dual to an EDGE F* of Delta*; its own dual-edge
      F* (length l(F*) = #points of X over F) is, on the mirror, a 2-FACE of
      Delta* -- we print that face's type.
  (3) Do the Hodge numbers swap, (h11,h21)(X^) <-> (h21,h11)(X^*), as the
      conifold-transition/mirror heuristic predicts?
"""
import os
import sys
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from batyrev_global import (analyze, facets, two_faces, face_lattice_polygon,
                            classify_polygon, dual_edge_length, vsub, vgcd)
from hodge_numbers import hodge_numbers
from toric_census import ccw_sort, verts_from_edges


def dual_vertices(V):
    """Vertices of Delta* = primitive facet normals of reflexive Delta."""
    facs = facets(V)
    assert all(c == 1 for _, c, _ in facs), "Delta not reflexive"
    return [u for u, _, _ in facs]


def face_report(V):
    """Classify every 2-face of conv(V): return (rows, all_unit_edges)."""
    facs = facets(V)
    rows = []
    for I, fpair in two_faces(V, facs):
        u1, u2 = facs[fpair[0]][0], facs[fpair[1]][0]
        coords, evs, lens = face_lattice_polygon(V, I, u1, u2)
        cl = classify_polygon(evs, lens)
        cl["npoints"] = dual_edge_length(u1, u2)
        rows.append(cl)
    # global edge lengths of conv(V)
    all_unit = True
    for I, fpair in two_faces(V, facs):
        u1, u2 = facs[fpair[0]][0], facs[fpair[1]][0]
        _, _, lens = face_lattice_polygon(V, I, u1, u2)
        if any(l >= 2 for l in lens):
            all_unit = False
    return rows, all_unit, facs


def summarize(rows):
    nonsm = [r for r in rows if "non-smoothable" in r["status"]]
    locsm = [r for r in rows if r["status"].startswith("smoothable")]
    noniso = [r for r in rows if "NON-ISOLATED" in r["status"]]
    return nonsm, locsm, noniso


def verdict(all_unit, nonsm, locsm, noniso):
    if noniso:
        return "NON-ISOLATED (X* has curves of A_n singularities)"
    if nonsm:
        kinds = sorted({("rigid" if r["status"].startswith("RIGID") else "def-only")
                        for r in nonsm})
        return f"NON-SMOOTHABLE (has {'/'.join(kinds)} face)"
    if locsm:
        return "locally smoothable at every point (global-defect question)"
    return "SMOOTH (no singular 2-faces)"


def face_line(r):
    tag = r["status"]
    return (f"k={r['k']} i={r['i']} 2A={r['A2']} npoints(l(F*))={r['npoints']}"
            f"  [{tag}]")


def run(name, V):
    print("=" * 74)
    print(f"### {name}")
    # --- original ---
    rows, unit, facs = face_report(V)
    h = hodge_numbers(V)
    nonsm, locsm, noniso = summarize(rows)
    print(f"  {name:9s}: {len(V)} vtx, {len(facs)} facets, unit-edges={unit}, "
          f"(h11,h21)={h}, chi={2*(h[0]-h[1])}")
    for r in nonsm + noniso:
        print(f"      singular face: {face_line(r)}")
    # --- mirror ---
    Vd = dual_vertices(V)
    rows_d, unit_d, facs_d = face_report(Vd)
    hd = hodge_numbers(Vd)
    nonsm_d, locsm_d, noniso_d = summarize(rows_d)
    print(f"  {name+'*':9s}: {len(Vd)} vtx, {len(facs_d)} facets, "
          f"unit-edges={unit_d}, (h11,h21)={hd}, chi={2*(hd[0]-hd[1])}")
    # tally the mirror's singular 2-faces by type
    from collections import Counter
    cnt = Counter()
    for r in rows_d:
        if r["status"] == "smooth":
            cnt["smooth"] += 1
        elif r["status"].startswith("RIGID"):
            cnt["RIGID"] += 1
        elif r["status"].startswith("def-only"):
            cnt["def-only"] += 1
        elif r["status"].startswith("smoothable"):
            cnt["smoothable(S)"] += 1
        elif "NON-ISOLATED" in r["status"]:
            cnt["NON-ISOLATED"] += 1
    print(f"      mirror 2-faces by type: {dict(cnt)}")
    # long edges of the mirror = transverse A_n curves on X*  (n = length-1)
    rd = analyze(name + "*", Vd, verbose=False)
    long_edges = sorted(l for l in rd["edges"].values() if l >= 2)
    from collections import Counter as _C
    Ln = _C(long_edges)
    maxA = (max(long_edges) - 1) if long_edges else 0
    print(f"      mirror edges >=2 (A_n curves): {dict(Ln)}  "
          f"=> transverse A_{maxA} worst; total {len(long_edges)} bad edges")
    # is the planted face's dual-edge length among the mirror's long edges?
    planted_l = max((r["npoints"] for r in (nonsm + noniso)), default=None)
    if planted_l and planted_l >= 2:
        print(f"      NB: planted l(F*)={planted_l} reappears as a mirror edge of "
              f"length {planted_l}: {'YES' if planted_l in long_edges else 'no'}")
    # Hodge swap check
    swap_ok = (hd == (h[1], h[0]))
    print(f"      Hodge swap (h11,h21)(X^*) == (h21,h11)(X^): {swap_ok}")
    print(f"  => X  : {verdict(unit, nonsm, locsm, noniso)}")
    print(f"  => X* : {verdict(unit_d, nonsm_d, locsm_d, noniso_d)}")
    return dict(name=name, h=h, hd=hd, unit=unit, unit_d=unit_d,
                nonsm=nonsm, nonsm_d=nonsm_d, locsm_d=locsm_d, noniso_d=noniso_d)


def examples():
    F = [(x, y, 1, 0) for (x, y) in [(0, 0), (-2, -1), (-2, -2), (-1, -2)]]
    A = F + [(0, 0, 0, 1), (1, 1, -1, -1)]
    P = [(x, y, 1, 0) for (x, y) in [(0, 0), (-2, -1), (-3, -2), (-2, -3), (-1, -2)]]
    B = P + [(0, 0, 0, 1), (1, 1, -1, -1)]
    w3 = [(0, 0, 0, 1), (1, 1, -1, 0), (-2, -2, 1, -1)]
    C = F + w3
    D = [(1, 0, 0, 0), (0, 1, 0, 0), (0, -1, 0, 0), (0, 0, 1, 0),
         (-3, 1, -1, 0), (0, 1, 1, 0), (0, 0, 0, 1), (-3, 2, 1, -1),
         (-2, 2, 2, -1), (1, 0, 1, 1)]
    sq = verts_from_edges(ccw_sort([(1, 0), (0, 1), (-1, 0), (0, -1)]))
    N = [(x, y, 1, 0) for (x, y) in sq] + [(0, 0, 0, 1), (0, 1, -1, 0),
                                           (-1, -2, 1, -1)]
    return [("Delta_A", A), ("Delta_B", B), ("Delta_C", C),
            ("Delta_D", D), ("Delta_N", N)]


if __name__ == "__main__":
    results = [run(name, V) for name, V in examples()]
    print("=" * 74)
    print("SUMMARY")
    print(f"  {'example':9s} {'X':30s} {'X* (mirror)':34s} Hodge")
    for r in results:
        from collections import Counter
        # short verdicts
        xv = ("non-smoothable" if r["nonsm"] else "smooth-ish")
        if r["noniso_d"]:
            xdv = "NON-ISOLATED (A_n curves)"
        elif r["nonsm_d"]:
            kinds = sorted({("rigid" if rr["status"].startswith("RIGID") else "def-only")
                            for rr in r["nonsm_d"]})
            xdv = "non-smoothable (" + "/".join(kinds) + ")"
        elif r["locsm_d"]:
            xdv = "loc. smoothable (defect Q)"
        else:
            xdv = "smooth"
        print(f"  {r['name']:9s} {xv:30s} {xdv:34s} "
              f"{r['h']}->{r['hd']}")
