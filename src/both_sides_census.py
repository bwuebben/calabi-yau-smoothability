#!/usr/bin/env python3
"""
The definitive both-sides census (paper 3): verify, classify, and assert.

Input: the scan results output/both_sides_v*.json (the complete chain over
all 473,800,776 reflexive 4-polytopes; see both_sides_fast.py) plus the
36-vertex polytope of missing_polytope.py.  Everything here is re-derived
with the exact reference toolkit (batyrev_global.analyze) — the scan is only
used as the list of candidates.

Asserts (= the theorems of paper 3):
  T1  Every 2-face of every both-sides-unit polytope is a triangle, a
      zonotope, or a reflexive polygon (exactly one interior point).
  T2  The non-smoothable faces occurring are exactly the cyclic-quotient
      triangles [3,1] (x6 on one polytope), [3,2] (x10 on one polytope) and
      the F1-cone [4,1] (exactly ONCE, on the unique 22-vertex polytope);
      no def-only face occurs anywhere.
  T3  The F1 polytope's polar (26 vertices) is both-sides with every face
      of type (S) — the unique doubly-isolated mirror pair with a
      non-quotient non-smoothable member — and the MPCP Hodge numbers of
      the pair are mirror-swapped.

Run:  python3 src/both_sides_census.py     (~2 min; everything asserted)
"""
from collections import Counter
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from batyrev_global import analyze, facets                                  # noqa: E402
from toric_census import equiv, gl2_bounded                                 # noqa: E402
from hodge_numbers import hodge_numbers                                     # noqa: E402
from missing_polytope import DELTA36                                        # noqa: E402

F1_EDGES = [(-2, -1), (0, -1), (1, 0), (1, 2)]        # census [4,1] = F1-cone


def is_zonotope(evs):
    c = Counter(tuple(e) for e in evs)
    return all(c[(-a, -b)] == m for (a, b), m in c.items())


def load_polytopes():
    polys, total = [], 0
    for f in sorted(glob.glob(os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..", "output", "both_sides_v*.json"))):
        for key, r in json.load(open(f)).items():
            v = int(key.split("-")[2])
            total += r["n"]
            for h in r["hits"]:
                polys.append((v, [tuple(x) for x in h["verts"]]))
    polys.append((36, [tuple(x) for x in DELTA36]))
    return polys, total + 1                      # +1: Delta_36, verified separately


def main():
    polys, total = load_polytopes()
    assert total == 473_800_776, total
    print(f"scanned: {total:,} reflexive 4-polytopes (complete classification)")
    print(f"both-sides-unit polytopes: {len(polys)}")
    assert len(polys) == 590, len(polys)

    per_v = Counter(v for v, _ in polys)
    print("distribution by vertex count:",
          {v: per_v[v] for v in sorted(per_v)})

    G = gl2_bounded(3)
    occur = Counter()
    nonsmooth_polys = []
    f1_polys = []
    for v, V in polys:
        rep = analyze("p", V, verbose=False)
        assert rep is not None
        assert all(l == 1 for l in rep["edges"].values()), (v, "non-unit edge")
        assert all(f["npoints"] == 1 for f in rep["faces"]), (v, "l(F*)>1")
        ns_here = []
        for f in rep["faces"]:
            k, i, st = f["k"], f["i"], f["status"]
            z = is_zonotope(f["edges2d"])
            occur[(k, i, z, st)] += 1
            # --- T1 ---
            assert k == 3 or z or i == 1, ("T1 violated", v, k, i, st)
            # --- T2 bookkeeping ---
            assert "def-only" not in st, ("def-only on both-sides!", v)
            if "non-smooth" in st:
                ns_here.append((k, i))
                if (k, i) == (4, 1):
                    assert equiv(tuple(sorted(f["edges2d"])),
                                 tuple(sorted(F1_EDGES)), G), "not the F1 class"
                    f1_polys.append((v, V))
        if ns_here:
            nonsmooth_polys.append((v, sorted(ns_here)))

    print("\nface-class occurrences (k, i, zonotope?, status):")
    for (k, i, z, st), n in sorted(occur.items()):
        print(f"  k={k} i={i} {'zonotope' if z else 'asym.':8s} {st:<26s} x{n}")

    # --- T2 ---
    assert sorted(nonsmooth_polys) == sorted([
        (5,  [(3, 2)] * 10),
        (9,  [(3, 1)] * 6),
        (22, [(3, 2)] * 0 + [(4, 1)]),
    ]), nonsmooth_polys
    assert len(f1_polys) == 1
    print("\nT1 PROVED: every face is a triangle, a zonotope, or reflexive (i=1).")
    print("T2 PROVED: non-smoothable both-sides germs = quotient triangles")
    print("           ([3,1] x6 on one 9-vertex polytope; [3,2] x10 on one")
    print("           5-vertex polytope) and the F1-cone (x1, 22 vertices);")
    print("           def-only faces: ZERO database-wide.")

    # --- T3: the unique F1 mirror pair ---
    v, V = f1_polys[0]
    assert v == 22
    polar = [u for (u, c, _) in facets(V)]
    repp = analyze("polar", polar, verbose=False)
    assert len(polar) == 26 and repp is not None
    assert all(l == 1 for l in repp["edges"].values())
    assert all(f["npoints"] == 1 for f in repp["faces"])
    assert all("non-smooth" not in f["status"] for f in repp["faces"]), \
        "polar carries a non-smoothable face?!"
    cX = Counter((f["k"], f["i"], f["status"]) for f in
                 analyze("X", V, verbose=False)["faces"])
    cM = Counter((f["k"], f["i"], f["status"]) for f in repp["faces"])
    print("\nT3: the unique F1 mirror pair (22 <-> 26 vertices):")
    print("  X  faces:", dict(cX))
    print("  X* faces:", dict(cM))
    h11, h21 = hodge_numbers(V, "F1 side")
    h11m, h21m = hodge_numbers(polar, "mirror side")
    assert (h11, h21) == (h21m, h11m), "Hodge numbers not mirror-swapped"
    print(f"  MPCP Hodge numbers: X-hat ({h11},{h21})  <->  mirror-hat "
          f"({h11m},{h21m})  [mirror-swapped, as they must be]")
    print("\nall both-sides census asserts pass")


if __name__ == "__main__":
    main()
