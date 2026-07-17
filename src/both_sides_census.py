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

# ---- literal data as printed in paper 3 (asserted against the scan) ----
D19_PAPER = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, -1, 0),
             (0, -1, 0, 0), (1, -1, -1, 0), (0, 0, 0, 1), (-1, 1, 1, -1),
             (0, -1, 0, 1), (0, 0, -1, 1), (-1, 1, 0, -1), (-1, 0, 1, -1),
             (-1, 0, 1, 0), (-1, 1, 0, 0), (0, -1, -1, 0), (0, -1, -1, 1),
             (-1, 0, 0, -1), (-2, 1, 1, -1), (-1, 0, 0, 1)]
DF1_PAPER = [(1, 0, 0, 0), (0, 1, 0, 0), (1, -1, 0, 0), (0, 0, 1, 0),
             (0, 0, 0, 1), (0, 0, 1, -1), (0, 0, -1, 1), (0, 0, 0, -1),
             (0, 0, -1, 0), (-1, 1, 0, 0), (0, -1, 0, 0), (-1, 1, -1, 1),
             (0, -1, 0, -1), (0, -1, -1, 0), (-1, 1, -1, 0), (-1, 0, 0, -1),
             (-1, 0, -1, 1), (-1, -1, 0, -1), (-2, 1, -1, 1), (-2, 1, -1, 0),
             (-1, -1, -1, 0), (-2, 0, -1, 0)]
S5_PAPER = [(1, 0, 0, 0), (0, 1, 0, 0), (2, 4, 5, 0), (3, 3, 0, 5),
            (-6, -8, -5, -5)]

RIG = "RIGID (non-smoothable)"
S1, S2, S3 = ("smoothable (1 comps)", "smoothable (2 comps)",
              "smoothable (3 comps)")
TAB1_PAPER = {5: 2, 6: 1, 8: 1, 9: 1, 10: 1, 12: 1, 15: 1, 16: 4, 17: 4,
              18: 7, 19: 12, 20: 24, 21: 37, 22: 69, 23: 84, 24: 98, 25: 85,
              26: 64, 27: 42, 28: 24, 29: 13, 30: 9, 31: 4, 33: 1, 36: 1}
TAB2_PAPER = {(3, 0, False, "smooth"): 25873, (3, 1, False, RIG): 6,
              (3, 2, False, RIG): 10, (4, 0, True, S1): 15652,
              (4, 1, False, RIG): 1, (4, 1, True, S1): 99,
              (4, 2, True, S1): 9, (5, 1, False, S1): 87,
              (6, 1, True, S2): 395, (6, 2, True, S1): 2,
              (6, 7, True, S1): 3, (8, 4, True, S3): 1}


def is_zonotope(evs):
    c = Counter(tuple(e) for e in evs)
    return all(c[(-a, -b)] == m for (a, b), m in c.items())


def load_polytopes():
    polys, total, scanned_by_v = [], 0, Counter()
    for f in sorted(glob.glob(os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..", "output", "both_sides_v*.json"))):
        for key, r in json.load(open(f)).items():
            v = int(key.split("-")[2])
            total += r["n"]
            scanned_by_v[v] += r["n"]
            for h in r["hits"]:
                polys.append((v, [tuple(x) for x in h["verts"]]))
    polys.append((36, [tuple(x) for x in DELTA36]))
    scanned_by_v[36] += 1
    return polys, total + 1, scanned_by_v       # +1: Delta_36, verified separately


def main():
    polys, total, scanned_by_v = load_polytopes()
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
    records = []                       # (v, V, per-poly face Counter, rep)
    for v, V in polys:
        rep = analyze("p", V, verbose=False)
        assert rep is not None
        assert all(l == 1 for l in rep["edges"].values()), (v, "non-unit edge")
        assert all(f["npoints"] == 1 for f in rep["faces"]), (v, "l(F*)>1")
        ns_here, cnt = [], Counter()
        for f in rep["faces"]:
            k, i, st = f["k"], f["i"], f["status"]
            z = is_zonotope(f["edges2d"])
            occur[(k, i, z, st)] += 1
            cnt[(k, i, z, st)] += 1
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
        records.append((v, V, cnt, rep))
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

    # ---- every remaining number of paper 3, asserted (adversarial round) ----
    assert dict(per_v) == TAB1_PAPER, per_v                       # Table 1
    assert dict(occur) == TAB2_PAPER, occur                       # Table 2
    assert sum(n for w, n in scanned_by_v.items() if w <= 18) == 395_406_329
    assert sum(n for w, n in scanned_by_v.items() if w <= 8) == 1_037_834

    # asymmetric (neither triangle nor centrally symmetric) faces: 88 on 55
    asym = [(v, V, cnt) for v, V, cnt, _ in records
            if any(k >= 4 and not z for (k, i, z, st) in cnt)]
    n_asym_faces = sum(m for v, V, cnt in asym
                       for (k, i, z, st), m in cnt.items() if k >= 4 and not z)
    assert n_asym_faces == 88 and len(asym) == 55, (n_asym_faces, len(asym))
    assert min(v for v, _, _ in asym) == 19       # conjecture TRUE for v <= 18
    dp7 = [(v, cnt) for v, V, cnt in asym if cnt[(5, 1, False, S1)]]
    assert sum(c[(5, 1, False, S1)] for _, c in dp7) == 87
    assert (min(v for v, _ in dp7), max(v for v, _ in dp7)) == (19, 31)

    # the smallest counterexample Delta_19: unique at 19 vertices,
    # 14 nodes + one dP7-cone point, matching the paper's coordinates
    v19 = [(V, cnt) for v, V, cnt in asym if v == 19]
    assert len(v19) == 1 and set(v19[0][0]) == set(D19_PAPER)
    c19 = v19[0][1]
    assert c19[(4, 0, True, S1)] == 14 and c19[(5, 1, False, S1)] == 1
    assert all((k, i) == (3, 0) for (k, i, z, st) in c19
               if (k, i, z, st) not in ((4, 0, True, S1), (5, 1, False, S1)))

    # the low-vertex stratum (paper 3, Example: two polar pairs + one
    # self-polar polytope among the 3,905,789 polytopes with <= 9 vertices)
    assert sum(n for w, n in scanned_by_v.items() if w <= 9) == 3_905_789
    low = [r for r in records if r[0] <= 9]
    assert Counter(r[0] for r in low) == Counter({5: 2, 6: 1, 8: 1, 9: 1})
    s5 = [r for r in low if set(r[1]) == set(S5_PAPER)]
    assert len(s5) == 1 and dict(s5[0][2]) == {(3, 2, False, RIG): 10}
    pol5 = [u for (u, c, _) in facets(s5[0][1])]
    rp5 = analyze("polar(S5)", pol5, verbose=False)
    assert len(pol5) == 5 and all(f["status"] == "smooth" for f in rp5["faces"])
    other5 = [r for r in low if r[0] == 5 and set(r[1]) != set(S5_PAPER)]
    assert len(other5) == 1 and all(st == "smooth"
                                    for (k, i, z, st) in other5[0][2])
    r6 = next(r for r in low if r[0] == 6)
    assert all(st == "smooth" for (k, i, z, st) in r6[2])
    assert all(len(idx) == 4 for (u, c, idx) in facets(r6[1]))    # simplicial
    pol6 = [u for (u, c, _) in facets(r6[1])]
    rp6 = analyze("polar(v6)", pol6, verbose=False)
    cp6 = Counter((f["k"], f["i"], is_zonotope(f["edges2d"]), f["status"])
                  for f in rp6["faces"])
    c9 = next(r for r in low if r[0] == 9)[2]
    assert len(pol6) == 9 and cp6 == c9
    assert c9[(3, 1, False, RIG)] == 6 and c9[(4, 2, True, S1)] == 9
    # the fifth polytope: 8 vertices, 16 facets, 32 smooth triangles; its
    # polar is a SIXTEEN-vertex both-sides polytope (24 diamond faces) --
    # the <=9 window is NOT closed under polarity (adversarial-round catch)
    r8 = next(r for r in low if r[0] == 8)
    pol8 = [u for (u, c, _) in facets(r8[1])]
    rp8 = analyze("polar(v8)", pol8, verbose=False)
    cp8 = Counter((f["k"], f["i"], is_zonotope(f["edges2d"]), f["status"])
                  for f in rp8["faces"])
    assert len(pol8) == 16
    assert dict(r8[2]) == {(3, 0, False, "smooth"): 32}
    assert dict(cp8) == {(4, 1, True, S1): 24}
    assert any(r[2] == cp8 for r in records if r[0] == 16)   # partner exists

    # the top of the classification: H x H and its polar free sum
    r36 = next(r for r in records if r[0] == 36)
    assert dict(r36[2]) == {(4, 0, True, S1): 36, (6, 1, True, S2): 12}
    pol36 = [u for (u, c, _) in facets(r36[1])]
    rp36 = analyze("polar(HxH)", pol36, verbose=False)
    assert len(pol36) == 12 and len(rp36["faces"]) == 72
    assert all(f["k"] == 3 and f["status"] == "smooth" for f in rp36["faces"])

    # smooth x smooth pairs (both hypersurfaces smooth, no resolution needed):
    # exactly ONE in the whole classification -- the 24-cell (Braun's (20,20)
    # self-mirror threefold), paper 3's Example "the unique smooth mirror pair"
    allsm = [r for r in records
             if all(st == "smooth" for (k, i, z, st) in r[2])]
    both_smooth = []
    for v, V, cnt, rep in allsm:
        fs = facets(V)
        rr = analyze("polar", [u for (u, c, _) in fs], verbose=False)
        if all(f["status"] == "smooth" for f in rr["faces"]):
            both_smooth.append((v, V, fs, rr))
    assert len(allsm) == 8 and len(both_smooth) == 1
    v24, V24, fs24, rp24 = both_smooth[0]
    assert v24 == 24 and len(fs24) == 24                 # self-dual f-vector
    assert all(len(idx) == 6 for (u, c, idx) in fs24)    # octahedral facets
    assert len(rp24["faces"]) == 96
    h24 = hodge_numbers(V24, "24-cell")
    h24m = hodge_numbers([u for (u, c, _) in fs24], "24-cell polar")
    assert h24 == (20, 20) == h24m                       # the self-mirror CY
    print(f"\nall-faces-smooth polytopes: {len(allsm)}; both-sides smooth: "
          f"exactly 1 -- the 24-cell, Hodge (20,20) x (20,20)")
    print("PAPER-ASSERTS: every number of paper 3 re-derived and asserted.")

    # --- T3: the unique F1 mirror pair ---
    v, V = f1_polys[0]
    assert v == 22 and set(V) == set(DF1_PAPER)   # matches the paper verbatim
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
    assert dict(cX) == {(3, 0, "smooth"): 50, (4, 0, S1): 20,
                        (4, 1, RIG): 1, (6, 1, S2): 1}
    assert dict(cM) == {(3, 0, "smooth"): 38, (4, 0, S1): 26,
                        (5, 1, S1): 2, (6, 1, S2): 2}
    assert sum(cX.values()) == 72 and sum(cM.values()) == 68  # 2-face totals
    h11, h21 = hodge_numbers(V, "F1 side")
    h11m, h21m = hodge_numbers(polar, "mirror side")
    assert (h11, h21) == (h21m, h11m), "Hodge numbers not mirror-swapped"
    assert (h11, h21) == (20, 26), (h11, h21)
    print(f"  MPCP Hodge numbers: X-hat ({h11},{h21})  <->  mirror-hat "
          f"({h11m},{h21m})  [mirror-swapped, as they must be]")
    print("\nall both-sides census asserts pass")


if __name__ == "__main__":
    main()
