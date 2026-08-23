import os
#!/usr/bin/env python3
"""The one 9-vertex census candidate: inside the Paper 4 framework AND with a
decomposable pentagon-bearing facet.

Found by dp7_facet_sweep.py on polytopes-4d-09-vertices.parquet: of 11 596
pentagon 2-faces, 63 are inside the framework and 2 350 sit in a facet that
decomposes -- exactly one is both.
"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper4", "certificates"))
from batyrev_global import (facets, two_faces, face_lattice_polygon,   # noqa
                            classify_polygon, dual_edge_length, dot)
from examples import run_example, analyze_faces                        # noqa

V = [tuple(v) for v in json.load(open(os.path.join(HERE, "v09_candidate.json")))["V"]]
facs = facets(V)
print(f"V = {V}")
print(f"{len(V)} vertices, {len(facs)} facets\n")

print("== every singular 2-face ==")
sing = []
for I, fp in two_faces(V, facs):
    _, evs, lens = face_lattice_polygon(V, I, facs[fp[0]][0], facs[fp[1]][0])
    cl = classify_polygon(evs, lens)
    if cl["status"] == "smooth":
        continue
    dl = dual_edge_length(facs[fp[0]][0], facs[fp[1]][0])
    sing.append((frozenset(I), cl, dl))
    kind = {4: "node", 5: "dP7 ", 6: "dP6 "}.get(cl["k"], f"{cl['k']}-gon")
    print(f"  {kind} {sorted(I)}: i={cl['i']} edges {lens} dual {dl}  {cl['status']}")
nodes = sum(dl for _, cl, dl in sing if cl["k"] == 4)
dp7 = sum(dl for _, cl, dl in sing if cl["k"] == 5)
dp6 = sum(dl for _, cl, dl in sing if cl["k"] == 6)
print(f"\n  hypersurface: {nodes} nodes, {dp7} dP7 points, {dp6} dP6 points")

print("\n== reachability at each vertex degree (one_facet.py's lemma) ==")
best = (-1, None)
for fi, (u, c, idx) in enumerate(facs):
    R = tuple(-x for x in u)
    vals = {i: dot(R, V[i]) for i in range(len(V))}
    reach = [S for S, _, _ in sing if all(vals[i] <= 0 for i in S)]
    bounded = [S for S, _, _ in sing if all(vals[i] < 0 for i in S)]
    if len(reach) > best[0]:
        best = (len(reach), fi)
    print(f"  facet {fi:>2} {sorted(idx)}: {len(reach)} of {len(sing)} reachable, "
          f"{len(bounded)} with a bounded cell")
print(f"  best: facet {best[1]} reaches {best[0]} of {len(sing)}")

print("\n== the Paper 4 necessity test ==")
run_example("the 9-vertex candidate", V)

print("\n== the two facets containing the pentagon ==")
PENT = [3, 4, 5, 6, 8]
import itertools
from fractions import Fraction as Fr
from examples import rank as _rank
pf = [f for I, fp in two_faces(V, facs) if sorted(I) == PENT for f in fp]
for fi in pf:
    tfl = [sorted(I) for I, fp in two_faces(V, facs) if fi in fp]
    eds = sorted({tuple(sorted(set(A) & set(B))) for A, B in itertools.combinations(tfl, 2)
                  if len(set(A) & set(B)) == 2})
    ei = {e: j for j, e in enumerate(eds)}
    cyc = {}
    for A in tfl:
        adj = {}
        for e in eds:
            if set(e) <= set(A):
                adj.setdefault(e[0], []).append(e[1]); adj.setdefault(e[1], []).append(e[0])
        c = [A[0]]; pv, cu = None, A[0]
        while len(c) < len(A):
            nx = [x for x in adj[cu] if x != pv][0]; c.append(nx); pv, cu = cu, nx
        cyc[tuple(A)] = c
    rows = []
    for A in tfl:
        cc = cyc[tuple(A)]
        for co in range(4):
            row = [Fr(0)] * len(eds)
            for j in range(len(cc)):
                a, b = cc[j], cc[(j + 1) % len(cc)]
                row[ei[tuple(sorted((a, b)))]] += V[b][co] - V[a][co]
            rows.append(row)
    dim = len(eds) - _rank(rows)
    types = sorted(len(A) for A in tfl)
    print(f"  facet {fi} {sorted(facs[fi][2])}: 2-faces {types}, {len(eds)} edges, "
          f"summand dim {dim}")
    if dim <= 1:
        print("     rigid")
        continue
    sols = [t for t in itertools.product((0, 1), repeat=len(eds))
            if not (all(x == 0 for x in t) or all(x == 1 for x in t))
            and all(sum(Fr(t[j]) * r[j] for j in range(len(eds))) == 0 for r in rows)]
    print(f"     {len(sols)} nontrivial lattice decompositions")
    for t in sols:
        desc = []
        for A in tfl:
            cc = cyc[tuple(A)]
            vals = [t[ei[tuple(sorted((cc[j], cc[(j + 1) % len(cc)])))]]
                    for j in range(len(cc))]
            if len(set(vals)) == 2:
                k = len(cc)
                kind = {4: "node", 5: "dP7"}.get(k, f"{k}-gon")
                if kind != f"{k}-gon":
                    desc.append(f"{kind}{tuple(A)} {sum(vals)}|{k - sum(vals)}")
        print(f"       t={t}: {'; '.join(desc) if desc else 'no singular face split'}")
