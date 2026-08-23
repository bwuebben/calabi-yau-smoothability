import os
"""Do any of the census candidates let the ambient reach the dP7 point?

dp7_facet_sweep.py finds reflexive 4-polytopes whose pentagon-bearing facet
is Minkowski-DECOMPOSABLE, unlike Delta_19's.  That is necessary but not
sufficient: by Ilten-Vollmert Definition 4.1(ii) the pentagon's summands are
faces of the facet's summands, so what matters is the decomposition INDUCED
on the pentagon.  Two questions per candidate:

  (a) does the facet's summand cone project onto more than the homothets of
      the pentagon (i.e. is the induced pentagon decomposition ever
      nontrivial)?
  (b) is the segment+triangle dilation -- the one Altmann's theory attaches
      to the smoothing of the dP7 cone (Tohoku (7.1.4)) -- attainable, with
      all dilations nonnegative on the whole facet?

(b) is an exact rational feasibility question, solved as a polyhedron.

Run:  sage candidates.sage     (from paper5/)
"""
import json, os, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper4", "certificates"))
from batyrev_global import (facets, two_faces, face_lattice_polygon,
                            classify_polygon, dual_edge_length)

cands = json.load(open(os.path.join(HERE, "dp7_sweep_v08.json")))
print(f"{len(cands)} candidates from the 8-vertex sweep\n")

def facet_complex(twofaces):
    edges = set()
    for A, B in itertools.combinations(twofaces, 2):
        c = tuple(sorted(set(A) & set(B)))
        if len(c) == 2:
            edges.add(c)
    edges = sorted(edges)
    cycles = []
    for A in twofaces:
        adj = {}
        for e in edges:
            if set(e) <= set(A):
                adj.setdefault(e[0], []).append(e[1])
                adj.setdefault(e[1], []).append(e[0])
        if sorted(adj) != sorted(A) or any(len(adj[v]) != 2 for v in A):
            return None, None
        cyc = [A[0]]; prev, cur = None, A[0]
        while len(cyc) < len(A):
            nxt = [x for x in adj[cur] if x != prev][0]
            cyc.append(nxt); prev, cur = cur, nxt
        cycles.append(cyc)
    return edges, cycles

summary = dict(nontrivial=0, seg_tri=0, dead=0, skipped=0)
winners = []
reasons = []
for k, c in enumerate(cands):
    Vp = [tuple(int(x) for x in v) for v in c["V"]]   # plain ints for batyrev_global
    V = [vector(QQ, v) for v in Vp]                   # Sage vectors for the algebra
    facs = facets(Vp)
    tf = list(two_faces(Vp, facs))
    by_facet = {}
    for I, fp in tf:
        for f in fp:
            by_facet.setdefault(f, []).append(sorted(I))
    PENT = c["pent"]
    for fi, dim in zip([f for I, fp in tf if sorted(I) == PENT for f in fp], c["dims"]):
        if dim is None or dim < 2:
            continue
        edges, cycles = facet_complex(by_facet[fi])
        if edges is None:
            summary["skipped"] += 1; continue
        eidx = {e: j for j, e in enumerate(edges)}
        # V(Q): edge dilations closing on every 2-face
        eqns = []
        for cyc in cycles:
            for co in range(4):
                row = [QQ(0)] * len(edges)
                for j in range(len(cyc)):
                    a, b = cyc[j], cyc[(j + 1) % len(cyc)]
                    row[eidx[tuple(sorted((a, b)))]] += (V[b] - V[a])[co]
                eqns.append([0] + row)
        # the pentagon's own edges, in its cyclic order
        pcyc = [cyc for cyc in cycles if sorted(cyc) == PENT][0]
        pedges = [tuple(sorted((pcyc[j], pcyc[(j + 1) % 5]))) for j in range(5)]
        # (a) does the summand SPACE restrict to more than constants on them?
        K = matrix(QQ, [r[1:] for r in eqns]).right_kernel().basis_matrix()
        proj = K.matrix_from_columns([eidx[e] for e in pedges])
        if proj.rank() <= 1:
            summary["dead"] += 1
            continue
        summary["nontrivial"] += 1
        # (b) is the segment+triangle pattern attainable with t >= 0 ?
        # the segment is the pair of antiparallel pentagon edges
        E = [V[pcyc[(j + 1) % 5]] - V[pcyc[j]] for j in range(5)]
        seg = [(i, j) for i in range(5) for j in range(i + 1, 5)
               if matrix(QQ, [list(E[i]), list(E[j])]).rank() == 1]
        hit = False
        for (i, j) in seg:
            pat = [1 if x in (i, j) else 0 for x in range(5)]
            extra = [[-pat[x]] + [1 if y == eidx[pedges[x]] else 0
                                  for y in range(len(edges))] for x in range(5)]
            P = Polyhedron(eqns=eqns + extra,
                           ieqs=[[0] + [1 if y == z else 0 for y in range(len(edges))]
                                 for z in range(len(edges))], base_ring=QQ)
            if not P.is_empty():
                hit = True
        if hit:
            summary["seg_tri"] += 1
            inv, bad = [], []
            for I, fp in tf:
                _, evs, lens = face_lattice_polygon(Vp, I,
                                                    facs[fp[0]][0], facs[fp[1]][0])
                cl = classify_polygon(evs, lens)
                dl = dual_edge_length(facs[fp[0]][0], facs[fp[1]][0])
                if cl["status"] == "smooth":
                    continue                      # a UNIMODULAR triangle only
                # everything else is a singular face and must be inside the
                # Paper 4 framework: unit edges, isolated germ, one point
                if any(l != 1 for l in lens) or dl != 1 or \
                   not ((cl["k"] == 4 and cl["i"] == 0) or
                        (cl["k"] in (5, 6) and cl["i"] == 1)):
                    bad.append((cl["k"], cl["i"], cl["status"][:14], dl, max(lens)))
                inv.append((cl["k"], cl["i"], cl["status"][:14], dl))
            clean = (not bad) and \
                    all("RIGID" not in st and "def-only" not in st
                        for _, _, st, _ in inv)
            reasons.append(tuple(sorted({("A_n edges" if ml > 1 else
                                          "dual length %d" % dl if dl > 1 else
                                          "%s" % st) for _, _, st, dl, ml in bad})))
            winners.append((k, fi, c["dual_len"], sorted(inv), clean, Vp))

print("of the pentagon-bearing facets that decompose:")
print(f"   {summary['dead']:>4} induce only homothets on the pentagon (still dead)")
print(f"   {summary['nontrivial']:>4} induce a nontrivial decomposition of the pentagon")
print(f"   {summary['seg_tri']:>4} of those attain the segment+triangle dilation")
print(f"   {summary['skipped']:>4} skipped (degenerate incidence)")
from collections import Counter
print("\nwhy the survivors fall outside the Paper 4 framework:")
for r, n in Counter(reasons).most_common():
    print(f"   {n:>4}  {'; '.join(r) if r else 'INSIDE THE FRAMEWORK'}")
clean = [w for w in winners if w[4]]
print(f"\nof the {len(winners)} that attain it, {len(clean)} have an "
      "inventory inside the Paper 4 framework\n(only nodes and dP6/dP7 points, all isolated), "
      "with every dual edge of length one.")
json.dump([{"V": w[5], "facet": w[1], "inv": [list(x) for x in w[3]]} for w in clean],
          open(os.path.join(HERE, "dp7_winners.json"), "w"), indent=1)
print(f"\n{len(clean)} written to dp7_winners.json")
for w in clean[:10]:
    print("   ", w[0], w[1], sorted(w[3]))
