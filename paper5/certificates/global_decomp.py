#!/usr/bin/env python3
"""Does Delta_19 itself decompose?  (The Mavlyutov route.)

route_b.md records Mavlyutov arXiv:0902.0967 Section 7 as the right
multi-face object: ONE global lattice Minkowski decomposition of the
reflexive polytope, inducing decompositions of every face at once.  Before
anything is proved about that construction, its combinatorial input can be
tested: does Delta_19 admit a nontrivial Minkowski decomposition at all,
and does any decomposition induce the segment+triangle splitting on the
pentagon that the dP7 point needs?

Method.  A weak Minkowski summand of a polytope assigns a dilation t_e >= 0
to each edge, and every 2-face must close: the cyclically traversed edge
vectors, weighted by t, sum to zero.  Any actual summand gives such a
datum, so a 1-dimensional solution space (only the homothets t = const)
PROVES indecomposability -- which is the direction used here.

Run:  python3 global_decomp.py     (from paper5/)
"""
import itertools
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper4", "certificates"))

from batyrev_global import facets, two_faces                       # noqa
from examples import V_19, V_20, rank                              # noqa

CHECKS = [0]
def ok(label, cond):
    CHECKS[0] += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {label}")
    assert cond, label


def skeleton(V):
    """2-faces (cyclically ordered) and edges of a 4-polytope."""
    facs = facets(V)
    tfs = [sorted(I) for I, _ in two_faces(V, facs)]
    edges = set()
    for A, B in itertools.combinations(tfs, 2):
        common = tuple(sorted(set(A) & set(B)))
        if len(common) == 2:
            edges.add(common)
    edges = sorted(edges)
    cycles = []
    for A in tfs:
        adj = {}
        for e in edges:
            if set(e) <= set(A):
                adj.setdefault(e[0], []).append(e[1])
                adj.setdefault(e[1], []).append(e[0])
        assert sorted(adj) == A and all(len(adj[v]) == 2 for v in A), (A, adj)
        cyc = [A[0]]; prev, cur = None, A[0]
        while len(cyc) < len(A):
            nxt = [x for x in adj[cur] if x != prev][0]
            cyc.append(nxt); prev, cur = cur, nxt
        cycles.append(cyc)
    return facs, tfs, edges, cycles


def summand_dim(V, tfs, edges, cycles, restrict=None):
    """dim of the space of edge dilations closing on every 2-face.
    restrict: optional dict edge -> fixed value (appended as equations)."""
    eidx = {e: k for k, e in enumerate(edges)}
    rows = []
    for cyc in cycles:
        for c in range(4):
            row = [Fr(0)] * len(edges)
            for j in range(len(cyc)):
                a, b = cyc[j], cyc[(j + 1) % len(cyc)]
                row[eidx[tuple(sorted((a, b)))]] += V[b][c] - V[a][c]
            rows.append(row)
    return len(edges) - rank(rows)


# ------------------------------------------------- validation in dimension 4
print("== the routine, on 4-polytopes with known summand spaces ==")
cube4 = [tuple(b) for b in itertools.product((0, 1), repeat=4)]
simp4 = [(0,0,0,0), (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1)]
tri = [(0,0), (1,0), (0,1)]
tritri = [a + b for a in tri for b in tri]
tetseg = [(a, b, c, d) for (a,b,c) in [(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
          for d in (0, 1)]
for name, P, want in [("4-cube = four segments", cube4, 4),
                      ("4-simplex (indecomposable)", simp4, 1),
                      ("triangle x triangle", tritri, 2),
                      ("tetrahedron x segment", tetseg, 2)]:
    facs, tfs, edges, cycles = skeleton(P)
    d = summand_dim(P, tfs, edges, cycles)
    ok(f"{name}: summand dim {d} = {want}", d == want)

def polar(V):
    return [tuple(-x for x in u) for u, c, _ in facets(V)]

print()
for name, VV in [("Delta_19", V_19), ("polar of Delta_19", polar([tuple(v) for v in V_19])),
                 ("Delta_20", V_20), ("polar of Delta_20", polar([tuple(v) for v in V_20]))]:
    V = [tuple(v) for v in VV]
    facs, tfs, edges, cycles = skeleton(V)
    f0, f1, f2, f3 = len(V), len(edges), len(tfs), len(facs)
    print(f"\n== {name} ==")
    ok(f"f-vector ({f0}, {f1}, {f2}, {f3}) satisfies Euler f0-f1+f2-f3 = 0",
       f0 - f1 + f2 - f3 == 0)
    d = summand_dim(V, tfs, edges, cycles)
    ok(f"summand space of the whole 4-polytope has dimension {d}"
       f"{' -- INDECOMPOSABLE' if d == 1 else ' -- decomposable, follow up'}",
       d >= 1)

print(f"\n{CHECKS[0]} checks passed.")
print("""
A 1-dimensional summand space means the polytope is Minkowski-indecomposable:
the only weak summands are its own homothets.  A global lattice decomposition
of the kind Mavlyutov's Section 7 takes as input then does not exist, so that
construction has no nontrivial input on this polytope -- independently of
whether its unproved statements can be proved.""")
