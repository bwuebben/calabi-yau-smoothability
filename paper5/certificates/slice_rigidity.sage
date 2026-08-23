#!/usr/bin/env sage
"""Independent Sage cross-check of slice_rigidity.py.

Uses Sage's convex hull (not the combinatorics of Delta_19) to build the
slice cells, and settles the all-degrees statement a second way: over the
function field Q(s), where s = 1/(1-t) for facet 21 and s = 1/t for facet
24 parametrises the rescaling of the off-pentagon vertices.

Run:  sage slice_rigidity.sage     (from paper5/)
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper4", "certificates"))
from batyrev_global import facets
from examples import V_19

V = [vector(ZZ, v) for v in V_19]
facs = facets(V_19)
PENT = [14, 15, 16, 17, 18]
CH = [0]
def ok(label, cond):
    CH[0] += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {label}")
    assert cond, label

def cell_combinatorics(Q):
    verts = [vector(QQ, v) for v in Q.vertices_list()]
    vidx = {tuple(v): k for k, v in enumerate(verts)}
    edges = sorted(tuple(sorted(vidx[tuple(vector(QQ, x))] for x in e.vertices()))
                   for e in Q.faces(1))
    eidx = {e: k for k, e in enumerate(edges)}
    faces = []
    for f in Q.faces(2):
        fv = set(vidx[tuple(vector(QQ, x))] for x in f.vertices())
        adj = {}
        for (a, b) in edges:
            if a in fv and b in fv:
                adj.setdefault(a, []).append(b); adj.setdefault(b, []).append(a)
        assert all(len(adj[x]) == 2 for x in fv)
        start = min(fv); cyc = [start]; prev = None; cur = start
        while len(cyc) < len(fv):
            nxt = [x for x in adj[cur] if x != prev][0]
            cyc.append(nxt); prev, cur = cur, nxt
        faces.append([(eidx[tuple(sorted((cyc[j], cyc[(j+1) % len(cyc)])))],
                       cyc[j], cyc[(j+1) % len(cyc)]) for j in range(len(cyc))])
    return verts, edges, faces

def summand_dim_from(verts, edges, faces, field=QQ):
    rows = []
    for cyc in faces:
        for c in range(4):
            row = [field(0)] * len(edges)
            for (e, a, b) in cyc:
                row[e] += (verts[b] - verts[a])[c]
            rows.append(row)
    return matrix(field, rows).right_kernel().dimension()

def summand_dim(Q):
    return summand_dim_from(*cell_combinatorics(Q))

print("== routine validated on polytopes with known summand spaces ==")
for name, P, want in [
        ("tetrahedron (indecomposable)",
         Polyhedron(vertices=[(0,0,0,0),(1,0,0,0),(0,1,0,0),(0,0,1,0)]), 1),
        ("cube = three segments",
         Polyhedron(vertices=[(a,b,c,0) for a in (0,1) for b in (0,1) for c in (0,1)]), 3),
        ("triangular prism = triangle + segment",
         Polyhedron(vertices=[(0,0,0,0),(1,0,0,0),(0,1,0,0),(0,0,1,0),(1,0,1,0),(0,1,1,0)]), 2)]:
    ok(f"{name}: summand dim {summand_dim(P)} = {want}", summand_dim(P) == want)

print("\n== agreement with ../paper4/necessity/staircase.py on the lattice facets ==")
for fi, want, note in [(21, 1, "pentagon facet, rigid"), (24, 1, "pentagon facet, rigid"),
                       (12, 2, "all-squares facet, splits all four squares")]:
    P = Polyhedron(vertices=[list(V[i]) for i in sorted(facs[fi][2])], base_ring=QQ)
    d = summand_dim(P)
    ok(f"facet {fi}: summand dim {d} ({note})", d == want)

print("\n== the slice cell for every degree, by convex hull and symbolically ==")
for fi, others in [(21, [3, 9, 10, 13]), (24, [4, 8, 11, 12])]:
    shape = None
    for s in [QQ(1), QQ((1,2)), QQ((1,3)), QQ((1,7)), QQ((1,50)), QQ((1,1000)),
              QQ((2,3)), QQ((3,7))]:
        Q = Polyhedron(vertices=[list(s * V[i]) for i in others]
                       + [list(V[j]) for j in PENT], base_ring=QQ)
        vv, ee, ff = cell_combinatorics(Q)
        sig = (Q.n_vertices(), len(ee), tuple(sorted(len(c) for c in ff)))
        shape = shape or sig
        d = summand_dim_from(vv, ee, ff)
        ok(f"facet {fi}, s = {s}: type {sig} unchanged, summand dim {d}",
           d == 1 and sig == shape)
    Fs = FractionField(PolynomialRing(QQ, 's')); s = Fs.gen()
    Q1 = Polyhedron(vertices=[list(QQ((1,2)) * V[i]) for i in others]
                    + [list(V[j]) for j in PENT], base_ring=QQ)
    vv, ee, ff = cell_combinatorics(Q1)
    lookup = {}
    for i in others: lookup[tuple(QQ((1,2)) * V[i])] = ('s', i)
    for j in PENT:   lookup[tuple(V[j])] = ('1', j)
    sym = []
    for v in vv:
        kind, i = lookup[tuple(v)]
        sym.append(vector(Fs, [s * x for x in V[i]]) if kind == 's'
                   else vector(Fs, [Fs(x) for x in V[i]]))
    d = summand_dim_from(sym, ee, ff, field=Fs)
    ok(f"facet {fi}: summand dim over Q(s) is {d}, i.e. rigid at every degree", d == 1)

print(f"\n{CH[0]} checks passed.")
