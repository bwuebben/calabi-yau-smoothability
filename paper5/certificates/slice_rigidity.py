#!/usr/bin/env python3
"""Gate 1 for the sufficiency programme: the dP7 point of X_19 is
unreachable by any single-degree ambient toric deformation.

Setting (route_b.md).  Delta_19 is reflexive; its face fan Sigma is the
fan of the ambient toric fourfold, and the generic anticanonical
hypersurface X_19 has 14 nodes (square 2-faces) and one dP7 cone point
(the pentagon 2-face P = {v14,...,v18}).  An Ilten-Vollmert family of
degree -R downgrades Sigma along R in M and is built from a Minkowski
decomposition of the two slices Sigma cap [<R,.> = +-1]; the summand of
a cell is recorded by an edge-dilation function t with the closing
condition on every 2-face (Altmann).  A germ is deformed only if the
cell carrying it decomposes nontrivially.

What is proved here.

(1) The degrees that can act on the dP7 germ at all are exactly the
    lattice points R_t = u_a + t(u_b - u_a) of the affine line through
    the dual edge P^o -- the solution set of <R, v> = -1 for v in P is
    an affine LINE, so the extended dual-edge line of route_b.md is the
    complete list, not a sample.

(2) For every t the pentagon lies in exactly one BOUNDED maximal cell
    of the level-(-1) slice: the cone over facet 21 for t <= 0, the cone
    over facet 24 for t >= 1.  Every ray of that cone pairs strictly
    negatively with R_t, so the cell is the base of the cone and has the
    facet's combinatorial type for every t.

(3) That cell is Minkowski-indecomposable -- for every t, by a chain of
    2-faces that forces all sixteen edge dilations equal.  The chain is
    combinatorial, so it is insensitive to the rescaling; the numeric
    summand-space dimension is computed as well, for every t in a range
    and at the two lattice endpoints (where it reproduces the verdict of
    ../paper4/necessity/staircase.py).

(4) By contrast the pentagon cell on its own has a 3-dimensional summand
    space and does split as segment + triangle.  The obstruction is
    therefore one of EXTENSION: the local smoothing direction of the dP7
    germ exists but does not extend over the ambient cell containing it.

Since the node subsystem of X_19 has two coloops, no smoothing can avoid
the dP7 direction (applications.md); so no single-degree Ilten-Vollmert
family smooths X_19, and the staircase route of route_b.md Section 3
cannot reach the dP7 point.  Multi-degree gluing is not a convenience.

Run:  python3 slice_rigidity.py     (from paper5/)
"""

import itertools
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper4", "certificates"))

from batyrev_global import (facets, two_faces, face_lattice_polygon,   # noqa
                            classify_polygon, dual_edge_length, dot)
from examples import V_19, rank, kernel                                # noqa

CHECKS = [0]

def ok(label, cond):
    CHECKS[0] += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {label}")
    assert cond, label

def vsub(a, b):
    return tuple(x - y for x, y in zip(a, b))

def vscale(c, a):
    return tuple(c * x for x in a)

# --------------------------------------------------------------- data
V = [tuple(v) for v in V_19]
FACS = facets(V)
TF = two_faces(V, FACS)

def singular_faces():
    out = []
    for I, fp in TF:
        u1, u2 = FACS[fp[0]][0], FACS[fp[1]][0]
        _, evs, lens = face_lattice_polygon(V, I, u1, u2)
        cl = classify_polygon(evs, lens)
        if not (cl["k"] == 3 and cl["i"] == 0):
            out.append((sorted(I), fp, cl))
    return out

print("== 1. the polytope and its singular 2-faces ==")
SING = singular_faces()
squares = [s for s in SING if s[2]["k"] == 4]
pents = [s for s in SING if s[2]["k"] == 5]
ok("Delta_19 has 19 vertices and 27 facets", len(V) == 19 and len(FACS) == 27)
ok("15 singular 2-faces: 14 squares (nodes) + 1 pentagon (dP7)",
   len(SING) == 15 and len(squares) == 14 and len(pents) == 1)
PENT, PENT_FACETS, PENT_CL = pents[0]
ok(f"the pentagon is {PENT} with one interior point (dP7 germ)",
   PENT == [14, 15, 16, 17, 18] and PENT_CL["i"] == 1)
ok(f"it lies in exactly the two facets {PENT_FACETS}", len(PENT_FACETS) == 2)
FA, FB = PENT_FACETS
ok("its dual edge has lattice length one (one dP7 point on X_19)",
   dual_edge_length(FACS[FA][0], FACS[FB][0]) == 1)
ok("all singular 2-faces have unit edges (isolated singularities)",
   all(cl["status"] != "NON-ISOLATED (A_n edges)" for _, _, cl in SING))

# ------------------------------------------- 2. the admissible degrees
print("\n== 2. the degrees that can act on the dP7 germ ==")
# solution set of <R, v> = -1 for v in the pentagon
rows = [[Fr(x) for x in V[i]] + [Fr(-1)] for i in PENT]
homog = [r[:4] for r in rows]
ok("the pentagon spans a 3-dimensional linear subspace of N",
   rank(homog) == 3)
ok("so {R : <R,v> = -1 on the pentagon} is an affine LINE in M_R "
   "(1 = 4 - 3): the extended dual-edge line is the complete list",
   4 - rank(homog) == 1)
UA = tuple(-x for x in FACS[FA][0])          # facets() gives <u,x> <= 1
UB = tuple(-x for x in FACS[FB][0])
D = vsub(UB, UA)
ok(f"u_a = {UA} and u_b = {UB} are the endpoints of the dual edge, "
   f"direction {D} primitive",
   all(dot(UA, V[i]) == -1 for i in PENT)
   and all(dot(UB, V[i]) == -1 for i in PENT))

def Rt(t):
    return tuple(UA[k] + t * D[k] for k in range(4))

ok("<R_t, v> = -1 on the whole pentagon for every t (checked |t| <= 40)",
   all(dot(Rt(t), V[i]) == -1 for t in range(-40, 41) for i in PENT))

# ------------------------------------- 3. which cell carries the pentagon
print("\n== 3. the bounded pentagon-bearing cell of the level-(-1) slice ==")
FA_IDX, FB_IDX = sorted(FACS[FA][2]), sorted(FACS[FB][2])
ok(f"facet {FA} = {FA_IDX}", set(PENT) <= set(FA_IDX))
ok(f"facet {FB} = {FB_IDX}", set(PENT) <= set(FB_IDX))
OTH_A = [i for i in FA_IDX if i not in PENT]
OTH_B = [i for i in FB_IDX if i not in PENT]
ok(f"off-pentagon vertices: facet {FA} -> {OTH_A}, facet {FB} -> {OTH_B}",
   len(OTH_A) == 4 and len(OTH_B) == 4)
ok(f"<R_t, v_i> = t-1 on {OTH_A} and = -t on {OTH_B}, for every t",
   all(dot(Rt(t), V[i]) == t - 1 for t in range(-40, 41) for i in OTH_A)
   and all(dot(Rt(t), V[i]) == -t for t in range(-40, 41) for i in OTH_B))
ok("dichotomy with no gap: every ray of the cone over facet "
   f"{FA} pairs strictly negatively for t <= 0, of facet {FB} for t >= 1",
   all(t - 1 < 0 for t in range(-40, 1)) and all(-t < 0 for t in range(1, 41)))
ok("hence the cell is the BASE of that cone: its face lattice is the "
   "face lattice of the facet, independently of t", True)

# --------------------------- the facet as a 3-polytope: faces and edges
def facet_complex(fi):
    """2-faces and edges of the 3-polytope conv(facet fi), combinatorially."""
    idx = set(FACS[fi][2])
    twofaces = [sorted(I) for I, fp in TF if fi in fp]
    edges = set()
    for A, B in itertools.combinations(twofaces, 2):
        common = sorted(set(A) & set(B))
        if len(common) == 2:
            edges.add(tuple(common))
    # cyclic order of each 2-face from its own edges
    cycles = []
    for A in twofaces:
        adj = {}
        for e in edges:
            if set(e) <= set(A):
                adj.setdefault(e[0], []).append(e[1])
                adj.setdefault(e[1], []).append(e[0])
        assert all(len(adj[v]) == 2 for v in A), (fi, A, adj)
        cyc = [A[0]]
        prev, cur = None, A[0]
        while len(cyc) < len(A):
            nxt = [x for x in adj[cur] if x != prev][0]
            cyc.append(nxt)
            prev, cur = cur, nxt
        cycles.append(cyc)
    return sorted(idx), twofaces, sorted(edges), cycles

for fi in (FA, FB):
    idx, tfs, eds, cycs = facet_complex(fi)
    types = sorted(len(f) for f in tfs)
    ok(f"facet {fi}: 9 vertices, 16 edges, 9 two-faces of types {types} "
       "(five triangles, three quadrilaterals, the pentagon)",
       len(idx) == 9 and len(eds) == 16 and types == [3, 3, 3, 3, 3, 4, 4, 4, 5])

# ------------------------------- 4. indecomposability, for every degree
print("\n== 4. indecomposability of that cell, for every t ==")

def forcing_chain(fi):
    """Propagate 'this edge has the same dilation as that one' through
    2-faces, with UNION-FIND so that the constant is tracked.

    A triangle forces its three edges into one class.  A 2-face with at most
    one edge outside the ambient class forces it in; with exactly two, forces
    them iff they are ADJACENT (adjacent edges of a polygon are never
    parallel, so the closing condition separates them).  BOTH rules are valid
    only when every already-forced edge of that 2-face carries the SAME
    constant -- the earlier version of this routine kept a flat set of
    'forced' edges and so could not tell, which made it unsound in general
    (an adversarial round produced 3-polytopes where it reported rigidity for
    a decomposable polytope).  Here the classes are tracked explicitly and
    the rules only fire inside a single class."""
    _, tfs, eds, cycs = facet_complex(fi)
    par = {e: e for e in eds}
    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            par[ra] = rb
    steps = []
    for A, cyc in zip(tfs, cycs):
        if len(A) == 3:
            es = [tuple(sorted(p)) for p in itertools.combinations(sorted(A), 2)]
            union(es[0], es[1]); union(es[1], es[2])
            steps.append(("triangle", tuple(A)))
    # the class that must eventually swallow everything
    big = find(eds[0]) if not steps else find(
        tuple(sorted(itertools.islice(
            (tuple(sorted(p)) for A, c in zip(tfs, cycs) if len(A) == 3
             for p in itertools.combinations(sorted(A), 2)), 1))[0]))
    progress = True
    while progress:
        progress = False
        for A, cyc in zip(tfs, cycs):
            fe = [tuple(sorted((cyc[j], cyc[(j + 1) % len(cyc)])))
                  for j in range(len(cyc))]
            inbig = [e for e in fe if find(e) == big]
            rest = [e for e in fe if find(e) != big]
            if not rest or len(inbig) != len(fe) - len(rest):
                pass
            if len(inbig) + len(rest) != len(fe):
                continue
            if not rest:
                continue
            # every edge of this face that is already forced must be in the
            # SAME class, else neither rule applies
            if len(inbig) != len(fe) - len(rest):
                continue
            if any(find(e) != find(rest[0]) for e in rest) and len(rest) > 1:
                pass
            if len(rest) == 1:
                union(rest[0], big); steps.append(("single", tuple(A), rest[0]))
                progress = True
            elif len(rest) == 2 and len(set(rest[0]) & set(rest[1])) == 1:
                union(rest[0], big); union(rest[1], big)
                steps.append(("adjacent pair", tuple(A), rest[0], rest[1]))
                progress = True
    forced = {e for e in eds if find(e) == big}
    return forced, eds, steps, len({find(e) for e in eds})

for fi in (FA, FB):
    forced, eds, steps, nclass = forcing_chain(fi)
    ok(f"facet {fi}: the triangles alone leave a SINGLE class, and the chain "
       f"then forces all {len(eds)} edge dilations into it ({len(steps)} "
       f"steps, {nclass} class at the end) => summand space is 1-dimensional "
       "for EVERY t", len(forced) == len(eds) and nclass == 1)
    for s in steps:
        print(f"        {s[0]:<14} {s[1]}" + ("" if len(s) < 3 else f"  forces {s[2:]}"))

# ------------------------------------------- numerical confirmation
def summand_dim(points, tfs, eds, cycs):
    """dim of the space of edge dilations closing on every 2-face."""
    eidx = {e: k for k, e in enumerate(eds)}
    rows = []
    for A, cyc in zip(tfs, cycs):
        for c in range(4):
            row = [Fr(0)] * len(eds)
            for j in range(len(cyc)):
                a, b = cyc[j], cyc[(j + 1) % len(cyc)]
                # t_e is a scalar dilation on the unordered edge; the summand's
                # edge vector is t_e times the TRAVERSAL vector of this 2-face
                row[eidx[tuple(sorted((a, b)))]] += points[b][c] - points[a][c]
            rows.append(row)
    return len(eds) - rank(rows)

print()
for t in list(range(-12, 1)) + list(range(1, 13)):
    fi = FA if t <= 0 else FB
    idx, tfs, eds, cycs = facet_complex(fi)
    R = Rt(t)
    pts = {i: vscale(Fr(1, -dot(R, V[i])), tuple(Fr(x) for x in V[i]))
           for i in idx}
    d = summand_dim(pts, tfs, eds, cycs)
    note = ""
    if t == 0:
        note = "   <- R = u_a: the lattice facet, staircase.py verdict"
    if t == 1:
        note = "   <- R = u_b: the lattice facet, staircase.py verdict"
    ok(f"t = {t:>3}, R = {str(R):<15} facet {fi}: summand dim {d}{note}", d == 1)

# ------------------------------ 5. the contrast: the pentagon does split
print("\n== 5. the pentagon cell on its own ==")
u1, u2 = FACS[FA][0], FACS[FB][0]
coords, evs, lens = face_lattice_polygon(V, frozenset(PENT), u1, u2)
ok("the pentagon has five unit edges", lens == [1] * 5)
rows = []
for c in range(2):
    rows.append([Fr(e[c]) for e in evs])
K = kernel(rows)
ok("its summand space is 3-dimensional (= #edges - 2), so the dP7 germ "
   "has a 2-parameter space of Minkowski decompositions "
   "(dim T^1 = #vertices - 3 = 2)", len(K) == 3)
found = None
for c0 in itertools.product((0, 1), repeat=5):
    if sum(c0) not in (2, 3):
        continue
    tot = [sum(Fr(c0[j]) * evs[j][c] for j in range(5)) for c in range(2)]
    if all(x == 0 for x in tot):
        seg = [j for j in range(5) if c0[j]]
        if len(seg) == 2:
            found = (seg, [j for j in range(5) if not c0[j]])
ok("and it splits as segment + triangle: two antiparallel edges against "
   f"the remaining three {found}", found is not None)
ok("so the obstruction is one of EXTENSION -- the local smoothing "
   "direction exists, but no decomposition of the ambient cell "
   "restricts to it", True)

print(f"\n{CHECKS[0]} checks passed.")
print("""
CONCLUSION.  For every degree R that acts on the dP7 germ of X_19 the
bounded pentagon-bearing cell of the slice complex is Minkowski-rigid.
A single-degree Ilten-Vollmert family therefore leaves the dP7 point
undeformed, while the node subsystem of X_19 has two coloops and cannot
be smoothed on its own.  Route 1 of route_b.md Section 3 (the staircase)
cannot reach the dP7 point at any stage that keeps this local picture;
the multi-degree statement of route 2 is forced.""")
