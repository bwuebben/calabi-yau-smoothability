"""Inventory the cells changed by the Delta_9 deformation.

This script does *not* infer orbit-stratum dimension from cell dimension.  The
Süß affine-locus dictionary cannot be used that way on complete-locus charts;
that was defect D-H.  The global singular-locus statement is certified by
charts.sage.  Here we only check which cells split nontrivially and verify the
three bounded germ cells against Ilten--Vollmert Corollary 2.12.

Run:  sage final_check.sage     (from paper5/)
"""
import os, sys, json, itertools
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper4", "certificates"))
from batyrev_global import facets, two_faces

CH = [0]
def ok(label, cond):
    CH[0] += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {label}")
    assert cond, label

W = json.load(open(os.path.join(HERE, "v09_candidate.json")))
Vt = [tuple(int(x) for x in v) for v in W["V"]]
V = [vector(ZZ, v) for v in Vt]
facs = facets(Vt); FI = 11
FSET = frozenset(facs[FI][2])
R = vector(ZZ, [-x for x in facs[FI][0]])
Tdil = (0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0)
root = sorted(FSET)[0]
w = [vector(ZZ, e) for e in identity_matrix(ZZ, 4).rows()
     if R.dot_product(vector(ZZ, e)) == 1][0]
K = matrix(ZZ, R.row().right_kernel().basis_matrix())
def coords(x): return vector(QQ, K.transpose().solve_right(vector(QQ, x)))
def s(x):
    x = vector(QQ, x); return x - R.dot_product(x) * vector(QQ, w)

pf = [sorted(I) for I, fp in two_faces(Vt, facs) if FI in fp]
edges = sorted({tuple(sorted(set(A) & set(B))) for A, B in itertools.combinations(pf, 2)
                if len(set(A) & set(B)) == 2})
eidx = {e: j for j, e in enumerate(edges)}
adj = {}
for (a, b) in edges:
    adj.setdefault(a, []).append(b); adj.setdefault(b, []).append(a)
phi0 = {root: vector(QQ, [0, 0, 0, 0])}
st = [root]
while st:
    a = st.pop()
    for b in adj[a]:
        val = phi0[a] + Tdil[eidx[tuple(sorted((a, b)))]] * (V[b] - V[a])
        if b in phi0: assert phi0[b] == val
        else: phi0[b] = val; st.append(b)
phi1 = {a: (V[a] - V[root]) - phi0[a] for a in phi0}
shift = vector(QQ, w) + vector(QQ, V[root])

n = len(Vt)
seen, frontier = {frozenset(range(n))}, [frozenset(range(n))]
while frontier:
    nxt = []
    for S in frontier:
        for _, _, I in facs:
            Tt = S & I
            if Tt and Tt not in seen: seen.add(Tt); nxt.append(Tt)
    frontier = nxt
cones = [S for S in seen if S != frozenset(range(n))]
Hperp = Polyhedron(eqns=[[0] + list(R)], base_ring=QQ)

def cone_smooth(P):
    """is closure(Q_{>=0} . ({1} x P)) a smooth toric cone?"""
    gens = [vector(QQ, [1] + list(v)) for v in P.vertices_list()]
    gens += [vector(QQ, [0] + list(r)) for r in P.rays()]
    prim = []
    for g in gens:
        d = lcm([QQ(x).denominator() for x in g]); g = d * g
        e = gcd([ZZ(x) for x in g])
        prim.append(vector(ZZ, [ZZ(x) / e for x in g]))
    M = matrix(ZZ, prim); r = M.rank()
    if len(set(map(tuple, prim))) != r:
        return False                                    # not simplicial
    return gcd(M.minors(r)) == 1

def slice_cells(level, restrict=None):
    out = {}
    Hl = Polyhedron(eqns=[[-level] + list(R)], base_ring=QQ)
    for S in cones:
        C = Polyhedron(rays=[list(V[i]) for i in sorted(S)], base_ring=QQ) & Hl
        if C.is_empty(): continue
        Ct = Polyhedron(rays=[list(V[i]) for i in sorted(S)], base_ring=QQ) & Hperp
        rs = [list(coords(r)) for r in Ct.rays()]
        if restrict is None:
            vs = [list(coords(s(vector(QQ, v)))) for v in C.vertices_list()]
        else:
            base = phi0 if restrict == 0 else phi1
            vs = [list(coords(base[a] + (shift if restrict == 0 else 0)))
                  for a in sorted(S & FSET)]
        if not vs: continue
        cell = Polyhedron(vertices=vs, rays=rs, base_ring=QQ)
        key = (tuple(sorted(map(tuple, cell.vertices_list()))),
               tuple(sorted(map(tuple, cell.rays_list()))))
        out[key] = (sorted(S), cell)
    return list(out.values())

def report(label, cells):
    bad = {1: [], 2: []}
    cnt = {0: 0, 1: 0, 2: 0, 3: 0}
    for S, C in cells:
        d = C.dim(); cnt[d] = cnt.get(d, 0) + 1
        if d in (1, 2) and not cone_smooth(C):
            bad[d].append((S, C.n_vertices(), C.n_rays()))
    print(f"    {label}: cells by dimension {cnt}; "
          f"singular cones on 1-cells {len(bad[1])}, on 2-cells {len(bad[2])}")
    for d in (1, 2):
        for b in bad[d]:
            print(f"        dim {d} cell from cone {b[0]}: "
                  f"{b[1]} vertices, {b[2]} rays -> SINGULAR")
    return bad

# ---------------------------------------------------------------------------
# The criterion above is only valid for p-divisors with AFFINE locus (Suess
# Thm 3.3).  A cell whose cone has rays on both sides of R^perp has COMPLETE
# locus, and there the germ is not the cone over one slice cell -- in the
# toric special fibre it is just the cone sigma itself.  So the test has to be
# split, and the useful question is narrower:
#
#   WHICH cells does the deformation actually change?
#
# A cell decomposes trivially when one summand is a single point p.  Then the
# coefficient {p} at p_lambda can be moved to {0} by translating slices along
# the divisor p([p_lambda] - [oo]), which has degree zero on P^1 and is
# therefore principal -- Suess Thm 1.8 -- so the p-divisor, and with it the
# germ, is unchanged.  Everything the family does is concentrated on the cells
# that split nontrivially.

print("== which cells does the deformation change? ==")
Hm = Polyhedron(eqns=[[1] + list(R)], base_ring=QQ)
changed, unchanged = [], 0
for S in cones:
    C = Polyhedron(rays=[list(V[i]) for i in sorted(S)], base_ring=QQ) & Hm
    if C.is_empty():
        continue
    onF = sorted(S & FSET)
    A = Polyhedron(vertices=[list(coords(phi0[a] + shift)) for a in onF], base_ring=QQ)
    B = Polyhedron(vertices=[list(coords(phi1[a])) for a in onF], base_ring=QQ)
    if A.n_vertices() == 1 or B.n_vertices() == 1:
        unchanged += 1
    else:
        Ct = Polyhedron(rays=[list(V[i]) for i in sorted(S)], base_ring=QQ) & Hperp
        rs = [list(coords(r)) for r in Ct.rays()]
        cell = Polyhedron(vertices=[list(coords(s(vector(QQ, v))))
                                    for v in C.vertices_list()],
                          rays=rs, base_ring=QQ)
        changed.append((sorted(S), cell.dim(), len(rs),
                        cone_smooth(Polyhedron(vertices=A.vertices_list(),
                                               rays=rs, base_ring=QQ)),
                        cone_smooth(Polyhedron(vertices=B.vertices_list(),
                                               rays=rs, base_ring=QQ))))
print(f"    {unchanged} cells keep a trivial decomposition -- germ unchanged")
print(f"    {len(changed)} cells split nontrivially:")
for S, d, nr, sa, sb in sorted(changed, key=lambda x: (x[1], x[0])):
    print(f"        cone {S}: cell dim {d}, {nr} tail rays, summand cones "
          f"{'smooth' if sa else 'SINGULAR'} / {'smooth' if sb else 'SINGULAR'}")

low = [(S, d, nr, sa, sb) for S, d, nr, sa, sb in changed if d <= 2]
SING = [[2, 3, 6, 7], [2, 4, 5, 7], [3, 4, 5, 6, 8]]
ok(f"the changed cells of dimension at most 2 are exactly the three germs of "
   f"X_9, {sorted(S for S, _, _, _, _ in low)}",
   sorted(S for S, _, _, _, _ in low) == sorted(SING))
ok("each of them is bounded, so its polyhedral divisor has affine locus and "
   "Ilten-Vollmert Corollary 2.12 applies", all(nr == 0 for _, _, nr, _, _ in low))
ok("and both summand cones are smooth at each of those three affine germ cells",
   all(sa and sb for _, _, _, sa, sb in low))

print("""
    No stratum conclusion is drawn for the four changed higher-dimensional
    cells.  They belong to complete-locus charts, and charts.sage performs the
    required maximal-chart singular-locus computation.""")

print(f"\n{CH[0]} checks passed.")
print("""
CONCLUSION.  Exactly seven cells change.  The three changed cells belonging to
the hypersurface germs have affine locus and smooth summand cones.  Use
charts.sage, not this inventory, for the global ambient singular locus.""")
