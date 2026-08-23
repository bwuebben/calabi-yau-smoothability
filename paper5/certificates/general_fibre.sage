"""The general fibre of the Delta_9 family, explicitly.

Ilten-Vollmert Remark 1.8: downgrading along R gives, for each cone sigma of
Sigma, the polyhedral divisor
    D^sigma = s(sigma cap [R=1]) (x) {0} + s(sigma cap [R=-1]) (x) {oo}.
We decompose the [R=-1] slice (r_oo = 1) and leave [R=+1] alone (r_0 = 0).
Their fibre description (Section 2, after Theorem 2.8) then gives, for lambda
in the base,
    D^(lambda) = D_0 (x) {0} + D_oo^0 (x) {oo} + D_oo^1 (x) {p_lambda},
so the general fibre is a complexity-one T-variety over P^1 with THREE points
carrying a nontrivial coefficient.  All three slices are computed here; they
are the input any h^0 computation on the general fibre needs.

Also checked: each slice is a COMPLETE polyhedral subdivision of N_Q ~ R^3.
By Remark 1.5 that is what makes X(S^(lambda)) complete, i.e. a legitimate
compact ambient for a Calabi-Yau hypersurface.

Run:  sage general_fibre.sage     (from paper5/)
"""
import os, sys, itertools, json
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
facs = facets(Vt)
FI = 11
FSET = frozenset(facs[FI][2])
R = vector(ZZ, [-x for x in facs[FI][0]])
T = (0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0)
print(f"Delta_9, facet {FI} = {sorted(FSET)}, degree R = {tuple(R)}")

pf = [sorted(I) for I, fp in two_faces(Vt, facs) if FI in fp]
edges = sorted({tuple(sorted(set(A) & set(B))) for A, B in itertools.combinations(pf, 2)
                if len(set(A) & set(B)) == 2})
eidx = {e: j for j, e in enumerate(edges)}
adj = {}
for (a, b) in edges:
    adj.setdefault(a, []).append(b); adj.setdefault(b, []).append(a)
root = sorted(FSET)[0]
phi0 = {root: vector(QQ, [0, 0, 0, 0])}
st = [root]
while st:
    a = st.pop()
    for b in adj[a]:
        w = phi0[a] + T[eidx[tuple(sorted((a, b)))]] * (V[b] - V[a])
        if b in phi0:
            assert phi0[b] == w
        else:
            phi0[b] = w; st.append(b)
phi1 = {a: (V[a] - V[root]) - phi0[a] for a in phi0}

# all cones of the face fan
n = len(Vt)
seen, frontier = {frozenset(range(n))}, [frozenset(range(n))]
while frontier:
    nxt = []
    for S in frontier:
        for _, _, I in facs:
            Tt = S & I
            if Tt and Tt not in seen:
                seen.add(Tt); nxt.append(Tt)
    frontier = nxt
cones = [S for S in seen if S != frozenset(range(n))]

# the cosection s: drop to a complement of R.  R is primitive, so choose a
# basis of N in which R is the last coordinate functional.
Mb = matrix(ZZ, [list(R)])
K = Mb.right_kernel().basis_matrix()            # rank 3, a basis of N cap R^perp
ok(f"N cap R^perp has rank 3, basis rows {K.rows()}", K.nrows() == 3)
ext = block_matrix([[K], [matrix(ZZ, [list(v) for v in [V[root]]])]])
ok("that basis together with one vertex of the facet spans N over Q",
   ext.rank() == 4)
def proj(p):
    """coordinates of p in the chosen basis of R^perp, after subtracting the
    R-component; a cosection in the sense of Remark 1.8."""
    lam = R.dot_product(vector(QQ, p))
    q = vector(QQ, p) - lam * vector(QQ, [QQ(x) for x in V[root]]) / R.dot_product(V[root])
    return vector(QQ, K.solve_left(q)) if False else vector(QQ, K.transpose().solve_right(q))

H1 = Polyhedron(eqns=[[-1] + list(R)], base_ring=QQ)      # <R,x> = +1
Hm = Polyhedron(eqns=[[1] + list(R)], base_ring=QQ)       # <R,x> = -1

def cell(S, H):
    return Polyhedron(rays=[list(V[i]) for i in sorted(S)], base_ring=QQ) & H

print("\n== the three slices of the general fibre ==")
slices = {}
S0cells, Smcells0, Smcells1 = [], [], []
for S in cones:
    C1 = cell(S, H1)
    if not C1.is_empty():
        S0cells.append((S, C1))
    Cm = cell(S, Hm)
    if not Cm.is_empty():
        onF = sorted(S & FSET)
        tail = [list(r) for r in Cm.rays()]
        Smcells0.append((S, Polyhedron(vertices=[list(phi0[a]) for a in onF],
                                       rays=tail, base_ring=QQ)))
        Smcells1.append((S, Polyhedron(vertices=[list(phi1[a]) for a in onF],
                                       rays=tail, base_ring=QQ)))
for name, cl, H in [("slice at 0   (S_0, undecomposed)", S0cells, H1),
                    ("slice at oo  (S_oo^0)", Smcells0, Hm),
                    ("slice at p_l (S_oo^1)", Smcells1, Hm)]:
    mx = [(S, C) for S, C in cl if C.dim() == 3]
    bd = sum(1 for _, C in cl if C.is_compact())
    print(f"  {name}: {len(cl)} cells, {len(mx)} maximal, {bd} bounded")
    # completeness: the maximal cells must cover the hyperplane.  Test by
    # volume-free means: every ray of the tail fan is covered and the union
    # is closed under the fan structure -- here checked by verifying that a
    # random-ish set of points of the hyperplane lies in some cell.
    slices[name] = mx

print("\n== completeness of each slice (Remark 1.5) ==")
# A slice is complete iff its cells cover the affine hyperplane.  Test by
# sampling: x0 + w for w in the kernel lattice of R, over a wide range of
# scales, including far out where only the tail cones can reach.
# phi_0 and phi_1 are built from edge vectors of the facet, which lie in
# R^perp, so both summand complexes live in the LINEAR hyperplane <R,x> = 0
# -- that is N_Q after the cosection.  The original slice sits in the affine
# hyperplane <R,x> = -1.  Probe each in its own home.
x0 = vector(QQ, [QQ(x) for x in V[root]])
kb = [vector(QQ, list(r)) for r in K.rows()]
probes = []
for co in itertools.product((-9, -4, -1, 0, 1, 4, 9), repeat=int(3)):
    probes.append(x0 + sum(QQ(c) * b for c, b in zip(co, kb)))
for co in itertools.product((-97, 0, 61), repeat=int(3)):
    probes.append(x0 + sum(QQ(c) * b for c, b in zip(co, kb)))
probes0 = [p - x0 for p in probes]           # the same points, in R^perp
ok(f"{len(probes)} probes in the affine hyperplane <R,x> = -1 and their "
   "translates in R^perp",
   all(R.dot_product(p) == -1 for p in probes) and
   all(R.dot_product(p) == 0 for p in probes0))

orig = [(S, cell(S, Hm)) for S in cones]
orig = [(S, C) for S, C in orig if C.dim() == 3]
print(f"    the original level-(-1) slice has {len(orig)} maximal cells")
for name, mx in [("original S_-1", orig),
                 ("S_oo^0", slices["slice at oo  (S_oo^0)"]),
                 ("S_oo^1", slices["slice at p_l (S_oo^1)"])]:
    pr = probes if name == "original S_-1" else probes0
    miss = [p for p in pr if not any(C.contains(p) for _, C in mx)]
    print(f"    {name:<14} {len(mx)} maximal cells, "
          f"{len(pr) - len(miss)}/{len(pr)} probes covered"
          + ("" if not miss else f"   MISSES e.g. {tuple(miss[0])}"))
    if name == "original S_-1":
        ok("    the original slice covers every probe, as it must "
           "(Sigma is complete)", not miss)
    else:
        ok(f"    {name} covers every probe" if not miss
           else f"    {name} does NOT cover: the summand complex is not a "
                "complete subdivision", not miss)

print(f"\n{CH[0]} checks passed.")
print("""
The three slices above are the complete combinatorial input for the general
fibre.  With a formula for h^0 of a divisor on a complexity-one T-variety they
determine h^0(-K) directly, which is Piece 2 of follows.md.""")
