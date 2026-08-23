"""h^0 of the anticanonical divisor on a complexity-one T-variety over P^1.

Petersen-Suess, Israel J. Math. 182 (2011):
  Thm 3.21   K_X = sum_{(P,v)} (mu(v) K_Y(P) + mu(v) - 1) D_{(P,v)} - sum_rho D_rho
  Cor 3.19   D_h = -sum_rho h_t(n_rho) D_rho - sum_{(P,v)} mu(v) h_P(v) D_{(P,v)}
  Def 3.22   Box_h = {u : <u,v> >= h_t(v) for all v in |tail S|},
             h^*(u) = sum_P minvert(u - h_P) . P
  Prop 3.23  L(h)_u = Gamma(Y, O(h^*(u)))         for u in Box_h
  Def 3.8    h is Cartier iff h|_D is principal for every D with complete locus
             -- on P^1, sum_P a_P(sigma) = 0.

Reading Cor 3.19 backwards for D = -K_X gives the vertex values directly:
    h_P(v) = K_Y(P) + 1 - 1/mu(v).
The linear part h_t is then solved for, one maximal tail cone at a time.

h^0(X,-K) = sum over u in Box cap M of max(0, deg floor(h^*(u)) + 1).

Validated on PS Example 3.25 (X = P(Omega_{P^2}), h^0 = 27) before use.
"""
import itertools

def hzero(name, tail_max, slices, KY, expect=None, verbose=False):
    """tail_max: list of maximal tail cones, each a list of primitive ray gens.
       slices:   dict P -> list of cells; a cell is (tailcone_or_None, [vertices]).
                 A cell whose tail is a maximal cone must be tagged with its
                 index in tail_max; bounded / lower-tail cells are tagged None.
       KY:       dict P -> coefficient of the chosen canonical divisor on P^1."""
    pts = sorted(slices)
    n = len(tail_max[0][0])
    def mu(v):
        return lcm([QQ(x).denominator() for x in v])
    # h_P at every vertex of every slice
    hval = {}
    for P in pts:
        for _, vs in slices[P]:
            for v in vs:
                hval[(P, tuple(v))] = QQ(KY[P]) + 1 - QQ(1) / mu(v)
    # solve for u_tau on each maximal tail cone
    U, XR = {}, {}
    for ti, tau in enumerate(tail_max):
        rows, rhs = [], []
        cellP = {}
        for P in pts:
            cells = [c for c in slices[P] if c[0] == ti]
            assert len(cells) == 1, (name, P, ti, len(cells))
            cellP[P] = cells[0][1]
            for v in cells[0][1]:
                rows.append(list(v) + [1 if Q == P else 0 for Q in pts])
                rhs.append(hval[(P, tuple(v))])
        rows.append([0] * n + [1] * len(pts))       # principal: sum a_P = 0
        rhs.append(0)
        # PS Def 3.18: rho is EXTREMAL iff it misses deg D = sum_P D_P.
        # Only extremal rays carry a horizontal divisor, and Cor 3.19 then
        # forces h_t(n_rho) = -coef_rho(-K) = -1.
        tailcone = Polyhedron(rays=[list(r) for r in tau], base_ring=QQ)
        degD = Polyhedron(vertices=[[0] * n], base_ring=QQ)
        for P in pts:
            degD = degD + Polyhedron(vertices=[list(v) for v in cellP[P]],
                                     rays=[list(r) for r in tau], base_ring=QQ)
        xr = []
        for r in tau:
            ray = Polyhedron(rays=[list(r)], base_ring=QQ)
            if (ray & degD).is_empty():
                xr.append(r)
                rows.append(list(r) + [0] * len(pts))
                rhs.append(-1)
        XR[ti] = xr
        A = matrix(QQ, rows); b = vector(QQ, rhs)
        sol = A.solve_right(b)
        assert A * sol == b, (name, ti, "inconsistent: not Cartier")
        U[ti] = vector(QQ, sol[:n])
    # h_t on the rays, and Box
    ht = {}
    for ti, tau in enumerate(tail_max):
        for r in tau:
            val = U[ti].dot_product(vector(QQ, r))
            if tuple(r) in ht:
                assert ht[tuple(r)] == val, (name, r, "h_t not well defined")
            ht[tuple(r)] = val
    Box = Polyhedron(ieqs=[[-ht[r]] + list(r) for r in ht], base_ring=QQ)
    pts_box = Box.integral_points()
    total = 0
    for u in pts_box:
        uu = vector(QQ, u)
        deg = 0
        for P in pts:
            m = min(uu.dot_product(vector(QQ, v)) - hval[(P, tuple(v))]
                    for _, vs in slices[P] for v in vs)
            deg += floor(m)
        total += max(0, deg + 1)
    nx = len(set(tuple(r) for v in XR.values() for r in v))
    print(f"  {name}: tail fan {len(tail_max)} maximal cones, {len(ht)} rays "
          f"({nx} extremal); "
          f"Box has {len(pts_box)} lattice points;  h^0(-K) = {total}"
          + ("" if expect is None else f"   (expected {expect}) "
             + ("OK" if total == expect else "*** MISMATCH ***")))
    if expect is not None:
        assert total == expect, name
    return total, ht, hval

# ---------------------------------------------------------- validation
print("== validation: Petersen-Suess Example 3.25, X = P(Omega_{P^2}) ==")
rays = [(0,1),(1,1),(1,0),(0,-1),(-1,-1),(-1,0)]
sig = [[(0,1),(1,1)], [(0,1),(-1,0)], [(-1,0),(-1,-1)],
       [(-1,-1),(0,-1)], [(0,-1),(1,0)], [(1,0),(1,1)]]
S0 = [(0,[(0,1)]), (1,[(0,1)]), (2,[(0,0),(0,1)]),
      (3,[(0,0)]), (4,[(0,0)]), (5,[(0,0),(0,1)])]
Si = [(0,[(0,0)]), (1,[(0,0),(-1,-1)]), (2,[(-1,-1)]),
      (3,[(-1,-1)]), (4,[(0,0),(-1,-1)]), (5,[(0,0)])]
S1 = [(0,[(0,0),(1,0)]), (1,[(0,0)]), (2,[(0,0)]),
      (3,[(0,0),(1,0)]), (4,[(1,0)]), (5,[(1,0)])]
hzero("P(Omega_P2)", sig, {"0": S0, "oo": Si, "1": S1},
      {"0": -2, "oo": 0, "1": 0}, expect=27)

# --------------------------------------------------- the Delta_9 family
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper4", "certificates"))
from batyrev_global import facets, two_faces

print("\n== Delta_9: special fibre (toric) and general fibre ==")
W = json.load(open(os.path.join(HERE, "v09_candidate.json")))
Vt = [tuple(int(x) for x in v) for v in W["V"]]
V = [vector(ZZ, v) for v in Vt]
facs = facets(Vt)
FI = 11
FSET = frozenset(facs[FI][2])
R = vector(ZZ, [-x for x in facs[FI][0]])
Tdil = (0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0)
root = sorted(FSET)[0]
w = None
for cand in [vector(ZZ, e) for e in identity_matrix(ZZ, 4).rows()]:
    if R.dot_product(cand) == 1:
        w = cand; break
assert w is not None
K = matrix(ZZ, R.row().right_kernel().basis_matrix())      # basis of R^perp
def coords(x):
    return vector(QQ, K.transpose().solve_right(vector(QQ, x)))
def s(x):
    x = vector(QQ, x)
    return x - R.dot_product(x) * vector(QQ, w)

# phi_0, phi_1 on the facet
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
        if b in phi0:
            assert phi0[b] == val
        else:
            phi0[b] = val; st.append(b)
phi1 = {a: (V[a] - V[root]) - phi0[a] for a in phi0}
shift = vector(QQ, w) + vector(QQ, V[root])           # keeps the slices summing right

# cones of the face fan and their tails
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
Hperp = Polyhedron(eqns=[[0] + list(R)], base_ring=QQ)
tails, tail_of = [], {}
for S in cones:
    C = Polyhedron(rays=[list(V[i]) for i in sorted(S)], base_ring=QQ) & Hperp
    if C.dim() == 3:
        rr = tuple(sorted(tuple(coords(r)) for r in C.rays()))
        if rr not in tail_of:
            tail_of[rr] = len(tails); tails.append([list(x) for x in rr])
        tail_of.setdefault(("cone", S), tail_of[rr])
tail_max = [[tuple(QQ(y) for y in x) for x in t] for t in tails]
print(f"    tail fan: {len(tail_max)} maximal cones")

def build(level, mapper):
    out = []
    Hl = Polyhedron(eqns=[[-level] + list(R)], base_ring=QQ)
    for S in cones:
        C = Polyhedron(rays=[list(V[i]) for i in sorted(S)], base_ring=QQ) & Hl
        if C.is_empty():
            continue
        Ct = Polyhedron(rays=[list(V[i]) for i in sorted(S)], base_ring=QQ) & Hperp
        tag = None
        if Ct.dim() == 3:
            tag = tail_of[tuple(sorted(tuple(coords(r)) for r in Ct.rays()))]
        vs = [tuple(mapper(S, vector(QQ, v))) for v in C.vertices_list()]
        if vs:
            out.append((tag, vs))
    return out

# special fibre: slices s(sigma cap [R = +-1]) at 0 and oo
S_0 = build(1, lambda S, v: coords(s(v)))
S_inf = build(-1, lambda S, v: coords(s(v)))
def dedupe(cells):
    best = {}
    for tag, vs in cells:
        if tag is None:
            continue
        best[tag] = vs
    return [(t, best[t]) for t in sorted(best)] + \
           [(None, vs) for tag, vs in cells if tag is None]
hzero("Delta_9 special fibre", tail_max,
      {"0": dedupe(S_0), "oo": dedupe(S_inf)}, {"0": -1, "oo": -1}, expect=162)

# general fibre: slice at 0 unchanged; the level-(-1) slice splits into the
# two summand complexes, placed at oo and at p_lambda.  The shift keeps
# (slice at oo) + (slice at p) = s(level -1 slice), so that the family
# degenerates to the special fibre as p_lambda -> oo.
def summand(which):
    out = []
    Hm = Polyhedron(eqns=[[1] + list(R)], base_ring=QQ)
    for S in cones:
        C = Polyhedron(rays=[list(V[i]) for i in sorted(S)], base_ring=QQ) & Hm
        if C.is_empty():
            continue
        Ct = Polyhedron(rays=[list(V[i]) for i in sorted(S)], base_ring=QQ) & Hperp
        tag = None
        if Ct.dim() == 3:
            tag = tail_of[tuple(sorted(tuple(coords(r)) for r in Ct.rays()))]
        onF = sorted(S & FSET)
        base = phi0 if which == 0 else phi1
        vs = [tuple(coords(base[a] + (shift if which == 0 else 0))) for a in onF]
        if vs:
            out.append((tag, sorted(set(vs))))
    return out

S_a, S_b = dedupe(summand(0)), dedupe(summand(1))
print()
for KYgen in ({"0": -1, "oo": -1, "p": 0},
              {"0": -2, "oo": 0, "p": 0},
              {"0": 0, "oo": -1, "p": -1}):
    hzero(f"Delta_9 GENERAL fibre, K_P1 = {KYgen}", tail_max,
          {"0": dedupe(S_0), "oo": S_a, "p": S_b}, KYgen)

# ---------------------------------------------------------------- Piece 3
# PS Thm 3.27: D_h is semiample iff every h_P is concave and -h|_sigma(0) is
# semiample.  We imposed sum_P a_P(sigma) = 0 (Cartier, PS Def 3.8), so
# h|_sigma(0) has degree 0 on P^1, hence is principal -- the second condition
# holds automatically.  So semiampleness reduces to concavity of each h_P, and
# Ilten-Suess (0910.5919) Section 5 then gives global generation, because
# Y = P^1.
print("\n== concavity of h_P, i.e. semiampleness (PS Thm 3.27) ==")

def cells_full(level, restrict=None):
    """3-dimensional cells of a slice: (vertices, tail rays), in N' coords."""
    out = []
    Hl = Polyhedron(eqns=[[-level] + list(R)], base_ring=QQ)
    for S in cones:
        C = Polyhedron(rays=[list(V[i]) for i in sorted(S)], base_ring=QQ) & Hl
        if C.is_empty() or C.dim() != 3:
            continue
        Ct = Polyhedron(rays=[list(V[i]) for i in sorted(S)], base_ring=QQ) & Hperp
        if restrict is None:
            vs = sorted(set(tuple(coords(s(vector(QQ, v)))) for v in C.vertices_list()))
        else:
            onF = sorted(S & FSET)
            base = phi0 if restrict == 0 else phi1
            vs = sorted(set(tuple(coords(base[a] + (shift if restrict == 0 else 0)))
                            for a in onF))
        rs = sorted(set(tuple(coords(r)) for r in Ct.rays()))
        # keep only cells that are still 3-dimensional AFTER taking the summand
        cell = Polyhedron(vertices=[list(v) for v in vs],
                          rays=[list(r) for r in rs], base_ring=QQ)
        if cell.dim() != 3:
            continue
        out.append((vs, rs))
    return out

def concave(label, cells, KYP, ht):
    """Build h_P piece by piece -- vertices give the values, the tail rays give
    the linear part -- then test h_P <= every affine piece at every vertex."""
    val = {}
    for vs, _ in cells:
        for v in vs:
            val[v] = QQ(KYP) + 1 - QQ(1) / lcm([QQ(x).denominator() for x in v])
    pieces, skipped = [], 0
    for vs, rs in cells:
        rows = [list(v) + [1] for v in vs] + [list(r) + [0] for r in rs]
        rhs = [val[v] for v in vs] + [ht[r] for r in rs]
        A = matrix(QQ, rows); b = vector(QQ, rhs)
        if A.rank() < 4:
            skipped += 1
            continue
        sol = A.solve_right(b)
        if A * sol != b:
            return False, 0, 0
        pieces.append(sol)
    bad = 0
    for sol in pieces:
        for v in val:
            if sum(sol[i] * QQ(v[i]) for i in range(3)) + sol[3] < val[v]:
                bad += 1
    print(f"    {label}: {len(pieces)} affine pieces determined, {skipped} "
          f"under-determined; {'CONCAVE' if bad == 0 else str(bad) + ' violations'}")
    return bad == 0, len(pieces), skipped

_, ht_sp, _ = hzero("(recompute special)", tail_max,
                    {"0": dedupe(S_0), "oo": dedupe(S_inf)}, {"0": -1, "oo": -1})
_, ht_gen, _ = hzero("(recompute general)", tail_max,
                     {"0": dedupe(S_0), "oo": S_a, "p": S_b},
                     {"0": -1, "oo": -1, "p": 0})
allc, tot_sk = True, 0
for lab, cl, ky, ht in [
        ("special, slice 0 ", cells_full(1), -1, ht_sp),
        ("special, slice oo", cells_full(-1), -1, ht_sp),
        ("general, slice 0 ", cells_full(1), -1, ht_gen),
        ("general, slice oo", cells_full(-1, 0), -1, ht_gen),
        ("general, slice p ", cells_full(-1, 1), 0, ht_gen)]:
    okc, npc, sk = concave(lab, cl, ky, ht)
    allc = allc and okc; tot_sk += sk
if allc and tot_sk == 0:
    print("    => every h_P is concave, so -K is SEMIAMPLE on both fibres;")
    print("       on Y = P^1 semiample implies globally generated "
          "(Ilten-Suess Sec. 5), which is Piece 3.")
else:
    print(f"    => NOT established ({tot_sk} cells under-determined)")
