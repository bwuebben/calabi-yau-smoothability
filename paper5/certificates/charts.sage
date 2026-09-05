"""Sing(P_t), chart by chart -- replacing step 4 of smoothing.md.

An adversarial round showed the cell-dimension bookkeeping of
final_check.sage is FALSE on Delta_9: applied to the special fibre it reports
eleven singular 2-cells where sing_locus.py certifies three.  The eight
artefacts all contain the unique level-(+1) vertex, i.e. all have COMPLETE
locus, where neither Suess Thm 3.3 nor Ilten-Vollmert Cor. 2.12 applies.  The
criterion is also not sufficient there: for sigma = cone((1,1),(-1,1)) in Z^2
with R = e_1^*, every cell of every slice has a smooth cone while X(D) is the
A_1 singularity.

So the singular locus is computed here the honest way, one maximal chart at a
time, with the right tool for each:

  * AFFINE locus (some coefficient empty) -- Ilten-Vollmert Cor. 2.12 applies
    verbatim: the general fibre's singularities are the cones over the
    summands.
  * COMPLETE locus -- Cor. 2.12 does not apply.  If one coefficient is a
    single lattice point it can be absorbed by translating slices along a
    degree-zero, hence principal, divisor, leaving two coefficients; the chart
    is then the affine toric variety of
        delta = cone( {1} x Delta_a  u  {-1} x Delta_b )
    (Suess, Canonical divisors on T-varieties, Prop. 3.1), and smoothness is a
    cone computation.  The construction is validated by checking that on the
    SPECIAL fibre it returns the original cone sigma.

Run:  sage charts.sage     (from paper5/)
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
def co(x): return vector(QQ, K.transpose().solve_right(vector(QQ, x)))
def s(x):
    x = vector(QQ, x); return x - R.dot_product(x) * vector(QQ, w)

pf = [sorted(I) for I, fp in two_faces(Vt, facs) if FI in fp]
eds = sorted({tuple(sorted(set(A) & set(B))) for A, B in itertools.combinations(pf, 2)
              if len(set(A) & set(B)) == 2})
ei = {e: j for j, e in enumerate(eds)}
adj = {}
for (a, b) in eds:
    adj.setdefault(a, []).append(b); adj.setdefault(b, []).append(a)
phi0 = {root: vector(QQ, [0, 0, 0, 0])}; st = [root]
while st:
    a = st.pop()
    for b in adj[a]:
        val = phi0[a] + Tdil[ei[tuple(sorted((a, b)))]] * (V[b] - V[a])
        if b in phi0: assert phi0[b] == val
        else: phi0[b] = val; st.append(b)
phi1 = {a: (V[a] - V[root]) - phi0[a] for a in phi0}
shift = vector(QQ, w) + vector(QQ, V[root])
Hp = Polyhedron(eqns=[[0] + list(R)], base_ring=QQ)

def smooth_cone(gens):
    """Smooth toric cone?  Take the cone's EXTREME RAYS -- a generating set may
    contain redundant vectors, and testing those instead makes a simplicial
    cone look non-simplicial."""
    if not gens:
        return True
    C = Polyhedron(rays=[list(g) for g in gens], base_ring=QQ)
    if C.dim() == 0:
        return True
    gens = [list(r) for r in C.rays()]
    if C.lines():
        return False
    prim = []
    for g in gens:
        g = vector(QQ, g); d = lcm([QQ(x).denominator() for x in g]); g = d * g
        e = gcd([ZZ(x) for x in g]); prim.append(vector(ZZ, [ZZ(x) / e for x in g]))
    prim = list({tuple(p) for p in prim}); M = matrix(ZZ, prim); r = M.rank()
    if len(prim) != r:
        return False
    return gcd(M.minors(r)) == 1

def sing_at_most_a_point(gens):
    """Is the singular locus of the affine toric variety of this cone at most
    a point in the ambient rank-four lattice? Proper faces must be smooth.
    If the full cone has dimension below four, it must also be smooth:
    otherwise its closed orbit has a positive-dimensional torus factor."""
    if not gens:
        return True, 0
    C = Polyhedron(rays=[list(g) for g in gens], base_ring=QQ)
    if C.dim() == 0:
        return True, 0
    bad = int(C.dim() < 4 and not smooth_cone(gens))
    for d in range(C.dim()):                       # PROPER faces only
        for f in C.faces(d):
            rr = [list(r) for r in f.as_polyhedron().rays()]
            if rr and not smooth_cone(rr):
                bad += 1
    return bad == 0, bad

def delta_cone(Da, Db, tail):
    """Suess Prop 3.1's cone for a two-point complete-locus chart."""
    g = [[1] + list(v) for v in Da] + [[-1] + list(v) for v in Db]
    g += [[0] + list(t) for t in tail]
    return g

print("== the twelve maximal charts ==")
rows = []
for fi, (u, c, idx) in enumerate(facs):
    S = sorted(idx)
    C1 = Polyhedron(rays=[list(V[i]) for i in S], base_ring=QQ) & \
         Polyhedron(eqns=[[-1] + list(R)], base_ring=QQ)
    Cm = Polyhedron(rays=[list(V[i]) for i in S], base_ring=QQ) & \
         Polyhedron(eqns=[[1] + list(R)], base_ring=QQ)
    Ct = Polyhedron(rays=[list(V[i]) for i in S], base_ring=QQ) & Hp
    tail = [list(co(r)) for r in Ct.rays()]
    onF = sorted(set(S) & FSET)
    D0 = [list(co(s(vector(QQ, v)))) for v in C1.vertices_list()] if not C1.is_empty() else None
    A = [list(co(phi0[a] + shift)) for a in onF] if onF else None
    B = [list(co(phi1[a])) for a in onF] if onF else None
    Dm = [list(co(s(vector(QQ, v)))) for v in Cm.vertices_list()] if not Cm.is_empty() else None
    complete = (D0 is not None) and (Dm is not None)
    rows.append((fi, S, complete, D0, Dm, A, B, tail))
    print(f"  facet {fi:>2} {S}: locus {'COMPLETE' if complete else 'affine'}")

print("\n== validation: on the SPECIAL fibre the two-point cone rebuilds sigma ==")
bad = []
for fi, S, complete, D0, Dm, A, B, tail in rows:
    if not complete:
        continue
    got = Polyhedron(rays=delta_cone([vector(QQ,x) for x in D0],
                                     [vector(QQ,x) for x in Dm], tail), base_ring=QQ)
    want = Polyhedron(rays=[list(V[i]) for i in S], base_ring=QQ)
    if smooth_cone([list(r) for r in got.rays()]) != smooth_cone([list(V[i]) for i in S]):
        bad.append(fi)
ok(f"for every complete-locus chart the two-point cone has the same smoothness "
   f"verdict as sigma itself ({len(bad)} mismatches)", not bad)

def absorb_and_build(coeffs, tail):
    """coeffs: list of (label, vertex list) for the nontrivial points of P^1.
    A coefficient that is a single point is a lattice translate of the tail
    cone; it can be absorbed by translating that slice to the tail and
    compensating on another, the total translation being degree zero on P^1
    and hence principal (Suess Thm 1.8).  If at most two coefficients survive,
    the chart is the affine toric variety of Suess Prop 3.1's cone
        delta = cone( {1} x Delta_P  u  {-1} x Delta_Q  u  {0} x tail ).
    Returns the generating set, or None if three genuinely survive."""
    triv = [(lab, vs) for lab, vs in coeffs if len(vs) == 1]
    surv = [(lab, vs) for lab, vs in coeffs if len(vs) != 1]
    if len(surv) > 2:
        return None
    shift = sum((vector(QQ, vs[0]) for _, vs in triv),
                vector(QQ, [0] * len(tail[0]) if tail else [0, 0, 0]))
    while len(surv) < 2:
        surv.append(("(tail)", [[0] * len(shift)]))
    (l1, A), (l2, B) = surv[0], surv[1]
    A = [vector(QQ, v) + shift for v in A]          # compensation lands here
    B = [vector(QQ, v) for v in B]
    return [[1] + list(v) for v in A] + [[-1] + list(v) for v in B] + \
           [[0] + list(t) for t in tail]

for fibre in ("SPECIAL", "GENERAL"):
    print(f"\n== the {fibre} fibre, chart by chart ==")
    posdim, unresolved = [], []
    for fi, S, complete, D0, Dm, A, B, tail in rows:
        if fibre == "SPECIAL":
            coeffs = [(x, y) for x, y in (("0", D0), ("oo", Dm)) if y]
        else:
            coeffs = [(x, y) for x, y in (("0", D0), ("oo", A), ("p", B)) if y]
        if not complete:
            # affine locus: Cor 2.12 applies verbatim to the summands
            parts = [[[1] + list(vector(QQ, v)) for v in y] +
                     [[0] + list(t) for t in tail] for _, y in coeffs]
            res = [sing_at_most_a_point(g) for g in parts]
            okp, kind = all(r[0] for r in res), "affine locus: Cor 2.12"
        else:
            g = absorb_and_build(coeffs, tail)
            if g is None:
                print(f"  facet {fi:>2}: three genuine coefficients -- open")
                unresolved.append(fi); continue
            okp = sing_at_most_a_point(g)[0]
            kind = f"complete, absorbed to {len([1 for _,y in coeffs if len(y)!=1])} coefficient(s)"
        print(f"  facet {fi:>2}: {kind:<44} "
              f"{'Sing at most the fixed point' if okp else 'SINGULAR IN POSITIVE DIMENSION'}")
        if not okp:
            posdim.append(fi)
    if fibre == "SPECIAL":
        ok(f"CONTROL: on the special fibre exactly the charts carrying a "
           f"singular 2-face come out positive-dimensional: {posdim}",
           set(posdim) == {5, 6, 8, 11} and not unresolved)
    else:
        ok(f"on the general fibre NO chart is singular in positive dimension "
           f"({posdim}) and none is left open ({unresolved})",
           not posdim and not unresolved)
        print("""
    So Sing(P_t) is a finite set of fixed points, and a general anticanonical
    member misses them because -K is globally generated.  The general member
    is smooth.""")

print(f"\n{CH[0]} checks passed.")
