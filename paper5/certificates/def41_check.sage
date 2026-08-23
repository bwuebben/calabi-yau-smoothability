"""Gate 1, completed against Ilten-Vollmert Definition 4.1 as written
(verbatim text from arXiv:0903.1393v3, source tarball).

Def. 2.1  : every summand carries the FULL tail cone of the polyhedron.
Def. 4.1(i) : (D cap E)^i = D^i cap E^i for cells D, E of the subdivision.
Def. 4.1(ii): sum_{i in I} cap_{D in I'} D^i  <  sum_{i in I} cap_{D in J} D^i,
              a FACE relation, for J subset I' subset C, I subset {0..r}.

With I' = {pentagon cell, Q}, J = {Q}, I = {i}, and (i) giving
(pentagon cap Q)^i = pentagon^i, condition (ii) says

        pentagon^i  is a FACE of  Q^i.

Q is Minkowski-rigid (slice_rigidity.py), so Q^i = lambda_i Q + v_i and its
faces are lambda_i (faces of Q) + v_i.  Hence the induced decomposition of
the pentagon is

    pentagon = lambda_0 F_0 + lambda_1 F_1 + c,   F_i faces of Q,
               lambda_0 + lambda_1 = 1,  lambda_i >= 0.

A Minkowski sum lies in the affine hull of its summands, so both F_i have
direction space inside the pentagon's plane W.  This script enumerates every
such face -- there are more than the pentagon's own -- and solves exactly,
by support functions, for every (F_0, F_1, lambda_0) making the sum a
translate of the pentagon.

Run:  sage def41_check.sage     (from paper5/)
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

def dirspace(pts):
    p0 = pts[0]
    return matrix(QQ, [list(p - p0) for p in pts[1:]]).row_space()

W = dirspace([vector(QQ, V[j]) for j in PENT])
BW = matrix(QQ, W.basis())                      # 2 x 4, a basis of W

def coords2(p, base):
    """coordinates of p - base in the basis BW (p - base must lie in W)."""
    return vector(QQ, BW.transpose().solve_right(vector(QQ, p) - vector(QQ, base)))

def support(pts2, u):
    return max(p.dot_product(u) for p in pts2)

def centred(P):
    vs = [vector(QQ, v) for v in P.vertices_list()]
    c = sum(vs) / len(vs)
    return sorted(tuple(v - c) for v in vs)

for fi, others in [(21, [3, 9, 10, 13]), (24, [4, 8, 11, 12])]:
    print(f"\n== cell over facet {fi} ==")
    for s in [QQ(1), QQ((1,2)), QQ((1,3)), QQ((1,5)), QQ((1,17)), QQ((2,3)),
              QQ((3,7)), QQ((1,1000))]:
        pts = {i: s * vector(QQ, V[i]) for i in others}
        pts.update({j: vector(QQ, V[j]) for j in PENT})
        Q = Polyhedron(vertices=[list(p) for p in pts.values()], base_ring=QQ)
        inv = {tuple(p): i for i, p in pts.items()}
        # every face of Q whose direction space lies in W, taken up to
        # translation (a Minkowski summand is only defined up to translation)
        cand = []
        for d in (0, 1, 2):
            for f in Q.faces(d):
                fv = [vector(QQ, x) for x in f.vertices()]
                if d == 0 or dirspace(fv).is_subspace(W):
                    lab = tuple(sorted(inv[tuple(x)] for x in fv))
                    cand.append((d, lab, [coords2(x, fv[0]) for x in fv]))
        pent2 = [coords2(vector(QQ, V[j]), V[PENT[0]]) for j in PENT]
        # directions to test: all edge normals of the pentagon and of every
        # candidate face, plus their negatives -- a spanning test set
        dirs = set()
        for _, _, P2 in cand + [(2, "pent", pent2)]:
            for a, b in [(P2[i], P2[j]) for i in range(len(P2))
                         for j in range(len(P2)) if i != j]:
                e = b - a
                if e:
                    dirs.add((e[1], -e[0])); dirs.add((-e[1], e[0]))
        dirs = [vector(QQ, u) for u in dirs]
        hits = []
        for (d0, l0, A) in cand:
            for (d1, l1, B) in cand:
                # unknowns: lam, c0, c1 with lam*hA + (1-lam)*hB = hP + <c,.>
                rows, rhs = [], []
                for u in dirs:
                    rows.append([support(A, u) - support(B, u), u[0], u[1]])
                    rhs.append(support(pent2, u) - support(B, u))
                M = matrix(QQ, rows); b = vector(QQ, rhs)
                try:
                    x = M.solve_right(b)
                except ValueError:
                    continue
                if M * x != b:
                    continue
                lam = x[0]
                if not (0 <= lam <= 1):
                    continue
                # exact confirmation: an actual Minkowski sum, compared to
                # the pentagon after centring (translation is free)
                SA = Polyhedron(vertices=[list(lam * p) for p in A], base_ring=QQ)
                SB = Polyhedron(vertices=[list((1 - lam) * p) for p in B], base_ring=QQ)
                PP = Polyhedron(vertices=[list(p) for p in pent2], base_ring=QQ)
                if centred(SA + SB) == centred(PP):
                    hits.append((l0, l1, lam))
        # A solution is TRIVIAL when each summand is a (possibly degenerate)
        # homothet of the pentagon: lambda_i = 0, or F_i is the pentagon.
        def triv(l0, l1, lam):
            return ((lam == 0 or l0 == tuple(PENT)) and
                    (lam == 1 or l1 == tuple(PENT)))
        nontrivial = [h for h in hits if not triv(*h)]
        ok(f"s = {str(s):<7} solutions of pentagon = lam*F0 + (1-lam)*F1 with "
           f"F_i faces of the cell: {len(hits)} in total, "
           f"{len(nontrivial)} not into homothets of the pentagon"
           + (f"  ->  {nontrivial[:4]}" if nontrivial else ""),
           len(nontrivial) == 0 and len(hits) > 0)

print(f"""
{CH[0]} checks passed.

The only decompositions Definition 4.1 permits at the pentagon are
pentagon^i = lambda_i * pentagon (translated), lambda_0 + lambda_1 = 1.  At
least one lambda_i is positive, and the cone over a positive multiple of the
pentagon IS the cone over the pentagon, so by Ilten-Vollmert Corollary 2.12
-- the general fibre of the deformation has exactly the analytic
singularities of the cones over the summands -- the general fibre still
carries a dP7 cone point.  The dP7 germ is not smoothed, at any degree.""")
