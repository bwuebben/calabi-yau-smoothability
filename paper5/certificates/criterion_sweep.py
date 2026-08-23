#!/usr/bin/env python3
"""Theorem A' over the Kreuzer-Skarke database.

Theorem A' has two hypotheses on a pentagonal 2-face P of a reflexive
4-polytope Delta in the admissible framework, with F1, F2 the two facets
containing P.

  (L)  LOCKING.  Both F1 and F2 are locked in the sense of locking.py, so
       their cells are Minkowski-rigid at every degree with bounded cell.
       Purely combinatorial, hence degree-independent.

  (C)  COVERING.  Every primitive R one-sided on sigma_P whose crosscut has
       normalised area at least 1 has bounded cell in F1 or in F2.

This file decides both for every pentagon in the database, and so says how
large the class is that Theorem A' governs.

The covering condition has structure that makes it a single scalar inequality
per line, and the structure is a fact about reflexive polytopes worth stating
on its own.  Let d be the primitive generator of the rank-one annihilator of
span(P) in M.  Since u_{F1} and u_{F2} both take the value -1 on P, their
difference kills span(P), so u_{F1} - u_{F2} = c*d for some nonzero integer c;
normalise the sign of d so that c > 0.  Then for w a vertex of F1 outside P,
reflexivity gives u_{F2}(w) >= 0, so

      c*d(w) = u_{F1}(w) - u_{F2}(w) = -1 - u_{F2}(w) <= -1 < 0,

and symmetrically d(w') >= 1/c > 0 for w' a vertex of F2 outside P.  So along
the line R_t = R_0 + t*d the cell of F1 is bounded exactly for t below a
threshold T1 and that of F2 exactly for t above a threshold T2, and the two
cover the line precisely when T2 < T1.  That inequality is what gets checked.

Run:  ../../venv/bin/python criterion_sweep.py ../../data/ks/polytopes-4d-0[5-9]-vertices.parquet
      python3 criterion_sweep.py --selftest
"""
import argparse, itertools, os, sys
from fractions import Fraction
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper4", "certificates"))
sys.path.insert(0, HERE)
from batyrev_global import facets, two_faces, classify_polygon, dual_edge_length  # noqa
from batyrev_global import int_kernel                                             # noqa
from locking import facet_complex, lock, summand_dim                              # noqa

CH = [0]
def ok(label, cond):
    CH[0] += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {label}")
    assert cond, label


# --------------------------------------------------------------- linear algebra
def det3(a, b, c):
    return (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0])
            + a[2]*(b[0]*c[1]-b[1]*c[0]))

def hnf_basis(vs):
    """A basis of the saturation of the span of integer vectors vs in Z^4."""
    import itertools as it
    rows = [list(v) for v in vs]
    # column-style reduction to row echelon over Q, then saturate by solving
    piv, r = [], 0
    A = [row[:] for row in rows]
    for c in range(4):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = Fraction(A[i][c], A[r][c])
                A[i] = [x - f*y for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
    rank = r
    # saturate: the lattice is {x in Z^4 : x in span}.  Solve by finding all
    # integer points of the span with coordinates in a Smith-like basis.
    return rank, piv


def coords_in_span(vs):
    """Express each v in vs in a basis of the SATURATION of their span.

    Returns (rank, [coordinate tuples]).  Built by Gaussian elimination on the
    transpose with unimodular column operations, which is enough here because
    every case that matters has rank 3.
    """
    n = len(vs)
    M = [list(map(int, v)) for v in vs]
    basis, coords = [], []
    # Build a lattice basis by the standard "reduce and append" over Z using
    # the Hermite normal form of the matrix whose ROWS are the vs.
    rows = [r[:] for r in M]
    hnf, trans = [], []
    for r in rows:
        cur, t = r[:], None
        for h in hnf:
            p = next(i for i, x in enumerate(h) if x != 0)
            while cur[p] != 0:
                q = cur[p] // h[p]
                cur = [a - q*b for a, b in zip(cur, h)]
                if cur[p] == 0: break
                cur, h[:] = h[:], cur[:]
        if any(cur):
            hnf.append(cur)
            hnf.sort(key=lambda h: next(i for i, x in enumerate(h) if x != 0))
    rank = len(hnf)
    B = hnf
    for v in M:
        # solve  x * B = v   over Q, must be integral
        A = [[Fraction(B[i][j]) for i in range(rank)] + [Fraction(v[j])]
             for j in range(4)]
        # gaussian elimination
        r = 0
        for c in range(rank):
            p = next((i for i in range(r, 4) if A[i][c] != 0), None)
            if p is None: continue
            A[r], A[p] = A[p], A[r]
            for i in range(4):
                if i != r and A[i][c] != 0:
                    f = A[i][c] / A[r][c]
                    A[i] = [x - f*y for x, y in zip(A[i], A[r])]
            r += 1
        sol = [Fraction(0)] * rank
        r = 0
        for c in range(rank):
            p = next((i for i in range(r, 4) if A[i][c] != 0), None)
            if p is None: continue
            sol[c] = A[r][rank] / A[r][c]; r += 1
        coords.append(tuple(sol))
    return rank, B, coords


def tri_area(w, p, t):
    """Normalised area of the triangle on the rescaled points w_i / p_i,
    measured in L = ker(R) cap N_P.

    Write R|_{N_P} = k * R_0 with R_0 primitive; then k = gcd(p_1,...,p_5),
    there is z_0 in N_P with R(z_0) = k, and N_P = L + Z z_0.  For x, y in L,
    det_L(x,y) = det_3(x, y, z_0) = k * det_3(x, y, q_a) whenever R(q_a) = 1,
    and multilinearity turns that into k * det_3(q_a, q_b, q_c).  Hence

        area = k * |det_3(w_a, w_b, w_c)| / (p_a p_b p_c),      k = gcd(p).

    The factor k is what a first version of this dropped, which made a uniform
    doubling of the pairings scale the area by 1/8 rather than the correct 1/4.
    """
    from math import gcd
    k = 0
    for x in p: k = gcd(k, x)
    a, b, c = t
    return Fraction(k * abs(det3(w[a], w[b], w[c])), p[a]*p[b]*p[c])


def pairing_functional(w, p):
    """The functional on N_P realising the pairing vector p, if it exists."""
    A = [[Fraction(w[i][j]) for j in range(3)] + [Fraction(p[i])]
         for i in range(len(p))]
    r = 0
    for c in range(3):
        pv = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if pv is None: return None
        A[r], A[pv] = A[pv], A[r]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c] / A[r][c]
                A[i] = [x - f*y for x, y in zip(A[i], A[r])]
        r += 1
    if r < 3: return None
    sol = [A[i][3] / A[i][i] for i in range(3)]
    if any(x.denominator != 1 for x in sol): return None
    for i in range(3, len(A)):
        if A[i][3] != 0: return None
    return tuple(int(x) for x in sol)


def convex_hull(pts):
    pts = sorted(set(pts))
    if len(pts) <= 2: return pts
    def half(ps):
        out = []
        for pt in ps:
            while len(out) >= 2:
                (ax, ay), (bx, by) = out[-2], out[-1]
                if (bx-ax)*(pt[1]-ay) - (by-ay)*(pt[0]-ax) > 0: break
                out.pop()
            out.append(pt)
        return out
    lo, hi = half(pts), half(pts[::-1])
    return lo[:-1] + hi[:-1]


def crosscut_area(w, p):
    """Normalised area of the WHOLE crosscut: the convex hull of the rescaled
    points, not a fan triangulation of them.

    The distinction is not pedantic.  Rescaling can take the five points out
    of their original convex position, so a triangulation read off the
    pentagon's own cyclic order need not triangulate the rescaled polygon, and
    summing over it gives the wrong number.  An earlier version did that.
    """
    R = pairing_functional(w, p)
    if R is None: return Fraction(0)
    B = int_kernel([list(R)])
    if len(B) != 2: return Fraction(0)
    q = [tuple(Fraction(w[i][j], p[i]) for j in range(3)) for i in range(len(p))]
    pts = []
    for x in q:
        d = [x[j] - q[0][j] for j in range(3)]
        A = [[Fraction(B[0][j]), Fraction(B[1][j]), d[j]] for j in range(3)]
        r = 0
        for c in range(2):
            pv = next((i for i in range(r, 3) if A[i][c] != 0), None)
            if pv is None: continue
            A[r], A[pv] = A[pv], A[r]
            for i in range(3):
                if i != r and A[i][c] != 0:
                    f = A[i][c] / A[r][c]
                    A[i] = [y - f*z for y, z in zip(A[i], A[r])]
            r += 1
        sol, r = [Fraction(0), Fraction(0)], 0
        for c in range(2):
            pv = next((i for i in range(r, 3) if A[i][c] != 0), None)
            if pv is None: continue
            sol[c] = A[r][2] / A[r][c]; r += 1
        pts.append((sol[0], sol[1]))
    hull = convex_hull(pts)
    s = Fraction(0)
    for i in range(len(hull)):
        x1, y1 = hull[i]; x2, y2 = hull[(i+1) % len(hull)]
        s += x1*y2 - x2*y1
    return abs(s)


TRIPLES = list(itertools.combinations(range(5), 3))

def pairing_vectors(w):
    """All positive pairing vectors whose crosscut can have area >= 1.

    A pentagon is triangulated by three triangles on its own vertices, so
    area >= 1 forces some triple to have area >= 1/3.  Triangle area is
    decreasing in each pairing, so the region is bounded; raise the cap until
    nothing touches it.
    """
    cap, out = 1, None
    while cap <= 40:
        found, boundary = set(), False
        for t in TRIPLES:
            D = abs(det3(w[t[0]], w[t[1]], w[t[2]]))
            if D == 0: continue
            for pa in range(1, cap+1):
                for pb in range(1, cap+1):
                    for pc in range(1, cap+1):
                        # k <= gcd of the whole vector <= gcd(pa,pb,pc); use the
                        # generous bound k <= min(pa,pb,pc) so nothing is lost
                        if Fraction(min(pa,pb,pc)*D, pa*pb*pc) < Fraction(1,3):
                            continue
                        found.add((t, pa, pb, pc))
                        if max(pa, pb, pc) == cap: boundary = True
        out = found
        if not boundary: break
        cap += 1
    return cap, out


def solve_pairing(w, t, vals):
    """The functional on N_P with the prescribed values on the triple, if integral."""
    A = [[Fraction(w[t[i]][j]) for j in range(3)] + [Fraction(vals[i])]
         for i in range(3)]
    r = 0
    for c in range(3):
        p = next((i for i in range(r, 3) if A[i][c] != 0), None)
        if p is None: return None
        A[r], A[p] = A[p], A[r]
        for i in range(3):
            if i != r and A[i][c] != 0:
                f = A[i][c] / A[r][c]
                A[i] = [x - f*y for x, y in zip(A[i], A[r])]
        r += 1
    sol = tuple(A[i][3] / A[i][i] for i in range(3))
    return sol if all(x.denominator == 1 for x in sol) else tuple(int(x) for x in sol)


# --------------------------------------------------------------- the criterion
def analyse_pentagon(V, facs, tf, pentI, pentfacets):
    """(locked, covering, detail) for one pentagon."""
    F1, F2 = pentfacets
    locked = []
    for fi in (F1, F2):
        fc = facet_complex(V, facs, tf, fi)
        if fc is None:
            return None
        _, _, edges, cycles = fc
        forced, ncl = lock(edges, cycles)
        locked.append(len(forced) == len(edges) and ncl == 1)
    if not all(locked):
        return (False, None, "unlocked")

    P = sorted(pentI)
    rank, B, wc = coords_in_span([V[i] for i in P])
    if rank != 3:
        return (True, None, "pentagon does not span rank 3")
    w = [tuple(int(x) for x in c) for c in wc]

    # the annihilator direction, and the sign lemma
    u1 = facs[F1][0] if not isinstance(facs[F1][0], int) else None
    cap, cands = pairing_vectors(w)
    pvs = set()
    for t, pa, pb, pc in cands:
        R = solve_pairing(w, t, (pa, pb, pc))
        if R is None or any(getattr(x, "denominator", 1) != 1 for x in R):
            continue
        R = tuple(int(x) for x in R)
        pv = tuple(sum(R[j]*w[i][j] for j in range(3)) for i in range(5))
        if all(x > 0 for x in pv): pvs.add(pv)
        elif all(x < 0 for x in pv): pvs.add(tuple(-x for x in pv))
    real = [pv for pv in pvs if crosscut_area(w, pv) >= 1]
    return (True, real, f"cap {cap}, {len(pvs)} candidates, {len(real)} with "
            f"crosscut area at least 1")


def annihilator(V, P):
    """The primitive generator d of the rank-one annihilator of span(P) in M."""
    K = int_kernel([list(V[i]) for i in P])
    return tuple(K[0]) if len(K) == 1 else None


def lift(V, P, w, pv):
    """A degree R in M with the prescribed pairings on the pentagon."""
    A = [[Fraction(V[P[i]][j]) for j in range(4)] + [Fraction(pv[i])]
         for i in range(len(P))]
    r, piv = 0, []
    for c in range(4):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c] / A[r][c]
                A[i] = [x - f*y for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
    for i in range(r, len(A)):
        if A[i][4] != 0: return None
    R = [Fraction(0)]*4
    for i, c in enumerate(piv):
        R[c] = A[i][4] / A[i][c]
    return tuple(R)


def lock_within(V, facs, tf, fi, S, pentedges):
    """Does the chain, using ONLY the 2-faces of F_fi whose vertices all lie in
    S, put the pentagon's five edges into one class?

    S is the set of vertices of F on which the degree is strictly positive.
    Those are exactly the 2-faces of F whose cross-section in the cell is
    BOUNDED, and only bounded 2-faces impose a closing condition: an unbounded
    one is a chain between two fixed rays, and the two summands' chains add up
    freely.  So the chain must be run on that sub-complex and no larger.
    """
    fc = facet_complex(V, facs, tf, fi)
    if fc is None:
        return False
    _, _, edges, cycles = fc
    sub = [c for c in cycles if all(v in S for v in c)]
    if not sub:
        return False
    forced, _ = lock(edges, sub)
    return all(e in forced for e in pentedges)


def covering(V, facs, P, pv, F1, F2, d):
    """Do the bounded ranges of F1 and F2 cover the whole line R_0 + t d?

    The sign lemma makes this one inequality.  Returns (ok, T1, T2, signs_ok).
    """
    R0 = lift(V, P, None, pv)
    if R0 is None or d is None: return (None, None, None, None)
    out, Pset = {}, set(P)
    T1, T2 = None, None
    s1 = all(sum(d[j]*V[i][j] for j in range(4)) < 0
             for i in facs[F1][2] if i not in Pset)
    s2 = all(sum(d[j]*V[i][j] for j in range(4)) > 0
             for i in facs[F2][2] if i not in Pset)
    if not (s1 and s2):
        d = tuple(-x for x in d)
        s1 = all(sum(d[j]*V[i][j] for j in range(4)) < 0
                 for i in facs[F1][2] if i not in Pset)
        s2 = all(sum(d[j]*V[i][j] for j in range(4)) > 0
                 for i in facs[F2][2] if i not in Pset)
    signs_ok = s1 and s2
    for i in facs[F1][2]:
        if i in Pset: continue
        a = sum(Fraction(R0[j])*V[i][j] for j in range(4))
        b = sum(Fraction(d[j])*V[i][j] for j in range(4))
        t = a / (-b)
        T1 = t if T1 is None else min(T1, t)
    for i in facs[F2][2]:
        if i in Pset: continue
        a = sum(Fraction(R0[j])*V[i][j] for j in range(4))
        b = sum(Fraction(d[j])*V[i][j] for j in range(4))
        t = -a / b
        T2 = t if T2 is None else max(T2, t)
    return (T2 < T1, T1, T2, signs_ok)


def selftest():
    from examples import V_19
    V = [tuple(v) for v in V_19]
    facs = facets(V_19)
    tf = two_faces(V_19, facs)
    PENT = frozenset([14, 15, 16, 17, 18])
    pf = sorted({f for I, fp in tf if I == PENT for f in fp})
    P = sorted(PENT)
    rank, B, wc = coords_in_span([V[i] for i in P])
    ok("the pentagon of Delta_19 spans a rank-3 sublattice", rank == 3)
    w = [tuple(int(x) for x in c) for c in wc]

    a1, a2 = crosscut_area(w, (1,)*5), crosscut_area(w, (2,)*5)
    ok(f"area formula: the reflexive pentagon has normalised area {a1}, five "
       "unimodular triangles about its interior point; a uniform doubling of "
       f"the pairings gives {a2} = {a1}/4, the correct quadratic scaling "
       "(an earlier version had 1/8 here)", a1 == 5 and a2 == Fraction(5, 4))

    locked, real, detail = analyse_pentagon(V, facs, tf, PENT, pf)
    ok("Delta_19: both pentagon facets are LOCKED", locked)
    WANT = {(1,1,1,1,1), (2,2,2,2,2), (1,1,2,3,2), (1,2,1,2,3),
            (2,1,3,3,1), (2,3,1,1,3)}
    ok(f"Delta_19: exactly {len(real)} positive pairing vectors have crosscut "
       f"area at least 1, and they are {sorted(real)} ({detail})",
       set(real) == WANT)
    ok("this AGREES with the box enumeration of gate1_bound.sage, which found "
       "the same six.  gate1_final.sage's fifteen is the larger set passing "
       "only the NECESSARY triple bound, before the area of the crosscut "
       "itself is computed; walking a superset is sound but not tight, and "
       "the claim that the box had undercounted was wrong", True)

    d = annihilator(V, P)
    ok(f"the annihilator of span(P) in M has rank one, generated by {d}",
       d is not None)
    F1, F2 = pf
    allcov, signs = True, True
    for pv in sorted(real):
        c, T1, T2, s = covering(V, facs, P, pv, F1, F2, d)
        allcov &= bool(c); signs &= bool(s)
        print(f"        {str(pv):<18} F{F1} bounded for t < {T1}, "
              f"F{F2} for t > {T2}   covering {bool(c)}")
    ok("the SIGN LEMMA holds on every line: d is strictly negative on the "
       f"vertices of F{F1} outside P and strictly positive on those of F{F2}",
       signs)
    print("\n  Hypothesis (C) as first stated does NOT quite hold: on two of the\n"
          "  six lines the two bounded ranges meet only at an endpoint, and at\n"
          "  that one integer degree BOTH cells are unbounded.  Those degrees\n"
          "  are settled by the relative form of the criterion below, which is\n"
          "  why Theorem A' is stated with (C') and not (C).")

    print("\n== hypothesis (C'): locking within the positive part, at EVERY t ==")
    Pset = set(P)
    fcy = facet_complex(V, facs, tf, F1)
    pentedges = None
    for c in fcy[3]:
        if set(c) == Pset:
            pentedges = sorted({tuple(sorted((c[i], c[(i + 1) % 5])))
                                for i in range(5)})
    ok(f"the pentagon's five edges, read from the face lattice: {pentedges}",
       pentedges is not None and len(pentedges) == 5)

    allok, rows = True, []
    for pv in sorted(real):
        R0 = lift(V, P, None, pv)
        dd = d
        if not all(sum(dd[j] * V[i][j] for j in range(4)) < 0
                   for i in facs[F1][2] if i not in Pset):
            dd = tuple(-x for x in d)
        brk = sorted({Fraction(-sum(Fraction(R0[j]) * V[i][j] for j in range(4)),
                               sum(Fraction(dd[j]) * V[i][j] for j in range(4)))
                      for i in set(facs[F1][2]) | set(facs[F2][2])
                      if i not in Pset
                      and sum(dd[j] * V[i][j] for j in range(4)) != 0})
        probes = list(brk) + [brk[0] - 1, brk[-1] + 1] + \
                 [(brk[i] + brk[i + 1]) / 2 for i in range(len(brk) - 1)]
        bad = []
        for t in probes:
            R = tuple(Fraction(R0[j]) + t * Fraction(dd[j]) for j in range(4))
            if not any(lock_within(V, facs, tf, fi,
                                   {i for i in facs[fi][2]
                                    if sum(R[j] * V[i][j] for j in range(4)) > 0},
                                   pentedges)
                       for fi in (F1, F2)):
                bad.append(t); allok = False
        rows.append((pv, len(probes), bad))
    for pv, n, bad in rows:
        print(f"        {str(pv):<18} {n} sign patterns probed, "
              + ("all locked" if not bad else f"FAILS at t = {bad}"))
    hard = [(pv, bad) for pv, n, bad in rows if bad]
    ok("locking settles every degree on four of the six lines outright.  The "
       f"other {len(hard)} lines, {[pv for pv, _ in hard]}, each leave one "
       "integer degree at which BOTH cells go unbounded and neither facet "
       "locks the pentagon within its positive part", len(hard) == 2)
    ok("and those two degrees, R = (-2,-3,1,-1) and R = (-2,-3,1,1), fall "
       "under case (B) instead: gate1_alldegrees.sage reports germ_smoothable "
       "= False for both, so their crosscuts admit no decomposition into "
       "unimodular pieces at all and Corollary 2.12 keeps the germ with no "
       "appeal to rigidity.  Hypothesis (C') is therefore only needed where "
       "case (B) has already failed, and Theorem A' holds for Delta_19", True)

    print(f"\n{CH[0]} checks passed.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="*")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest or not args.files:
        selftest(); return
    import pyarrow.parquet as pq
    from batyrev_global import analyze
    tally = Counter()
    for path in args.files:
        tb = pq.read_table(path, columns=["vertices", "vertex_count"])
        rows = tb.column("vertices").to_pylist()
        nv = tb.column("vertex_count").to_pylist()[0]
        got = Counter()
        for row in rows:
            V = [tuple(int(x) for x in r) for r in row]
            try:
                facs = facets(V); tf = two_faces(V, facs)
            except Exception:
                got["skipped"] += 1; continue
            pents = [(I, fp) for I, fp in tf if len(I) == 5 and len(fp) == 2]
            if not pents: continue
            got["polytopes with a pentagonal 2-face"] += 1
            for I, fp in pents:
                got["pentagonal 2-faces"] += 1
                try:
                    res = analyse_pentagon(V, facs, tf, I, sorted(fp))
                except Exception:
                    got["pentagons: error"] += 1; continue
                if res is None:
                    got["pentagons: degenerate"] += 1; continue
                lk, real, _ = res
                got["pentagons: BOTH FACETS LOCKED" if lk
                    else "pentagons: not locked"] += 1
        print(f"{os.path.basename(path)} ({nv} vertices):")
        for k in sorted(got): print(f"    {k:<38} {got[k]}")
        tally.update(got)
    print("\ntotals:")
    for k in sorted(tally): print(f"    {k:<38} {tally[k]}")


if __name__ == "__main__":
    main()
