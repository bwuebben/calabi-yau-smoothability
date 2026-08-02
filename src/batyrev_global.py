#!/usr/bin/env python3
"""
Globalization (Batyrev): compact CY threefolds with prescribed non-smoothable
canonical toric singularities.

Setup and facts used (all standard):
  * Delta in Z^4 reflexive => the generic anticanonical hypersurface X in the
    Gorenstein Fano toric 4-fold P_Delta (face fan over Delta) is a Calabi-Yau
    THREEFOLD with Gorenstein canonical singularities [Batyrev, J. Alg. Geom. 3
    (1994)].  Generic X is Delta-regular, so its singularities are exactly the
    ones inherited from the ambient toric strata.
  * A 2-face F of Delta spans the 3-dim cone sigma_F; V(sigma_F) is a toric
    curve (P^1).  The transverse germ of P_Delta along it is the 3-fold toric
    singularity C(F) = cone over the lattice polygon F, taken in the lattice
    induced on aff(F); it is Gorenstein because the two facets of Delta
    containing F have integral supporting functionals equal to 1 there.
  * The restriction of -K to V(sigma_F) is the toric line bundle whose moment
    polytope is the DUAL EDGE F* of the dual polytope Delta*; hence
        #(X  intersect  V(sigma_F)) = lattice length l(F*)  >= 1,
    and at each such point X is analytically the singularity C(F).
  * Edges e of Delta of lattice length m >= 2 similarly give X curves of
    transverse A_{m-1} singularities; edges of length 1 give nothing.  So if
    ALL edges of Delta are unit, Sing X = finitely many points, one C(F)-type
    packet per non-smooth 2-face F.
  * LOCAL-TO-GLOBAL: if some 2-face F is NON-SMOOTHABLE as an isolated
    Gorenstein toric 3-fold singularity — i.e. F has unit edges and admits NO
    Minkowski decomposition into unit segments and standard triangles
    (Altmann Invent. Math. 128 (1997) + Corti-Filip-Petracci arXiv:2006.16885
    Sec 5.1; our census trichotomy: rigid or def-only) — then X admits no
    global smoothing: a smoothing of X would induce at a C(F)-point a local
    deformation with smooth generic fibre, classified by a map to the versal
    base of the isolated germ, which has no smoothing component.

So the hunt is combinatorial: find reflexive Delta with a rigid (or def-only)
2-face.  Anchors: the quintic simplex (all 2-faces smooth => X smooth) and
P(1,1,1,3,3) (2-face = 1/3(1,1,1) rigid triangle, dual edge length 3 => the
classical 3 x 1/3(1,1,1) points on X_9 — a KNOWN non-smoothable CY).  The new
targets are polytopes whose bad 2-face is a rigid NON-simplicial polygon
(non-quotient singularity) or a def-only polygon.
"""
from fractions import Fraction
from itertools import combinations
from math import atan2, gcd
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from toric_census import rigid, smoothing_components   # noqa: E402

# ---------------- exact linear algebra helpers ----------------
def dot(u, v):
    return sum(a * b for a, b in zip(u, v))

def vsub(u, v):
    return tuple(a - b for a, b in zip(u, v))

def vgcd(xs):
    g = 0
    for x in xs:
        g = gcd(g, abs(x))
    return g

def det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))

def hyperplane_normal(p0, p1, p2, p3):
    """Integer normal of the affine span of 4 points in Z^4 (zero if degenerate)."""
    d = [vsub(p, p0) for p in (p1, p2, p3)]
    u = []
    for i in range(4):
        cols = [c for c in range(4) if c != i]
        m = [[d[r][c] for c in cols] for r in range(3)]
        u.append(((-1) ** i) * det3(m))
    return tuple(u)

def affine_rank(points):
    """Rank over Q of {p - points[0]}."""
    if len(points) < 2:
        return 0
    rows = [[Fraction(x) for x in vsub(p, points[0])] for p in points[1:]]
    rank, ncols = 0, len(rows[0])
    for c in range(ncols):
        piv = next((r for r in range(rank, len(rows)) if rows[r][c] != 0), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for r in range(len(rows)):
            if r != rank and rows[r][c] != 0:
                f = rows[r][c] / rows[rank][c]
                rows[r] = [a - f * b for a, b in zip(rows[r], rows[rank])]
        rank += 1
    return rank

def int_kernel(A):
    """Basis of the (saturated) lattice {x in Z^n : A x = 0}, A integer m x n.
    Column-HNF with unimodular column ops tracked in T."""
    m, n = len(A), len(A[0])
    M = [list(row) for row in A]
    T = [[1 if i == j else 0 for j in range(n)] for i in range(n)]

    def colop(j, k, a, b, c, d):        # (col_j, col_k) <- (a*col_j + b*col_k, c*col_j + d*col_k)
        for row in M + T:
            row[j], row[k] = a * row[j] + b * row[k], c * row[j] + d * row[k]

    pivot_col = 0
    for r in range(m):
        # find a column >= pivot_col with M[r][c] != 0, gcd-reduce into pivot_col
        nz = [c for c in range(pivot_col, n) if M[r][c] != 0]
        if not nz:
            continue
        c0 = nz[0]
        if c0 != pivot_col:
            colop(pivot_col, c0, 0, 1, 1, 0)            # swap
        for c in range(pivot_col + 1, n):
            while M[r][c] != 0:
                q = M[r][pivot_col] // M[r][c]
                colop(pivot_col, c, 1, -q, 0, 1)        # col_p -= q*col_c ... then swap roles
                colop(pivot_col, c, 0, 1, 1, 0)         # swap so smaller lands in col_c
            # after loop col_c entry is 0
        if M[r][pivot_col] == 0:                        # possible if all became zero
            continue
        pivot_col += 1
    ker = [tuple(T[i][c] for i in range(n)) for c in range(pivot_col, n)]
    for k in ker:                                       # safety
        assert all(dot(row, k) == 0 for row in A)
    return ker

def solve_int_coords(basis, target):
    """Write integer vector target as integer combination of basis (list of Z^4
    vectors, linearly independent).  Exact rational solve + integrality check."""
    rows = [[Fraction(b[i]) for b in basis] + [Fraction(target[i])] for i in range(len(target))]
    ncols = len(basis)
    rank = 0
    for c in range(ncols):
        piv = next((r for r in range(rank, len(rows)) if rows[r][c] != 0), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        pr = rows[rank]
        rows[rank] = [x / pr[c] for x in pr]
        for r in range(len(rows)):
            if r != rank and rows[r][c] != 0:
                f = rows[r][c]
                rows[r] = [a - f * b for a, b in zip(rows[r], rows[rank])]
        rank += 1
    assert rank == ncols, "basis not independent"
    for r in range(rank, len(rows)):
        assert rows[r][-1] == 0, "target not in span"
    sol = [rows[i][-1] for i in range(ncols)]
    assert all(x.denominator == 1 for x in sol), "non-integral coordinates (lattice not saturated?)"
    return tuple(int(x) for x in sol)

# ---------------- polytope machinery ----------------
def facets(V):
    """All facets of conv(V) (V in Z^4, full-dim).  Returns list of
    (u primitive int normal, c Fraction, frozenset vertex-indices) with
    <u, x> <= c on V and equality exactly on the facet."""
    assert affine_rank(V) == 4, "polytope not full-dimensional"
    out = {}
    for S in combinations(range(len(V)), 4):
        u = hyperplane_normal(*[V[i] for i in S])
        if u == (0, 0, 0, 0):
            continue
        c = dot(u, V[S[0]])
        vals = [dot(u, v) for v in V]
        if all(x <= c for x in vals):
            pass
        elif all(x >= c for x in vals):
            u = tuple(-x for x in u); c = -c; vals = [-x for x in vals]
        else:
            continue
        idx = frozenset(i for i, x in enumerate(vals) if x == c)
        g = vgcd(u)
        key = (tuple(x // g for x in u), Fraction(c, g))
        out[key] = idx
    return [(u, c, idx) for (u, c), idx in out.items()]

def is_reflexive(facs):
    """Reflexive <=> every facet supporting functional is <u,x> = 1 with u in M
    primitive (0 automatically interior)."""
    return all(c == 1 for _, c, _ in facs)

def origin_interior(facs):
    return all(c > 0 for _, c, _ in facs)

def two_faces(V, facs):
    """2-faces as (vertex-index-frozenset, [the 2 facets containing it])."""
    found = {}
    for (f1, f2) in combinations(range(len(facs)), 2):
        I = facs[f1][2] & facs[f2][2]
        if len(I) >= 3 and affine_rank([V[i] for i in I]) == 2:
            found.setdefault(I, set()).update((f1, f2))
    return [(I, sorted(fs)) for I, fs in found.items()]

def face_lattice_polygon(V, face_idx, u1, u2):
    """2-face -> lattice polygon in Z^2 via the induced (saturated) lattice on
    aff(F): basis = int_kernel([u1, u2]); vertices ordered convex-CCW.
    Returns (verts2d, edge_vectors, edge_lengths)."""
    pts = [V[i] for i in face_idx]
    ker = int_kernel([list(u1), list(u2)])
    assert len(ker) == 2
    p0 = pts[0]
    coords = [solve_int_coords(ker, vsub(p, p0)) for p in pts]
    cx = Fraction(sum(c[0] for c in coords), len(coords))
    cy = Fraction(sum(c[1] for c in coords), len(coords))
    coords.sort(key=lambda c: atan2(c[1] - cy, c[0] - cx))
    k = len(coords)
    evs, lens = [], []
    for i in range(k):
        d = vsub(coords[(i + 1) % k], coords[i])
        g = vgcd(d)
        evs.append((d[0] // g, d[1] // g))
        lens.append(g)
    return coords, evs, lens

def classify_polygon(evs, lens):
    """Classify the cone over the 2-face as a 3-fold singularity."""
    k = len(evs)
    A2 = 0
    # 2*area from ordered full edge vectors
    x = y = 0
    for e, l in zip(evs, lens):
        A2 += x * (l * e[1]) - y * (l * e[0])
        x += l * e[0]; y += l * e[1]
    A2 = abs(A2)
    b = sum(lens)
    i = (A2 - b + 2) // 2                       # Pick
    if k == 3 and A2 == 1:
        return dict(status="smooth", k=k, i=i, A2=A2, sc=None)
    if any(l >= 2 for l in lens):
        return dict(status="NON-ISOLATED (A_n edges)", k=k, i=i, A2=A2, sc=None)
    sc = smoothing_components(evs)
    if rigid(evs):
        st = "RIGID (non-smoothable)"
    elif sc == 0:
        st = "def-only (non-smoothable)"
    else:
        st = f"smoothable ({sc} comps)"
    return dict(status=st, k=k, i=i, A2=A2, sc=sc)

def dual_edge_length(u1, u2):
    """Lattice length of the dual edge F* = conv{u1, u2} in Delta*
    = number of C(F)-points on the generic anticanonical X."""
    return vgcd(vsub(u1, u2))

def analyze(name, V, verbose=True):
    """Full report for Delta = conv(V).  Returns dict or None if not reflexive."""
    facs = facets(V)
    if not origin_interior(facs):
        if verbose:
            print(f"{name}: 0 not interior — not a Fano polytope")
        return None
    refl = is_reflexive(facs)
    if verbose:
        print(f"{name}: {len(V)} points, {len(facs)} facets, reflexive: {refl}")
    if not refl:
        return None
    tf = two_faces(V, facs)
    rows, edge_set = [], {}
    for I, fpair in tf:
        assert len(fpair) == 2, "2-face not on exactly 2 facets?"
        u1, u2 = facs[fpair[0]][0], facs[fpair[1]][0]
        coords, evs, lens = face_lattice_polygon(V, I, u1, u2)
        cl = classify_polygon(evs, lens)
        cl["npoints"] = dual_edge_length(u1, u2)
        cl["verts"] = sorted(I)
        cl["edges2d"] = evs
        rows.append(cl)
    # global edges: from all 2-faces, via consecutive vertices in 2d order
    for I, fpair in tf:
        u1, u2 = facs[fpair[0]][0], facs[fpair[1]][0]
        pts = [V[i] for i in I]
        ker = int_kernel([list(u1), list(u2)])
        p0 = pts[0]
        cs = [(solve_int_coords(ker, vsub(p, p0)), p) for p in pts]
        cx = Fraction(sum(c[0][0] for c in cs), len(cs))
        cy = Fraction(sum(c[0][1] for c in cs), len(cs))
        cs.sort(key=lambda c: atan2(c[0][1] - cy, c[0][0] - cx))
        k = len(cs)
        for i in range(k):
            p, q = cs[i][1], cs[(i + 1) % k][1]
            key = tuple(sorted((p, q)))
            edge_set[key] = vgcd(vsub(p, q))
    if verbose:
        for cl in sorted(rows, key=lambda r: (r["status"], r["k"])):
            print(f"   2-face {cl['verts']}: k={cl['k']} i={cl['i']} 2A={cl['A2']}"
                  f"  -> {cl['status']};  X gets {cl['npoints']} such point(s)")
        long_edges = {e: l for e, l in edge_set.items() if l >= 2}
        print(f"   edges: {len(edge_set)} total, {len(long_edges)} of length >= 2"
              + (f"  {sorted(long_edges.values())}" if long_edges else "  (all unit => isolated sings only)"))
    bad = [cl for cl in rows if "non-smoothable" in cl["status"]]
    return dict(reflexive=True, faces=rows, edges=edge_set, bad=bad, facets=facs)


# ================= ANCHORS =================
def anchors():
    print("ANCHOR 1 — quintic simplex (P^4): expect all 2-faces smooth, X smooth")
    e = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
    quintic = e + [(-1, -1, -1, -1)]
    r = analyze("  quintic", quintic)
    assert r and not r["bad"] and all(cl["status"] == "smooth" for cl in r["faces"])
    # NB: Delta edges have length 1 here; X_5 smooth.  (The famous 'degree 5'
    # lives on the dual side.)

    print("\nANCHOR 2 — P(1,1,1,3,3): expect a rigid 1/3(1,1,1) 2-face, 3 points on X_9")
    p11133 = e + [(-1, -1, -3, -3)]
    r = analyze("  P(1,1,1,3,3)", p11133)
    assert r is not None
    rig = [cl for cl in r["bad"] if cl["status"].startswith("RIGID")]
    assert rig, "expected the rigid 1/3(1,1,1) face"
    tri = [cl for cl in rig if cl["k"] == 3 and cl["A2"] == 3 and cl["i"] == 1]
    assert tri and tri[0]["npoints"] == 3, f"expected 3 points of 1/3(1,1,1), got {rig}"
    print("  => X_9 in P(1,1,1,3,3): CY threefold, 3 x 1/3(1,1,1) rigid points"
          " => NON-SMOOTHABLE (recovers the classical example)")
    return True


# ================= SEARCH: rigid quadrilateral 2-face =================
# Target polygon: the smallest rigid non-triangle, edges
# [(-2,-1),(0,-1),(1,0),(1,2)], vertices {(0,0),(-2,-1),(-2,-2),(-1,-2)}.
# Embed in the plane {x3=1, x4=0} (a saturated affine sublattice); translations
# within the plane and the choice w1=(0,0,0,1) are unimodular gauge (see log).
def search_2face(Q2verts, name, B=3, QMAX=3, verbose_hits=2):
    F = [(x, y, 1, 0) for (x, y) in Q2verts]
    w1 = (0, 0, 0, 1)
    hits = []
    tried = 0
    for b1 in range(-B, B + 1):
        for b2 in range(-B, B + 1):
            for b3 in range(-B, 0):                     # some vertex needs x3 < 0
                for q in range(1, QMAX + 1):
                    w2 = (b1, b2, b3, -q)
                    V = F + [w1, w2]
                    if affine_rank(V) != 4:
                        continue
                    tried += 1
                    facs = facets(V)
                    if not origin_interior(facs) or not is_reflexive(facs):
                        continue
                    # F must be a 2-face (vertex set exactly the embedded Q)
                    tf = two_faces(V, facs)
                    fidx = frozenset(range(len(F)))
                    hit = next(((I, fp) for I, fp in tf if I == fidx), None)
                    if hit is None:
                        continue
                    hits.append((V, w2))
                    if len(hits) <= verbose_hits:
                        print(f"\nHIT for {name}: extra vertices w1={w1}, w2={w2}")
                        analyze("  Delta", V)
    print(f"\nsearch[{name}]: {tried} candidates in box, {len(hits)} reflexive hits "
          f"with the target 2-face")
    return hits


# ================= HEADLINE EXAMPLES (found 2026-07-10, asserted) =================
def headline_examples():
    """The two prize polytopes.  Both: reflexive, ALL edges unit (so X has only
    isolated singular points), and every 2-face other than the planted one is
    smooth.  Hence Sing X = exactly 3 points of the planted type, and X is a
    compact CY threefold with NO global smoothing."""
    print("HEADLINE A — Delta_A: only singularities of X = 3 x rigid quadrilateral")
    F = [(x, y, 1, 0) for (x, y) in [(0, 0), (-2, -1), (-2, -2), (-1, -2)]]
    VA = F + [(0, 0, 0, 1), (1, 1, -1, -1)]
    rA = analyze("  Delta_A", VA)
    assert rA is not None
    assert all(l == 1 for l in rA["edges"].values())
    assert len(rA["bad"]) == 1
    bA = rA["bad"][0]
    assert bA["status"].startswith("RIGID") and bA["k"] == 4 and bA["i"] == 1 \
        and bA["npoints"] == 3
    assert all(cl["status"] == "smooth" for cl in rA["faces"] if cl is not bA)
    print("  => X_A: compact CY threefold, Sing = 3 rigid NON-QUOTIENT points C(Q),"
          " NOT smoothable\n")

    print("HEADLINE B — Delta_B: only singularities of X = 3 x def-only pentagon")
    P = [(x, y, 1, 0) for (x, y) in [(0, 0), (-2, -1), (-3, -2), (-2, -3), (-1, -2)]]
    VB = P + [(0, 0, 0, 1), (1, 1, -1, -1)]
    rB = analyze("  Delta_B", VB)
    assert rB is not None
    assert all(l == 1 for l in rB["edges"].values())
    assert len(rB["bad"]) == 1
    bB = rB["bad"][0]
    assert bB["status"].startswith("def-only") and bB["k"] == 5 and bB["i"] == 2 \
        and bB["npoints"] == 3
    assert all(cl["status"] == "smooth" for cl in rB["faces"] if cl is not bB)
    print("  => X_B: compact CY threefold whose singular points all admit local"
          " deformations,\n     yet X_B admits NO smoothing (def-only, globalized)\n")

    print("HEADLINE C — Delta_C: a SINGLE non-smoothable point (found day 3, "
          "src/plant_search.py)")
    # Third completing vertex at x3 = 1 forces the pencil facets (0,0,1,0) and
    # (0,0,1,1), so the dual edge of the planted face has lattice length 1:
    # exactly ONE C(F)-point on X.  (In the two-vertex family l(F*) >= 2.)
    F = [(x, y, 1, 0) for (x, y) in [(0, 0), (-2, -1), (-2, -2), (-1, -2)]]
    VC = F + [(0, 0, 0, 1), (1, 1, -1, 0), (-2, -2, 1, -1)]
    rC = analyze("  Delta_C", VC)
    assert rC is not None
    assert all(l == 1 for l in rC["edges"].values())
    assert len(rC["bad"]) == 1
    bC = rC["bad"][0]
    assert bC["status"].startswith("RIGID") and bC["k"] == 4 and bC["i"] == 1 \
        and bC["npoints"] == 1
    assert all(cl["status"] == "smooth" for cl in rC["faces"] if cl is not bC)
    assert len(rC["facets"]) == 10 and len(rC["faces"]) == 21
    print("  => X_C: compact CY threefold, Sing = ONE anticanonical-cone-over-F1"
          " point,\n     NOT smoothable — the single-point example (paper Sec. 5.3)\n")

    print("HEADLINE C' — same completions on the 1/3(1,1,1) triangle: single"
          " quotient point")
    T = [(x, y, 1, 0) for (x, y) in [(0, 0), (-2, -1), (-1, -2)]]
    VP = T + [(0, 0, 0, 1), (1, 1, -1, 0), (-2, -2, 1, -1)]
    rP = analyze("  Delta_C'", VP, verbose=False)
    assert rP is not None
    assert all(l == 1 for l in rP["edges"].values())
    assert len(rP["bad"]) == 1 and rP["bad"][0]["npoints"] == 1 \
        and rP["bad"][0]["status"].startswith("RIGID") and rP["bad"][0]["k"] == 3
    print("  => X: compact CY threefold with a SINGLE 1/3(1,1,1) point"
          " (X_9 classically has three)\n")

    print("HEADLINE D — Delta_D (found by the KS scan, src/ks_sweep.py, v10):"
          " a SINGLE DEF-ONLY point")
    VD = [(1, 0, 0, 0), (0, 1, 0, 0), (0, -1, 0, 0), (0, 0, 1, 0),
          (-3, 1, -1, 0), (0, 1, 1, 0), (0, 0, 0, 1), (-3, 2, 1, -1),
          (-2, 2, 2, -1), (1, 0, 1, 1)]
    rD = analyze("  Delta_D", VD, verbose=False)
    assert rD is not None
    assert all(l == 1 for l in rD["edges"].values()) and len(rD["edges"]) == 30
    assert len(rD["facets"]) == 15 and len(rD["faces"]) == 35
    assert len(rD["bad"]) == 1
    bD = rD["bad"][0]
    assert bD["status"].startswith("def-only") and bD["k"] == 5 and bD["i"] == 2 \
        and bD["A2"] == 7 and bD["npoints"] == 1
    assert all(cl["status"] == "smooth" for cl in rD["faces"] if cl is not bD)
    # the unique bad face is GL2(Z)-equivalent to the Theorem B pentagon T+s:
    from toric_census import equiv, ccw_sort
    QB = ccw_sort([(-2, -1), (-1, -1), (1, -1), (1, 1), (1, 2)])
    assert equiv(ccw_sort(bD["edges2d"]), QB)
    print("  => X_D: compact CY threefold, Sing = ONE point of the def-only"
          " pentagon type C(Q_B):\n     the point deforms non-trivially, yet"
          " X_D admits NO smoothing (paper Thm D)\n")


if __name__ == "__main__":
    anchors()
    print("\n" + "=" * 72)
    headline_examples()
    print("=" * 72)
    print("SEARCH — rigid quadrilateral [(-2,-1),(0,-1),(1,0),(1,2)] as a 2-face")
    quad = [(0, 0), (-2, -1), (-2, -2), (-1, -2)]
    hits_quad = search_2face(quad, "rigid-quad", verbose_hits=0)
    print("\n" + "=" * 72)
    print("SEARCH — rigid 1/3(1,1,1) triangle as a 2-face (control: should be findable)")
    tri = [(0, 0), (-2, -1), (-1, -2)]
    hits_tri = search_2face(tri, "P2-triangle", verbose_hits=0)
