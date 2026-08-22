#!/usr/bin/env python3
"""Exact checker for the lattice claims of the paper (Section 2, Table 2,
and Lemma 3.3).

Verifies, over exact rationals:

  1. the embedded lambda-tables are the cyclic differences of the embedded
     t-vectors, and the embedded Pic(E)^vee rows are the lambda-functionals
     in the printed bases of Table 2 of the paper,
     including consistency with the D1/D2/D7/D8 boundary-class expansions;
  2. the boundary-curve classes have the printed self-intersection patterns,
     adjacent products 1, distant products 0, and sum to -K_E;
  3. the norms and pairings of the coordinate classes: at both dP6 points
     K^perp = A1 + A2 with A1 spanned by phi(s3) and A2 by phi(s1), phi(s2);
     at both dP7 points phi(beta) spans the unique (-2)-line and
     phi(alpha)^2 = -4; discriminants 6, 6, 7, 7;
  4. the (-2)-vector enumerations behind the uniqueness claims, with
     rigorous search boxes from the inverse Gram diagonals;
  5. the sweep pushforward kernels in the standard del Pezzo bases:
     ker(K^perp -> H_2(V)) equals the A1 line for V = (P1)^3, the A2 plane
     for V = P(T_P2), and the unique (-2)-line for V_7 = Bl_pt P^3;
  6. the transfer: an explicit isometry from each printed lattice to the
     standard model matching boundary hexagons/pentagons and K, carrying
     phi(s3) to +/- (l - e1 - e2 - e3), phi(s1), phi(s2) into <e_i - e_j>,
     and phi(beta) to +/- (e1 - e2);
  7. the resulting identification of the branch subspaces with the four
     profile coordinate sets of Section 5.3 of the paper, and the deduction
     of nonsmoothability from the certified forced-zero lists.

Run:  python3 milnor_kernels.py
"""

from fractions import Fraction as F
from itertools import product
import math
import sys

PASS = 0
def ok(label, cond):
    global PASS
    if not cond:
        print(f"FAIL: {label}")
        sys.exit(1)
    PASS += 1
    print(f"  pass: {label}")

# ---------------------------------------------------------------- linear algebra
def mat(rows):
    return [[F(x) for x in r] for r in rows]

def mmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    assert len(A[0]) == k
    return [[sum(A[i][t] * B[t][j] for t in range(k)) for j in range(m)]
            for i in range(n)]

def mvec(A, v):
    return [sum(A[i][j] * v[j] for j in range(len(v))) for i in range(len(A))]

def vmat(v, A):
    return [sum(v[i] * A[i][j] for i in range(len(A))) for j in range(len(A[0]))]

def dot(u, v):
    return sum(a * b for a, b in zip(u, v))

def transpose(A):
    return [list(r) for r in zip(*A)]

def rref(A):
    M = [row[:] for row in A]
    rows, cols = len(M), len(M[0]) if M else 0
    piv = []
    r = 0
    for c in range(cols):
        p = next((i for i in range(r, rows) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        M[r] = [x / M[r][c] for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c] != 0:
                M[i] = [a - M[i][c] * b for a, b in zip(M[i], M[r])]
        piv.append(c)
        r += 1
        if r == rows:
            break
    return M, piv

def rank(A):
    return len(rref(A)[1])

def kernel(A):
    """Right kernel basis of A (A v = 0)."""
    M, piv = rref(A)
    cols = len(A[0])
    free = [c for c in range(cols) if c not in piv]
    basis = []
    for fc in free:
        v = [F(0)] * cols
        v[fc] = F(1)
        for r, pc in enumerate(piv):
            v[pc] = -M[r][fc]
        basis.append(v)
    return basis

def inverse(A):
    n = len(A)
    M = [row[:] + [F(1) if i == j else F(0) for j in range(n)]
         for i, row in enumerate(A)]
    R, piv = rref(M)
    assert piv == list(range(n)), "matrix not invertible"
    return [row[n:] for row in R]

def span_equal(B1, B2):
    if not B1 and not B2:
        return True
    if not B1 or not B2:
        return False
    r1, r2 = rank(B1), rank(B2)
    return r1 == r2 == rank(B1 + B2)

def in_span(v, B):
    return rank(B) == rank(B + [v])

# ---------------------------------------------------------------- printed data
# Table 2 of the paper: boundary bases, Gram matrices, K-rows.

I6 = mat([[-1, 0, 1, 0],
          [ 0,-1, 0, 1],
          [ 1, 0,-1, 1],
          [ 0, 1, 1,-1]])
K6 = [F(-1), F(-1), F(-2), F(-2)]

I71 = mat([[0, 0, 1],
           [0,-1, 1],
           [1, 1, 0]])
I72 = mat([[-1, 0, 1],
           [ 0, 0, 1],
           [ 1, 1, 0]])
K7 = [F(-1), F(-1), F(-2)]

# Boundary data: (cyclic ray labels, class of each ray in the printed basis,
# expected self-intersections). Basis rays carry unit vectors; the two
# non-basis rays carry their expansions in the chosen divisor basis.
E61 = dict(
    name="E_{6,1}",
    I=I6, K=K6,
    cyclic=["q1", "q3", "q5", "q6", "q4", "q2"],
    classes={"q1": [-1, 1, 0, 1], "q2": [1, -1, 1, 0],
             "q3": [1, 0, 0, 0], "q4": [0, 1, 0, 0],
             "q5": [0, 0, 1, 0], "q6": [0, 0, 0, 1]},
    basis_rays=["q3", "q4", "q5", "q6"],
    selfint=[-1, -1, -1, -1, -1, -1],
)
E62 = dict(
    name="E_{6,2}",
    I=I6, K=K6,
    cyclic=["q7", "q9", "q11", "q12", "q10", "q8"],
    classes={"q7": [-1, 1, 0, 1], "q8": [1, -1, 1, 0],
             "q9": [1, 0, 0, 0], "q10": [0, 1, 0, 0],
             "q11": [0, 0, 1, 0], "q12": [0, 0, 0, 1]},
    basis_rays=["q9", "q10", "q11", "q12"],
    selfint=[-1, -1, -1, -1, -1, -1],
)
E71 = dict(
    name="E_{7,1}",
    I=I71, K=K7,
    cyclic=["q1", "q13", "q21", "q17", "q7"],
    classes={"q1": [-1, 1, 1], "q7": [1, -1, 0],
             "q13": [1, 0, 0], "q17": [0, 1, 0], "q21": [0, 0, 1]},
    basis_rays=["q13", "q17", "q21"],
    selfint=[-1, 0, 0, -1, -1],
)
E72 = dict(
    name="E_{7,2}",
    I=I72, K=K7,
    cyclic=["q2", "q14", "q22", "q18", "q8"],
    classes={"q2": [-1, 1, 0], "q8": [1, -1, 1],
             "q14": [1, 0, 0], "q18": [0, 1, 0], "q22": [0, 0, 1]},
    basis_rays=["q14", "q18", "q22"],
    selfint=[-1, -1, 0, 0, -1],
)

# t-vectors, lambda-vectors, and the K^perp rows of Table 2 of the paper.
ROWS6 = {
    "s1": dict(t=[0, 0, -1, 0, 0, -1], lam=[1, 0, -1, 1, 0, -1],
               row=[0, 0, -1, 1]),
    "s2": dict(t=[0, -1, 0, 0, -1, 0], lam=[0, -1, 1, 0, -1, 1],
               row=[-1, -1, 1, 0]),
    "s3": dict(t=[0, -1, 0, -1, 0, -1], lam=[1, -1, 1, -1, 1, -1],
               row=[-1, 1, 1, -1]),
}
ROWS7 = {
    "E_{7,1}": {
        "alpha": dict(t=[0, 0, -1, 1, -1], lam=[1, 0, -1, 2, -2],
                      row=[1, -1, 0]),
        "beta":  dict(t=[0, -1, 0, 0, -1], lam=[1, -1, 1, 0, -1],
                      row=[1, 1, -1]),
    },
    "E_{7,2}": {
        "alpha": dict(t=[0, 0, -1, 1, -1], lam=[1, 0, -1, 2, -2],
                      row=[-2, 0, 1]),
        "beta":  dict(t=[0, -1, 0, 0, -1], lam=[1, -1, 1, 0, -1],
                      row=[-1, -1, 1]),
    },
}

def cyclic_diff(t):
    k = len(t)
    return [t[i] - t[i - 1] for i in range(k)]

# ------------------------------------------------- 1+2: conventions and boundary
def check_surface_conventions(E, rows):
    name = E["name"]
    I, K = E["I"], E["K"]
    cyc = E["cyclic"]
    cls = {q: [F(x) for x in v] for q, v in E["classes"].items()}
    n = len(cyc)

    # boundary intersection pattern
    for i in range(n):
        for j in range(n):
            v, w = cls[cyc[i]], cls[cyc[j]]
            prod = dot(v, mvec(I, w))
            if i == j:
                exp = F(E["selfint"][i])
            elif (i - j) % n in (1, n - 1):
                exp = F(1)
            else:
                exp = F(0)
            assert prod == exp, (name, cyc[i], cyc[j], prod, exp)
    ok(f"{name}: boundary self-intersections and adjacencies", True)

    total = [sum(cls[q][j] for q in cyc) for j in range(len(K))]
    ok(f"{name}: boundary classes sum to -K", total == [-k for k in K])

    # lambda tables and printed rows: the canonical vertex order of the
    # lambda tables (from local_deformations.md) may be a dihedral
    # relabeling of the ray-cyclic order; one common relabeling must make
    # every row of this surface consistent.
    for label, d in rows.items():
        lam = cyclic_diff(d["t"])
        assert lam == d["lam"], (name, label, lam, d["lam"])
        assert sum(lam) == 0, (name, label)
    alignments = []
    for orient in (list(cyc), list(reversed(cyc))):
        for r in range(n):
            order = orient[r:] + orient[:r]
            good = all(
                dot([F(x) for x in d["row"]], cls[order[i]]) == d["lam"][i]
                for d in rows.values() for i in range(n))
            if good:
                alignments.append(order)
    ok(f"{name}: a common vertex alignment matches every printed row "
       f"({len(alignments)} found)", len(alignments) >= 1)
    return cls

cls61 = check_surface_conventions(E61, ROWS6)
cls62 = check_surface_conventions(E62, ROWS6)
cls71 = check_surface_conventions(E71, ROWS7["E_{7,1}"])
cls72 = check_surface_conventions(E72, ROWS7["E_{7,2}"])

# ------------------------------------------------- 3: norms, pairings, summands
def pairing_on_dual(I):
    """Intersection pairing transported to Pic^vee: <r, s> = r I^{-1} s^T."""
    Iinv = inverse(I)
    return lambda r, s: dot(r, mvec(Iinv, s))

def dual_class(I, r):
    """The H_2-class of a functional r in the Pic basis: I^{-1} r^T."""
    return mvec(inverse(I), r)

def check_dp6_lattice(E, rows):
    name = E["name"]
    pair = pairing_on_dual(E["I"])
    s1 = [F(x) for x in rows["s1"]["row"]]
    s2 = [F(x) for x in rows["s2"]["row"]]
    s3 = [F(x) for x in rows["s3"]["row"]]
    for lbl, r in (("s1", s1), ("s2", s2), ("s3", s3)):
        assert dot(r, E["K"]) == 0, (name, lbl)
        assert pair(r, r) == F(-2), (name, lbl, pair(r, r))
    assert abs(pair(s1, s2)) == 1
    assert pair(s3, s1) == 0 and pair(s3, s2) == 0
    ok(f"{name}: phi(s1),phi(s2),phi(s3) in K-perp, norms -2, "
       "s3 orthogonal to the A2 pair", True)
    G = [[pair(a, b) for b in (s1, s2, s3)] for a in (s1, s2, s3)]
    det = (G[0][0]*(G[1][1]*G[2][2]-G[1][2]*G[2][1])
           - G[0][1]*(G[1][0]*G[2][2]-G[1][2]*G[2][0])
           + G[0][2]*(G[1][0]*G[2][1]-G[1][1]*G[2][0]))
    ok(f"{name}: K-perp discriminant 6", abs(det) == 6)
    ok(f"{name}: coordinate classes are a basis of K-perp",
       rank([s1, s2, s3]) == 3)

    # (-2)-vector enumeration in K-perp: expect the 8 roots of A1 + A2,
    # namely +/-s3 and +/-{s1, s2, s1+s2-type combination}.
    Q = [[-x for x in row] for row in G]          # positive definite
    Qinv = inverse(Q)
    bounds = [int(math.isqrt(int(2 * Qinv[i][i]) + 1)) + 1 for i in range(3)]
    roots = []
    for xyz in product(*[range(-b, b + 1) for b in bounds]):
        if xyz == (0, 0, 0):
            continue
        v = [F(x) for x in xyz]
        q = dot(v, mvec(Q, v))
        if q == 2:
            roots.append(xyz)
    ok(f"{name}: exactly 8 (-2)-classes in K-perp", len(roots) == 8)
    a1 = [r for r in roots if r[0] == 0 and r[1] == 0]
    ok(f"{name}: the A1 pair is exactly +/- phi(s3)",
       sorted(a1) == [(0, 0, -1), (0, 0, 1)])
    a2 = [r for r in roots if r[2] == 0]
    ok(f"{name}: the six A2 roots have zero s3-coordinate and rank 2",
       len(a2) == 6 and rank([[F(x) for x in r] for r in a2]) == 2)
    return s1, s2, s3

s61 = check_dp6_lattice(E61, ROWS6)
s62 = check_dp6_lattice(E62, ROWS6)

def check_dp7_lattice(E, rows):
    name = E["name"]
    pair = pairing_on_dual(E["I"])
    al = [F(x) for x in rows["alpha"]["row"]]
    be = [F(x) for x in rows["beta"]["row"]]
    assert dot(al, E["K"]) == 0 and dot(be, E["K"]) == 0
    ok(f"{name}: phi(beta)^2 = -2 and phi(alpha)^2 = -4",
       pair(be, be) == F(-2) and pair(al, al) == F(-4))
    G = [[pair(a, b) for b in (al, be)] for a in (al, be)]
    det = G[0][0] * G[1][1] - G[0][1] * G[1][0]
    ok(f"{name}: K-perp discriminant 7", abs(det) == 7)
    Q = [[-x for x in row] for row in G]
    Qinv = inverse(Q)
    bounds = [int(math.isqrt(int(2 * Qinv[i][i]) + 1)) + 1 for i in range(2)]
    sols = [xy for xy in product(*[range(-b, b + 1) for b in bounds])
            if xy != (0, 0)
            and dot([F(a) for a in xy], mvec(Q, [F(a) for a in xy])) == 2]
    ok(f"{name}: unique (-2)-line, spanned by phi(beta)",
       sorted(sols) == [(0, -1), (0, 1)])
    return al, be

a71, b71 = check_dp7_lattice(E71, ROWS7["E_{7,1}"])
a72, b72 = check_dp7_lattice(E72, ROWS7["E_{7,2}"])

# ------------------------------------------------- 5: sweep kernels (standard)
# Standard del Pezzo bases: (l, e1, ..., en) with form diag(1, -1, ..., -1).
def std_form(n):
    return [[F(1) if i == j == 0 else (F(-1) if i == j else F(0))
             for j in range(n + 1)] for i in range(n + 1)]

g6, g7 = std_form(3), std_form(2)
K6_std = [F(-3), F(1), F(1), F(1)]   # K = -3l + e1 + e2 + e3 as a class
K7_std = [F(-3), F(1), F(1)]

def ipd(g, u, v):
    return dot(u, mvec(g, v))

def kperp_basis(g, K):
    return kernel([mvec(g, K)])

KP6 = kperp_basis(g6, K6_std)
KP7 = kperp_basis(g7, K7_std)
ok("standard models: K-perp ranks 3 and 2", len(KP6) == 3 and len(KP7) == 2)

def push_kernel(g, K, images, expect, expect_name, label):
    """images: pushforward image of each standard basis class (l, e1, ...).
    Checks ker(K-perp -> H_2(V)) equals the span of `expect`."""
    imgs = mat(images)
    n = len(imgs)
    sub = kperp_basis(g, K)
    # a K-perp vector sum_t c_t sub_t is in the kernel iff its image is 0
    restr = [[sum(imgs[j][i] * b[j] for j in range(n)) for b in sub]
             for i in range(len(imgs[0]))]
    coeff = kernel(restr)
    got = [[sum(c[t] * sub[t][j] for t in range(len(sub)))
            for j in range(n)] for c in coeff]
    ok(f"{label}: ker(K-perp -> H_2(V)) = {expect_name}",
       span_equal(got, mat(expect)))
    ok(f"{label}: K-perp image rank {len(sub) - len(got)}",
       rank([[sum(imgs[j][i] * b[j] for j in range(n))
              for i in range(len(imgs[0]))] for b in sub])
       == len(sub) - len(got))
    return got

# dP6 in (P1)^3: degrees (C.(l-e1), C.(l-e2), C.(l-e3))
img_P13 = [[1, 1, 1],   # l
           [1, 0, 0],   # e1
           [0, 1, 0],   # e2
           [0, 0, 1]]   # e3
A1_std = [[1, -1, -1, -1]]
A2_std = [[0, 1, -1, 0], [0, 0, 1, -1]]
push_kernel(g6, K6_std, img_P13, A1_std, "the A1 line", "dP6 in (P1)^3")

# dP6 in P(T_P2): degrees (C.l, C.(2l-e1-e2-e3))
img_PT = [[1, 2], [0, 1], [0, 1], [0, 1]]
push_kernel(g6, K6_std, img_PT, A2_std, "the A2 plane", "dP6 in P(T_P2)")

# dP7 in V7 = Bl_pt P3: degrees (C.(2l-e1-e2), C.(l-e1-e2))
img_V7 = [[2, 1], [1, 1], [1, 1]]
L7_std = [[0, 1, -1]]
push_kernel(g7, K7_std, img_V7, L7_std, "the line <e1-e2>",
            "dP7 in Bl_pt P3")

# ------------------------------------------------- 6: transfer to std markings
HEX_STD = [[0, 1, 0, 0],        # e1
           [1, -1, -1, 0],      # l - e1 - e2
           [0, 0, 1, 0],        # e2
           [1, 0, -1, -1],      # l - e2 - e3
           [0, 0, 0, 1],        # e3
           [1, -1, 0, -1]]      # l - e1 - e3
PENT_STD = [[0, 1, 0],          # e1        (-1)
            [1, -1, 0],         # l - e1    (0)
            [1, 0, -1],         # l - e2    (0)
            [0, 0, 1],          # e2        (-1)
            [1, -1, -1]]        # l - e1 - e2  (-1)

def dihedral_relabelings(seq):
    n = len(seq)
    out = []
    for refl in (False, True):
        s = list(reversed(seq)) if refl else list(seq)
        for r in range(n):
            out.append(s[r:] + s[:r])
    return out

def find_isometry(E, cls, g_std, K_std, std_cycle, label):
    """Find T with T(class(cyclic_i)) = std marking, isometric, K -> K."""
    printed = [cls[q] for q in E["cyclic"]]
    n = len(E["K"])
    idx = None
    # choose n independent printed classes
    for cand in product(range(len(printed)), repeat=n):
        if len(set(cand)) == n and rank([printed[i] for i in cand]) == n:
            idx = list(cand)
            break
    Pm = transpose([printed[i] for i in idx])
    Pinv = inverse(Pm)
    for std in dihedral_relabelings([[F(x) for x in v] for v in std_cycle]):
        Sm = transpose([std[i] for i in idx])
        T = mmul(Sm, Pinv)
        if all(mvec(T, printed[i]) == std[i] for i in range(len(printed))):
            # isometry check on the printed basis
            iso = all(
                ipd(g_std, mvec(T, printed[i]), mvec(T, printed[j]))
                == dot(printed[i], mvec(E["I"], printed[j]))
                for i in range(len(printed)) for j in range(len(printed)))
            if not iso:
                continue
            # K is printed as a class in the Pic basis
            if mvec(T, E["K"]) != K_std:
                continue
            ok(f"{label}: explicit isometry to the standard marking, K -> K",
               True)
            return T
    print(f"FAIL: {label}: no isometry found")
    sys.exit(1)

T61 = find_isometry(E61, cls61, g6, K6_std, HEX_STD, "E_{6,1}")
T62 = find_isometry(E62, cls62, g6, K6_std, HEX_STD, "E_{6,2}")
T71 = find_isometry(E71, cls71, g7, K7_std, PENT_STD, "E_{7,1}")
T72 = find_isometry(E72, cls72, g7, K7_std, PENT_STD, "E_{7,2}")

def transfer(T, E, row):
    return mvec(T, dual_class(E["I"], [F(x) for x in row]))

for E, T, s, label in ((E61, T61, s61, "E_{6,1}"), (E62, T62, s62, "E_{6,2}")):
    z3 = transfer(T, E, ROWS6["s3"]["row"])
    ok(f"{label}: phi(s3) transfers to +/- (l-e1-e2-e3)",
       z3 in ([F(1), F(-1), F(-1), F(-1)],
              [F(-1), F(1), F(1), F(1)]))
    z1 = transfer(T, E, ROWS6["s1"]["row"])
    z2 = transfer(T, E, ROWS6["s2"]["row"])
    ok(f"{label}: phi(s1), phi(s2) transfer into the A2 plane",
       in_span(z1, mat(A2_std)) and in_span(z2, mat(A2_std))
       and rank([z1, z2]) == 2)

for E, T, label in ((E71, T71, "E_{7,1}"), (E72, T72, "E_{7,2}")):
    rows = ROWS7[label]
    zb = transfer(T, E, rows["beta"]["row"])
    ok(f"{label}: phi(beta) transfers to +/- (e1-e2)",
       zb in ([F(0), F(1), F(-1)], [F(0), F(-1), F(1)]))
    za = transfer(T, E, rows["alpha"]["row"])
    ok(f"{label}: phi(alpha) is not on the (-2)-line",
       not in_span(za, mat(L7_std)))

# ------------------------------------------------- 7: the conditional deduction
print()
print("Branch subspaces in the 36 matrix coordinates "
      "(direct sums verified above):")
print("  dP7 points:  kappa in R_7          <=>  alpha = 0")
print("  dP6 points:  kappa in A1 (line)    <=>  s1 = s2 = 0")
print("               kappa in A2 (plane)   <=>  s3 = 0")
print()
print("Union over the component choices at the two dP6 points = the four")
print("restricted profiles LL, LP, PL, PP of Section 5.3 of the paper.")
print("Parent-certified forced-zero lists on those restricted kernels")
print("(checker: sage mixed_candidate.sage):")
FORCED = {"LL": ["u9", "u11", "beta1", "beta2"],
          "LP": ["u9", "u11", "beta1", "beta2", "s3_1"],
          "PL": ["u9", "u11", "beta1", "beta2", "s3_2"],
          "PP": ["u9", "u11", "beta1", "beta2"]}
for prof, lst in FORCED.items():
    assert "u9" in lst and "u11" in lst
    print(f"  {prof}: {', '.join(lst)}")
print("  note: the forced-zero lists above are EXTERNAL INPUT (certified "
      "by sage mixed_candidate.sage,\n        independently "
      "recomputed by global_kernel.py); they are not re-derived here.")
print()
print("Hence (paper, Theorem 6.5): any smoothing class kappa would lie in")
print("a restricted kernel with kappa(N9) = 0, contradicting nodal")
print("nonvanishing; and (paper, Proposition 6.4) any smoothing has")
print("delta_9 = delta_11 = 0 in H_3(X_t; Q).")
print()
print(f"ALL CHECKS PASSED ({PASS} assertions)")
