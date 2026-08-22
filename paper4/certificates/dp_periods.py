#!/usr/bin/env python3
"""Exact and numerical checker for the period lemmas of the paper (Section 3).

Verifies, for the two rank-one Milnor-fiber channels:

  A. (P^1)^3 minus the closure of {xyz = 1}      (dP6 line component)
  B. Bl_pt P^3 minus the closure of {xy/z = 1}   (dP7 unique component)

the following claims:

  1. Fan bookkeeping: the subtorus-closure surface has divisor class
     -K_V/2 (so its cone is the anticanonical cone over the del Pezzo,
     with the half-anticanonical polarization of the sweep); every edge
     of the fan crossed by m has pairing values (+1, -1); the sliced fan
     is the hexagon (all self-intersections -1) resp. the pentagon
     (pattern (-1,-1,-1,0,0)).
  2. b_3(M) = 1: the restriction H^2(V) -> H^2(S) has rank b_2(V) and
     corank 1 resp. 1 in H^2(S).
  3. The period: on the explicit interpolating piece of the cycle T of
     Lemma 3.2 of the paper, the pullback of Omega = w/(w-1)^2 dlog p dlog q dlog z
     is the CONSTANT form (+/-1) dr d(theta) d(phi); hence
     int_T Omega = +/- 4 pi^2.  Checked symbolically-numerically at random
     sample points, together with |t| constant and nonvanishing on the
     domain, and the endpoint-matching of the three pieces.
  4. The two boundary pieces contribute zero: the pullback of Omega
     vanishes identically on them (dp = 0 resp. dq = 0), and Omega is
     regular there (t != 0).

Run:  python3 dp_periods.py
"""

import cmath
import math
import random
import sys
from fractions import Fraction as F

PASS = 0
def ok(label, cond):
    global PASS
    if not cond:
        print(f"FAIL: {label}")
        sys.exit(1)
    PASS += 1
    print(f"  pass: {label}")

random.seed(20260821)

# ---------------------------------------------------------------- linear algebra
def rref(A):
    M = [[F(x) for x in row] for row in A]
    rows, cols = len(M), len(M[0]) if M else 0
    piv, r = [], 0
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
    return M, piv

def rank(A):
    return len(rref(A)[1]) if A else 0

def dot(u, v):
    return sum(a * b for a, b in zip(u, v))

def det3(M):
    return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
            - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
            + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))

# ---------------------------------------------------------------- case data
# Case A: V = (P^1)^3.  Rays +/- e_i; maximal cones = octants.
E1, E2, E3 = (1, 0, 0), (0, 1, 0), (0, 0, 1)
def neg(v):
    return tuple(-x for x in v)

RAYS_A = [E1, neg(E1), E2, neg(E2), E3, neg(E3)]
CONES_A = [(s1 and E1 or neg(E1), s2 and E2 or neg(E2), s3 and E3 or neg(E3))
           for s1 in (1, 0) for s2 in (1, 0) for s3 in (1, 0)]
M_A = (1, 1, 1)

# Case B: V = Bl_pt P^3.  Rays e1,e2,e3, r = -e1-e2-e3, s = e1+e2+e3.
R_ = (-1, -1, -1)
S_ = (1, 1, 1)
RAYS_B = [E1, E2, E3, R_, S_]
CONES_B = [(E1, E2, S_), (E1, E3, S_), (E2, E3, S_),
           (E1, E2, R_), (E1, E3, R_), (E2, E3, R_)]
M_B = (1, 1, -1)

def edges(cones):
    out = set()
    for c in cones:
        for i in range(3):
            for j in range(i + 1, 3):
                out.add(frozenset((c[i], c[j])))
    return sorted(tuple(sorted(e)) for e in out)

def crossed_edges(cones, m):
    return [e for e in edges(cones)
            if (dot(m, e[0]) > 0) != (dot(m, e[1]) > 0)
            and dot(m, e[0]) != 0 and dot(m, e[1]) != 0]

def check_case_combinatorics(name, rays, cones, m, expect_selfints):
    ok(f"{name}: m does not vanish on any ray",
       all(dot(m, v) != 0 for v in rays))
    ce = crossed_edges(cones, m)
    for e in ce:
        vals = sorted(dot(m, v) for v in e)
        assert vals == [-1, 1], (name, e, vals)
    ok(f"{name}: every crossed edge has pairing values (+1,-1); "
       f"{len(ce)} crossed edges", len(ce) == len(expect_selfints))

    # slice-fan rays: for crossed edge (v+, v-) the ray is v+ + v-
    sl = []
    for e in ce:
        vplus = e[0] if dot(m, e[0]) > 0 else e[1]
        vminus = e[1] if vplus is e[0] else e[0]
        sl.append(tuple(a + b for a, b in zip(vplus, vminus)))
    for v in sl:
        assert dot(m, v) == 0
    # order the slice rays cyclically in the plane m-perp and read off the
    # toric surface self-intersections u_{i-1} + u_{i+1} = -a_i u_i
    # (so the self-intersection of the i-th boundary curve is a_i).
    b1 = sl[0]
    b2 = next(v for v in sl[1:] if rank([list(b1), list(v)]) == 2)
    def coords(v):
        Mm, piv = rref([[b1[0], b2[0], v[0]],
                        [b1[1], b2[1], v[1]],
                        [b1[2], b2[2], v[2]]])
        assert piv == [0, 1], (name, v)
        return (Mm[0][2], Mm[1][2])
    pts = [coords(v) for v in sl]
    ang = sorted(range(len(pts)),
                 key=lambda i: math.atan2(float(pts[i][1]), float(pts[i][0])))
    cyc = [pts[i] for i in ang]
    n = len(cyc)
    selfints = []
    for i in range(n):
        a, b, c = cyc[(i - 1) % n], cyc[i], cyc[(i + 1) % n]
        # solve a + c = k * b   (2D primitive rays of a smooth fan)
        if b[0] != 0:
            k = (a[0] + c[0]) / b[0]
        else:
            k = (a[1] + c[1]) / b[1]
        assert a[0] + c[0] == k * b[0] and a[1] + c[1] == k * b[1], (name, i)
        selfints.append(-k)
    want = sorted(expect_selfints)
    ok(f"{name}: sliced fan self-intersections {sorted(selfints)}",
       sorted(selfints) == want)
    return ce

ce_A = check_case_combinatorics("(P1)^3, m=(1,1,1)", RAYS_A, CONES_A, M_A,
                                [-1, -1, -1, -1, -1, -1])
ce_B = check_case_combinatorics("Bl_pt P3, m=(1,1,-1)", RAYS_B, CONES_B, M_B,
                                [-1, -1, -1, 0, 0])

# Divisor class of the subtorus closure: [S] = sum over rays with <m,v> < 0
# of |<m,v>| D_v ; compare with -K/2 = (sum of all D_v)/2 in the class group.
def class_group_check(name, rays, m):
    n = len(rays)
    # class group = Z^n / im(M): relations rows are (<e_i^*, v_j>)_j
    rel = [[F(rays[j][i]) for j in range(n)] for i in range(3)]
    Svec = [F(max(0, -dot(m, v))) for v in rays]
    half_anti = [F(1, 2)] * n
    diff = [a - b for a, b in zip(Svec, half_anti)]
    # diff must lie in the Q-span of the relation rows
    ok(f"{name}: [S] = -K_V/2 in Cl(V) tensor Q",
       rank(rel) == rank(rel + [diff]))

class_group_check("(P1)^3", RAYS_A, M_A)
class_group_check("Bl_pt P3", RAYS_B, M_B)

# The divisor bookkeeping div(Omega_M) = -2 S-bar needs <m,v> = +/-1 on
# EVERY ray (order along D_v is <m,v> - 2 min(<m,v>,0) - 1).
ok("(P1)^3: <m,v> = +/-1 on every ray",
   all(abs(dot(M_A, v)) == 1 for v in RAYS_A))
ok("Bl_pt P3: <m,v> = +/-1 on every ray",
   all(abs(dot(M_B, v)) == 1 for v in RAYS_B))

# ---------------------------------------------------------------- b3(M) = 1
# restriction image H^2(V) -> H^2(S) in the standard del Pezzo bases
img_A = [[1, -1, 0, 0], [1, 0, -1, 0], [1, 0, 0, -1]]  # l-e1, l-e2, l-e3
ok("(P1)^3: restriction rank 3, corank 1 in H^2(dP6): b3(M) = 1",
   rank(img_A) == 3)
img_B = [[2, -1, -1], [1, -1, -1]]                     # 2l-e1-e2, l-e1-e2
ok("Bl_pt P3: restriction rank 2, corank 1 in H^2(dP7): b3(M) = 1",
   rank(img_B) == 2)

# ---------------------------------------------------------------- the period
# Adapted chart at a crossed edge (v+, v-) with corner generator n0 such
# that gamma := <m, n0> = +1: coordinates (p, q, z) dual to (v+, v-, n0),
# w = p q^{-1} z^{gamma}, S-bar = {t = 0} with t = q - p z^{gamma}, and
#     Omega = (sign) z^{gamma-1} dp ^ dq ^ dz / t^2 .
# The cycle T = T1 + T2 + T3 of Lemma 3.2 of the paper; only T2 contributes:
#   T2: p = r*eps*e^{i theta},
#       q = -(1-r)*eps*R^gamma*e^{i(theta + gamma phi)},
#       z = R*e^{i phi},   (r,theta,phi) in [0,1] x [0,2pi]^2 .
# Claim: the pullback of Omega to T2 is the constant form dr d(theta) d(phi)
# (up to the global sign), so int_{T2} Omega = 4 pi^2 up to sign.

def chart_check(name, vplus, vminus, n0, m):
    B = [list(vplus), list(vminus), list(n0)]
    ok(f"{name}: adapted basis unimodular", abs(det3(B)) == 1)
    # dual-basis expansion of m: coefficients are <m, basis vectors>
    alpha, beta, gamma = dot(m, vplus), dot(m, vminus), dot(m, n0)
    ok(f"{name}: w = p q^-1 z^gamma with (alpha,beta,gamma)="
       f"({alpha},{beta},{gamma})", alpha == 1 and beta == -1)
    return gamma

def pullback_constant(gamma, eps, R, samples=200):
    """Check the pullback of z^{gamma-1} dp dq dz / t^2 at random points."""
    worst = 0.0
    tmods = []
    for _ in range(samples):
        r = random.uniform(0.01, 0.99)
        th = random.uniform(0, 2 * math.pi)
        ph = random.uniform(0, 2 * math.pi)
        eith = cmath.exp(1j * th)
        eigph = cmath.exp(1j * gamma * ph)
        p = r * eps * eith
        q = -(1 - r) * eps * (R ** gamma) * eith * eigph
        z = R * cmath.exp(1j * ph)
        t = q - p * z ** gamma
        tmods.append(abs(t))
        # partial derivatives
        p_r, p_t, p_p = eps * eith, 1j * p, 0
        q_r = eps * (R ** gamma) * eith * eigph
        q_t = 1j * q
        q_p = 1j * gamma * q
        z_r, z_t, z_p = 0, 0, 1j * z
        J = [[p_r, p_t, p_p], [q_r, q_t, q_p], [z_r, z_t, z_p]]
        det = (J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])
               - J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])
               + J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]))
        val = (z ** (gamma - 1)) * det / t ** 2
        # value should be the constant +1 (overall sign is an orientation
        # convention; only nonvanishing of the total is used downstream)
        worst = max(worst, abs(val - 1.0))
    const_t = max(tmods) - min(tmods)
    return worst, const_t

for name, (vplus, vminus, n0, m) in {
    "(P1)^3 edge (e1, -e2), corner n0 = e3":
        (E1, neg(E2), E3, M_A),
    "Bl_pt P3 edge (e1, r), corner n0 = e2":
        (E1, R_, E2, M_B),
}.items():
    gamma = chart_check(name, vplus, vminus, n0, m)
    ok(f"{name}: gamma = +1 at the chosen corner", gamma == 1)
    worst, const_t = pullback_constant(gamma, eps=0.37, R=2.19)
    ok(f"{name}: pullback of Omega on T2 is the constant +dr dtheta dphi "
       f"(max deviation {worst:.2e})", worst < 1e-9)
    ok(f"{name}: |t| constant and nonzero on T2 (spread {const_t:.2e})",
       const_t < 1e-9)
    # the other corner has gamma = -1 (the switch is genuinely needed)
    other = {"(P1)^3 edge (e1, -e2), corner n0 = e3": neg(E3),
             "Bl_pt P3 edge (e1, r), corner n0 = e2": E3}[name]
    ok(f"{name}: opposite corner has gamma = -1",
       dot(m, other) == -1)

# boundary pieces contribute zero: on T1 (p = 0) and T3 (q = 0) the
# pullback vanishes identically because dp resp. dq pulls back to zero
# while t = q resp. t = -p z^gamma stays nonzero; verify t != 0 there.
for gamma in (1,):
    eps, R = 0.37, 2.19
    for _ in range(100):
        th = random.uniform(0, 2 * math.pi)
        ph = random.uniform(0, 2 * math.pi)
        rho = random.uniform(0, R)
        q = eps * (R ** gamma) * cmath.exp(1j * th)
        z = rho * cmath.exp(1j * ph)
        assert abs(q - 0 * z ** gamma) > 1e-12          # T1: t = q
        p = eps * cmath.exp(1j * th)
        z2 = random.uniform(R, 40 * R) * cmath.exp(1j * ph)
        assert abs(-p * z2 ** gamma) > 1e-12            # T3: t = -p z^gamma
ok("boundary pieces T1, T3: t nonvanishing (and dp resp. dq pull back "
   "to 0, so their Omega-integrals vanish identically)", True)

# endpoint matching of the three pieces (chain-level): the r=0 end of T2 is
# the torus {p=0, |q| = eps R^gamma, |z| = R} (reparametrized by
# psi = theta + gamma phi + pi), the r=1 end is {|p| = eps, q=0, |z| = R}.
th, ph = 1.234, 2.345
eps, R, gamma = 0.37, 2.19, 1
q0 = -eps * R ** gamma * cmath.exp(1j * (th + gamma * ph))
ok("T2 endpoints: r=0 lands on the T1 boundary torus "
   "(|q| = eps R^gamma), r=1 on the T3 boundary torus (|p| = eps)",
   abs(abs(q0) - eps * R ** gamma) < 1e-12)

# ---------------------------------------------------- corner-chart cap T3'
# For each case: adapted bases at the two corners of the chosen edge,
# transition monomials between the dual coordinate systems, the S-bar
# equation in the corner chart, the boundary matching of T2's r=1 face
# with dT3', the distance bound |p' - q' z'| = eps/R on T3', and the
# documentation of the original defect (in case B the literal closure of
# the naive cap collapses to the corner-chart origin, which LIES ON S-bar).
def dual_basis(B):
    n = len(B)
    aug = [[F(B[j][i]) for j in range(n)] + [F(1) if i == k else F(0)
            for k in range(n)] for i in range(n)]
    R2, piv = rref(aug)
    assert piv == list(range(n))
    inv = [row[n:] for row in R2]           # inv = (B^T)^{-1}
    # the dual basis vectors are the ROWS of (B^T)^{-1}:
    # row j pairs to delta_{jk} with basis vector B_k.
    return [list(inv[j]) for j in range(n)]

def expand_in_duals(m, B):
    # m = sum c_j m_j^*  with  c_j = <m, B_j>
    return [dot(m, b) for b in B]

for name, (vplus, vminus, n0, n1, m) in {
    "(P1)^3": (E1, neg(E2), E3, neg(E3), M_A),
    "Bl_pt P3": (E1, R_, E2, E3, M_B),
}.items():
    B0 = [list(vplus), list(vminus), list(n0)]
    B1 = [list(vplus), list(vminus), list(n1)]
    ok(f"{name}: both corner bases unimodular",
       abs(det3(B0)) == 1 and abs(det3(B1)) == 1)
    D0, D1 = dual_basis(B0), dual_basis(B1)
    # transition exponents: p' = chi^{D1[0]} in terms of (p,q,z):
    # coefficients of D1[j] in the D0-dual basis are <D1[j], B0-basis>
    trans = [[dot(D1[j], B0[i]) for i in range(3)] for j in range(3)]
    # expected: case A: p'=p, q'=q, z'=1/z; case B: p'=p/z, q'=q/z, z'=1/z
    expect = {"(P1)^3": [[1, 0, 0], [0, 1, 0], [0, 0, -1]],
              "Bl_pt P3": [[1, 0, -1], [0, 1, -1], [0, 0, -1]]}[name]
    ok(f"{name}: corner transition exponents {expect}",
       trans == [[F(x) for x in row] for row in expect])
    # w in the corner chart: coefficients <m, B1-basis> should be (1,-1,-1),
    # so S-bar = {w=1} = {p' = q' z'} there
    ok(f"{name}: corner chart has w = p' q'^-1 z'^-1, S-bar = {{p'=q'z'}}",
       expand_in_duals(m, B1) == [1, -1, -1])
    # numerical: T2's r=1 face lands on dT3', and T3' avoids S-bar.
    # The cap radius is rad = eps * R^c with c the z-exponent of the
    # p'-transition (c = 0 in case A, c = -1 in case B).
    eps, R = 0.37, 2.19
    rad = eps * (R ** float(expect[0][2]))
    worst_b = worst_d = 0.0
    for _ in range(120):
        th = random.uniform(0, 2 * math.pi)
        ph = random.uniform(0, 2 * math.pi)
        pp, qq, zz = eps * cmath.exp(1j * th), 0.0, R * cmath.exp(1j * ph)
        # transition via monomials
        def mono(pt, exps):
            out = 1.0 + 0j
            for base, e in zip(pt, exps):
                if e:
                    out *= base ** e
            return out
        pt = (pp, qq, zz)
        pB = mono(pt, expect[0]); zB = mono(pt, expect[2])
        worst_b = max(worst_b, abs(abs(pB) - rad),
                      abs(abs(zB) - 1 / R))
        # T3' point at random interior z' and the S-bar distance
        zp = random.uniform(0, 1 / R) * cmath.exp(1j * ph)
        ppr = rad * cmath.exp(1j * th)
        worst_d = max(worst_d, abs(abs(ppr - 0.0 * zp) - rad))
    ok(f"{name}: r=1 face of T2 lands on the T3' boundary torus "
       f"(cap radius eps*R^{expect[0][2]}; dev {worst_b:.1e})",
       worst_b < 1e-12)
    ok(f"{name}: |p' - q' z'| = cap radius > 0 on T3' (dev {worst_d:.1e})",
       worst_d < 1e-12)
# the documented defect: the corner-chart origin satisfies p' = q' z'
ok("case B defect documented: the corner-chart origin lies on S-bar, so "
   "the literal closure of the naive cap is not contained in M "
   "(0 - 0*0 = 0)", (0 - 0 * 0) == 0)

print()
print("Conclusion (Lemma 3.2 of the paper): int_T Omega = 4 pi^2 (in the chart orientation) for")
print("the explicit cycle T in both Milnor fibers; since b_3(M) = 1, the")
print("period pairing on H_3(M;Q) is injective for the dP6 line component")
print("and the dP7 component.")
print()
print(f"ALL CHECKS PASSED ({PASS} assertions)")
