#!/usr/bin/env python3
"""Independent rebuild of the mixed obstruction matrix and its forced zeros.

This script re-derives, in pure Python over exact rationals and WITHOUT the
any Sage code, every finite input consumed by the nonsmoothability
theorem of the paper (Theorem 6.5):

  1. the 26 exceptional-curve rows from the printed rays and diagonals of
     Table 1 of the paper (the square completions are re-derived from the
     embedded maximal cones and re-checked against the unit-square identity
     q_a + q_c = q_b + q_d);
  2. the node-subsystem invariants: rank 18 with coloops
     exactly N_9 and N_11;
  3. the ten divisorial rows from the printed lambda-functionals of
     Table 2 of the paper and the embedded restriction columns;
  4. the full 36 x 26 matrix: rank 21, kernel dimension 15, and the unique
     unrestricted row coloop D_{7,1}:beta;
  5. the four branch-restricted kernels (dimensions 9, 10, 10, 12) and the
     forced-zero lists: u9, u11, beta1, beta2 vanish identically on all
     four, plus s3 on the mixed profiles;
  6. sanity: the forced-zero test is not vacuous (u1 is not forced).

Coordinate order (matching Section 5.3 of the paper):
  (u_1..u_26; s_{1,1}, s_{2,1}, s_{3,1}; s_{1,2}, s_{2,2}, s_{3,2};
   alpha_1, beta_1; alpha_2, beta_2).

Run:  python3 global_kernel.py
"""

from fractions import Fraction as F
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
        if r == rows:
            break
    return M, piv

def rank(A):
    return len(rref(A)[1]) if A else 0

def kernel(A):
    M, piv = rref(A)
    cols = len(A[0])
    free = [c for c in range(cols) if c not in piv]
    out = []
    for fc in free:
        v = [F(0)] * cols
        v[fc] = F(1)
        for r, pc in enumerate(piv):
            v[pc] = -M[r][fc]
        out.append(v)
    return out

# ---------------------------------------------------------------- printed data
# Table 1 of the paper: the 30 rays (1-based index -> vector).
Q = {
    1: (-1, -1, -1, -1), 2: (-1, -1, -1, 0), 3: (-1, -1, 0, -1),
    4: (-1, -1, 0, 1), 5: (-1, -1, 1, 0), 6: (-1, -1, 1, 1),
    7: (-1, 0, -1, -1), 8: (-1, 0, -1, 0), 9: (-1, 0, 0, -1),
    10: (-1, 0, 0, 1), 11: (-1, 0, 1, 0), 12: (-1, 0, 1, 1),
    13: (0, -1, -1, -1), 14: (0, -1, -1, 0), 15: (0, -1, 0, 1),
    16: (0, -1, 0, 0), 17: (0, 1, -1, -1), 18: (0, 1, -1, 0),
    19: (0, 1, 0, -1), 20: (0, 1, 0, 0), 21: (1, 1, -1, -1),
    22: (1, 0, -1, 0), 23: (0, 0, 0, -1), 24: (0, 0, 0, 1),
    25: (0, 0, 1, 0), 26: (0, 0, 1, 1), 27: (-1, -1, 0, 0),
    28: (-1, 0, 0, 0), 29: (0, 0, -1, -1), 30: (0, 0, -1, 0),
}

MAXIMAL_CONES = (
    (1, 2, 7, 27), (1, 2, 7, 29), (1, 2, 13, 27), (1, 2, 13, 29),
    (1, 3, 7, 27), (1, 3, 7, 29), (1, 3, 13, 27), (1, 3, 13, 29),
    (2, 4, 8, 27), (2, 4, 8, 30), (2, 4, 14, 27), (2, 4, 14, 30),
    (2, 7, 8, 27), (2, 7, 8, 29), (2, 8, 29, 30), (2, 13, 14, 27),
    (2, 13, 14, 29), (2, 14, 29, 30), (3, 5, 9, 23), (3, 5, 9, 27),
    (3, 5, 13, 23), (3, 5, 13, 27), (3, 7, 9, 27), (3, 7, 9, 29),
    (3, 9, 23, 29), (3, 13, 23, 29), (4, 6, 10, 15), (4, 6, 10, 27),
    (4, 6, 15, 27), (4, 8, 10, 27), (4, 8, 10, 30), (4, 10, 15, 30),
    (4, 14, 15, 27), (4, 14, 15, 30), (5, 6, 11, 25), (5, 6, 11, 27),
    (5, 6, 16, 25), (5, 6, 16, 27), (5, 9, 11, 23), (5, 9, 11, 27),
    (5, 11, 23, 25), (5, 13, 16, 23), (5, 13, 16, 27), (5, 16, 23, 25),
    (6, 10, 12, 15), (6, 10, 12, 27), (6, 11, 12, 25), (6, 11, 12, 27),
    (6, 12, 15, 26), (6, 12, 25, 26), (6, 15, 16, 26), (6, 15, 16, 27),
    (6, 16, 25, 26), (7, 8, 17, 28), (7, 8, 17, 29), (7, 8, 27, 28),
    (7, 9, 17, 28), (7, 9, 17, 29), (7, 9, 27, 28), (8, 10, 18, 28),
    (8, 10, 18, 30), (8, 10, 27, 28), (8, 17, 18, 28), (8, 17, 18, 29),
    (8, 18, 29, 30), (9, 11, 19, 23), (9, 11, 19, 28), (9, 11, 27, 28),
    (9, 17, 19, 28), (9, 17, 19, 29), (9, 19, 23, 29), (10, 12, 15, 24),
    (10, 12, 18, 24), (10, 12, 18, 28), (10, 12, 27, 28), (10, 15, 24, 30),
    (10, 18, 24, 30), (11, 12, 19, 25), (11, 12, 19, 28), (11, 12, 27, 28),
    (11, 19, 23, 25), (12, 15, 24, 26), (12, 18, 20, 24), (12, 18, 20, 28),
    (12, 19, 20, 25), (12, 19, 20, 28), (12, 20, 24, 26), (12, 20, 25, 26),
    (13, 14, 16, 22), (13, 14, 16, 27), (13, 14, 22, 29), (13, 16, 22, 23),
    (13, 21, 22, 23), (13, 21, 22, 29), (13, 21, 23, 29), (14, 15, 16, 22),
    (14, 15, 16, 27), (14, 15, 22, 30), (14, 22, 29, 30), (15, 16, 22, 26),
    (15, 22, 24, 26), (15, 22, 24, 30), (16, 22, 23, 25), (16, 22, 25, 26),
    (17, 18, 19, 21), (17, 18, 19, 28), (17, 18, 21, 29), (17, 19, 21, 29),
    (18, 19, 20, 21), (18, 19, 20, 28), (18, 20, 21, 24), (18, 21, 22, 24),
    (18, 21, 22, 30), (18, 21, 29, 30), (18, 22, 24, 30), (19, 20, 21, 25),
    (19, 21, 23, 25), (19, 21, 23, 29), (20, 21, 24, 26), (20, 21, 25, 26),
    (21, 22, 23, 25), (21, 22, 24, 26), (21, 22, 25, 26), (21, 22, 29, 30),
)

NODE_DIAGONALS = (
    (2, 7), (2, 13), (3, 7), (4, 8), (4, 14), (5, 9), (5, 13), (5, 23),
    (9, 23), (6, 10), (10, 15), (6, 11), (6, 16), (6, 25), (8, 17),
    (9, 17), (12, 18), (12, 24), (12, 19), (12, 25), (14, 16), (16, 23),
    (16, 26), (18, 19), (20, 24), (20, 25),
)

# pivot (non-basis) rays.
PIVOTS = {15, 22, 24, 26}
BASIS = [i for i in range(1, 31) if i not in PIVOTS]     # 26 rays
assert len(BASIS) == 26
ok("(q15,q22,q24,q26) is a maximal cone (unimodular pivot cone)",
   (15, 22, 24, 26) in MAXIMAL_CONES)

# restrictions of ambient divisors to the
# four exceptional surfaces (all unlisted divisors restrict to zero).
K6 = (-1, -1, -2, -2)
K7 = (-1, -1, -2)
SURFACES = {
    "E61": {1: (-1, 1, 0, 1), 2: (1, -1, 1, 0), 3: (1, 0, 0, 0),
            4: (0, 1, 0, 0), 5: (0, 0, 1, 0), 6: (0, 0, 0, 1), 27: K6},
    "E62": {7: (-1, 1, 0, 1), 8: (1, -1, 1, 0), 9: (1, 0, 0, 0),
            10: (0, 1, 0, 0), 11: (0, 0, 1, 0), 12: (0, 0, 0, 1), 28: K6},
    "E71": {1: (-1, 1, 1), 7: (1, -1, 0), 13: (1, 0, 0), 17: (0, 1, 0),
            21: (0, 0, 1), 29: K7},
    "E72": {2: (-1, 1, 0), 8: (1, -1, 1), 14: (1, 0, 0), 18: (0, 1, 0),
            30: K7},
}

# the ten local rows in Pic(E)^vee (the phi-rows of Table 2 of the paper).
DIV_ROWS = [
    ("E61", "s1_1", (0, 0, -1, 1)),
    ("E61", "s2_1", (-1, -1, 1, 0)),
    ("E61", "s3_1", (-1, 1, 1, -1)),
    ("E62", "s1_2", (0, 0, -1, 1)),
    ("E62", "s2_2", (-1, -1, 1, 0)),
    ("E62", "s3_2", (-1, 1, 1, -1)),
    ("E71", "alpha1", (1, -1, 0)),
    ("E71", "beta1", (1, 1, -1)),
    ("E72", "alpha2", (-2, 0, 1)),
    ("E72", "beta2", (-1, -1, 1)),
]

# ---------------------------------------------------------------- node rows
def cone_faces_contain(triple):
    s = set(triple)
    return any(s <= set(c) for c in MAXIMAL_CONES)

node_rows_30 = []
squares = []
for idx, (a, c) in enumerate(NODE_DIAGONALS, start=1):
    target = tuple(Q[a][k] + Q[c][k] for k in range(4))
    pairs = []
    for b in range(1, 27):
        if b in (a, c):
            continue
        for d in range(b + 1, 27):
            if d in (a, c):
                continue
            if tuple(Q[b][k] + Q[d][k] for k in range(4)) == target:
                pairs.append((b, d))
    good = [(b, d) for (b, d) in pairs
            if cone_faces_contain((a, b, c)) and cone_faces_contain((a, c, d))]
    assert len(good) == 1, (idx, a, c, pairs, good)
    b, d = good[0]
    squares.append((a, b, c, d))
    row = [F(0)] * 30
    row[a - 1] = row[c - 1] = F(-1)
    row[b - 1] = row[d - 1] = F(1)
    # circuit annihilates the lattice relations: sum row_i q_i = 0
    assert all(sum(row[i] * Q[i + 1][k] for i in range(30)) == 0
               for k in range(4))
    node_rows_30.append(row)
ok("all 26 unit squares re-derived from diagonals, cones, and the "
   "square identity q_a+q_c=q_b+q_d", len(node_rows_30) == 26)

# functionals on Pic in the dual basis of {D_j : j in BASIS}
def to_basis(row30):
    return [row30[j - 1] for j in BASIS]

node_rows = [to_basis(r) for r in node_rows_30]

ok("node subsystem has rank 18 (Paper 3, NODE-1)", rank(node_rows) == 18)
colo = [i for i in range(26)
        if rank(node_rows[:i] + node_rows[i + 1:]) == 17]
ok("node coloops are exactly N_9 and N_11 (Paper 3, NODE-2)",
   colo == [8, 10])

# ---------------------------------------------------------------- divisorial rows
div_rows = []
for surf, label, phi in DIV_ROWS:
    rest = SURFACES[surf]
    row = []
    for j in BASIS:
        col = rest.get(j)
        row.append(F(0) if col is None
                   else sum(F(p) * F(c) for p, c in zip(phi, col)))
    div_rows.append(row)
ok("ten divisorial rows assembled from printed restriction columns",
   len(div_rows) == 10 and rank(div_rows) == 10)

B = node_rows + div_rows          # 36 x 26
LABELS = [f"u{i}" for i in range(1, 27)] + [lb for _, lb, _ in DIV_ROWS]

ok("rank B = 21", rank(B) == 21)
ker_full = kernel([[B[r][c] for r in range(36)] for c in range(26)])  # B^T x=0
ok("dim ker(B^T) = 15", len(ker_full) == 15)

def forced_zero(ker_basis, coord):
    return all(v[coord] == 0 for v in ker_basis)

coloops = [LABELS[r] for r in range(36)
           if forced_zero(ker_full, r)]
ok("the single unrestricted row coloop is D_{7,1}:beta",
   coloops == ["beta1"])

# ---------------------------------------------------------------- profiles
IDX = {lb: i for i, lb in enumerate(LABELS)}
PROFILES = {
    "LL": ["s1_1", "s2_1", "s1_2", "s2_2", "alpha1", "alpha2"],
    "LP": ["s1_1", "s2_1", "s3_2", "alpha1", "alpha2"],
    "PL": ["s3_1", "s1_2", "s2_2", "alpha1", "alpha2"],
    "PP": ["s3_1", "s3_2", "alpha1", "alpha2"],
}
EXPECT_DIM = {"LL": 9, "LP": 10, "PL": 10, "PP": 12}
EXPECT_FORCED = {
    "LL": {"u9", "u11", "beta1", "beta2"},
    "LP": {"u9", "u11", "beta1", "beta2", "s3_1"},
    "PL": {"u9", "u11", "beta1", "beta2", "s3_2"},
    "PP": {"u9", "u11", "beta1", "beta2"},
}
# Label semantics (as in Section 5.3 of the paper): the profile letter
# is the component at each dP6 point, L = line {s1=s2=0}, P = plane
# {s3=0}; PROFILES lists the coordinates set to zero.  On LP the FREE line
# coordinate s3_1 of the first point is the one forced to vanish, matching
# the printed table of Section 5.3.

results = {}
for prof, prohibited in PROFILES.items():
    pro = [IDX[lb] for lb in prohibited]
    allowed = [i for i in range(36) if i not in pro]
    # restricted kernel: kernel of the restricted transpose matrix
    Bt_restricted = [[B[r][c] for r in allowed] for c in range(26)]
    kb = kernel(Bt_restricted)
    dim = len(kb)
    forced = {LABELS[allowed[k]] for k in range(len(allowed))
              if all(v[k] == 0 for v in kb)}
    results[prof] = (dim, forced)
    print(f"  {prof}: allowed {len(allowed)}, kernel dim {dim}, "
          f"forced zero: {sorted(forced)}")

for prof in PROFILES:
    dim, forced = results[prof]
    ok(f"{prof}: kernel dimension {EXPECT_DIM[prof]}",
       dim == EXPECT_DIM[prof])
    need = EXPECT_FORCED[prof]
    # the decisive facts: u9 and u11 forced; full expected list contained
    ok(f"{prof}: u9 and u11 forced to zero",
       "u9" in forced and "u11" in forced)
    ok(f"{prof}: full expected forced list holds",
       need <= forced)
    ok(f"{prof}: test not vacuous (u1 not forced)",
       "u1" not in forced)

# ---------------------------------------------------------------- robustness
# The enlarged profile {alpha1 = alpha2 = 0} (both dP6 points fully
# released): the forced zeros survive, so Theorem 1 is immune to every
# dP6-side question; and the forcing is knife-edge in the dP7 data:
# releasing either alpha destroys it.
def restricted_forced(prohibited_labels):
    pro = [IDX[lb] for lb in prohibited_labels]
    allowed = [i for i in range(36) if i not in pro]
    kb = kernel([[B[r][c] for r in allowed] for c in range(26)])
    forced = {LABELS[allowed[k]] for k in range(len(allowed))
              if all(v[k] == 0 for v in kb)}
    return len(kb), forced

dim13, forced13 = restricted_forced(["alpha1", "alpha2"])
ok("enlarged profile {alpha1 = alpha2 = 0}: kernel dimension 13",
   dim13 == 13)
ok("enlarged profile: u9, u11, beta1, beta2 all forced to zero "
   "(dP6 apparatus not needed for Theorem 1)",
   {"u9", "u11", "beta1", "beta2"} <= forced13)
dA, fA = restricted_forced(["alpha2"])
dB, fB = restricted_forced(["alpha1"])
ok("knife-edge: releasing either alpha loses the forcing of u9",
   "u9" not in fA and "u9" not in fB)

print()
print("Independently rebuilt from printed rays, diagonals, cones, "
      "restriction columns, and lambda-rows:")
print("  rank B = 21, dim ker(B^T) = 15, unique coloop D_{7,1}:beta;")
print("  on all four branch-restricted kernels the coordinates u9 and u11")
print("  vanish identically.  This is the decisive finite input of the")
print("  nonsmoothability theorem of the paper (Theorem 6.5).")
print()
print(f"ALL CHECKS PASSED ({PASS} assertions)")
