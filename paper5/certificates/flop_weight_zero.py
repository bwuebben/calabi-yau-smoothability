#!/usr/bin/env python3
"""Weight-zero Hodge cohomology across the degree-seven flop (Lemma 4.4).

Near a degree-seven point the global resolution Y is the flop of
A = Tot(K_E) along c = h - e_1 - e_2.  Write A' for the flopped
neighbourhood and O = A minus c = A' minus c'.  The proof of Lemma 4.4 uses
that the local tangent space T^1, which after contraction with the toric
volume form lies in torus weight zero, injects into H^2_{E minus c}(O, Omega^2).
The kernel is the image of H^1(O, Omega^2), so it suffices that

    H^1(O, Omega^2)_0 = 0.

On the affine chart of a smooth cone sigma, the torus-invariant regular
p-forms are wedge^p of (sigma^perp tensor Q); Cech cohomology on the
maximal cones computes H^q(-, Omega^p)_0 for every fan below.  The program
also checks H^1(U, Omega^2)_0 = 2 = dim T^1 on the punctured cone, and the
Hodge numbers of P^2, P^1 x P^1 and P^3 as controls.

Run:  sage -python flop_weight_zero.py      (or python3 with sympy installed)
"""
import itertools
from sympy import Matrix, eye, zeros

CHECKS = [0]


def ok(label, cond):
    CHECKS[0] += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {label}")
    assert cond, label


def perp_basis(rays, n):
    if not rays:
        return [list(r) for r in eye(n).tolist()]
    return [list(v) for v in Matrix(rays).nullspace()]


def wedge_space(vecs, p, n):
    idx = list(itertools.combinations(range(n), p))
    dim = len(idx)
    if len(vecs) < p:
        return zeros(dim, 0)
    cols = []
    for comb in itertools.combinations(range(len(vecs)), p):
        M = Matrix([vecs[c] for c in comb])
        cols.append([M[:, list(I)].det() for I in idx])
    basis = Matrix(cols).T.columnspace()
    return Matrix.hstack(*basis) if basis else zeros(dim, 0)


def weight_zero_cohomology(maxcones, p, n=3):
    """dims of H^q(X_Sigma, Omega^p)_0, q = 0..m-1, by Cech on maximal cones."""
    m = len(maxcones)
    space = {}
    for k in range(1, m + 1):
        for S in itertools.combinations(range(m), k):
            common = set(maxcones[S[0]])
            for s in S[1:]:
                common &= set(maxcones[s])
            space[S] = wedge_space(perp_basis([list(r) for r in sorted(common)], n), p, n)

    def cdim(k):
        return sum(space[S].shape[1] for S in itertools.combinations(range(m), k + 1))

    def dmat(k):
        src = list(itertools.combinations(range(m), k + 1))
        tgt = list(itertools.combinations(range(m), k + 2))
        soff, o = {}, 0
        for S in src:
            soff[S] = o; o += space[S].shape[1]
        toff, o2 = {}, 0
        for T in tgt:
            toff[T] = o2; o2 += space[T].shape[1]
        D = zeros(o2, o)
        for T in tgt:
            BT = space[T]
            if BT.shape[1] == 0:
                continue
            for j in range(len(T)):
                S = T[:j] + T[j + 1:]
                BS = space[S]
                if BS.shape[1] == 0:
                    continue
                sol = (BT.T * BT).inv() * BT.T * BS
                assert BT * sol == BS
                for a in range(BT.shape[1]):
                    for b in range(BS.shape[1]):
                        D[toff[T] + a, soff[S] + b] += (-1) ** j * sol[a, b]
        return D

    out = []
    for k in range(m):
        dk = dmat(k) if k < m - 1 else zeros(0, cdim(k))
        dkm1 = dmat(k - 1) if k > 0 else zeros(cdim(0), 0)
        rk = dk.rank() if dk.shape[0] * dk.shape[1] else 0
        rkm = dkm1.rank() if dkm1.shape[0] * dkm1.shape[1] else 0
        out.append(cdim(k) - rk - rkm)
    return out


# The degree-seven cone: the pentagon of Proposition 3.7 at height one.
P = [(0, 0, 1), (1, 0, 1), (2, 1, 1), (1, 2, 1), (0, 1, 1)]
c = (1, 1, 1)
star = [[c, P[i], P[(i + 1) % 5]] for i in range(5)]                     # A = Tot(K_E)
flop = [[P[0], P[1], P[4]], [c, P[1], P[4]],
        [c, P[1], P[2]], [c, P[2], P[3]], [c, P[3], P[4]]]             # A'
punct = [[P[i], P[(i + 1) % 5]] for i in range(5)]                     # U_p
# O: remove the cone of c = (c, P0) from the star fan (equivalently c' = (P1, P4)
# from the flopped fan); the remaining maximal cones:
common_open = [[c, P[1], P[2]], [c, P[2], P[3]], [c, P[3], P[4]],
               [P[0], P[1]], [P[4], P[0]]]

def triangulates_pentagon(cones):
    """Unimodularity, total area and boundary of the marked subdivision."""
    triangles = {frozenset(cone) for cone in cones}
    if len(triangles) != 5 or any(len(t) != 3 for t in triangles):
        return False
    determinants = [abs(Matrix(sorted(t)).det()) for t in triangles]
    polygon_area = abs(sum(P[i][0] * P[(i + 1) % 5][1]
                           - P[i][1] * P[(i + 1) % 5][0]
                           for i in range(5)))
    edges = {}
    for t in triangles:
        for edge in itertools.combinations(t, 2):
            edge = frozenset(edge)
            edges[edge] = edges.get(edge, 0) + 1
    boundary = {frozenset((P[i], P[(i + 1) % 5])) for i in range(5)}
    return (all(d == 1 for d in determinants)
            and sum(determinants) == polygon_area
            and {edge for edge, count in edges.items() if count == 1} == boundary
            and all(count in (1, 2) for count in edges.values()))


def complement_fan(cones, removed_face):
    """Maximal cones remaining after removing an invariant curve."""
    remaining = set()
    for cone in cones:
        for size in range(len(cone) + 1):
            for face in itertools.combinations(cone, size):
                face = frozenset(face)
                if not removed_face <= face:
                    remaining.add(face)
    return {face for face in remaining
            if not any(face < other for other in remaining)}


print("== the marked fan subdivisions ==")
ok("flop circuit: c + P0 = P1 + P4",
   all(c[i] + P[0][i] == P[1][i] + P[4][i] for i in range(3)))
ok("star fan triangulates the entire pentagon", triangulates_pentagon(star))
ok("flopped fan triangulates the entire pentagon", triangulates_pentagon(flop))
marked_open = {frozenset(cone) for cone in common_open}
ok("common open is the star fan minus the flopping curve",
   complement_fan(star, frozenset((c, P[0]))) == marked_open)
ok("common open is the flopped fan minus the flopped curve",
   complement_fan(flop, frozenset((P[1], P[4]))) == marked_open)

print("== controls ==")
p2 = [[(1, 0, 0), (0, 1, 0)], [(0, 1, 0), (-1, -1, 0)], [(-1, -1, 0), (1, 0, 0)]]
p2 = [[r[:2] for r in cone] for cone in p2]
ok("P^2: h^{1,1} = 1", weight_zero_cohomology(p2, 1, 2)[1] == 1)
p1p1 = [[(1, 0), (0, 1)], [(0, 1), (-1, 0)], [(-1, 0), (0, -1)], [(0, -1), (1, 0)]]
ok("P^1 x P^1: h^{1,1} = 2", weight_zero_cohomology(p1p1, 1, 2)[1] == 2)
rays3 = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, -1, -1)]
p3 = [list(x) for x in itertools.combinations(rays3, 3)]
ok("P^3: h^{2,2} = 1", weight_zero_cohomology(p3, 2, 3)[2] == 1)

print("== the degree-seven cone ==")
ok("punctured cone: H^1(U, Omega^2)_0 has dimension 2 = dim T^1",
   weight_zero_cohomology(punct, 2)[1] == 2)
ok("star resolution: H^1(A, Omega^2)_0 = 0", weight_zero_cohomology(star, 2)[1] == 0)
ok("flopped resolution: H^1(A', Omega^2)_0 = 0", weight_zero_cohomology(flop, 2)[1] == 0)
ok("common open set: H^1(O, Omega^2)_0 = 0", weight_zero_cohomology(common_open, 2)[1] == 0)
print(f"\n{CHECKS[0]} checks passed.")
