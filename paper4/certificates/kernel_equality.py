#!/usr/bin/env python3
"""The certified kernel IS the topological kernel, on the whole framework.

Certifies Proposition 8.3 of the paper.

The necessity theorem is stated with the topological kernel

    K(X) = ker( sum_p H_2(L_p; Q)  ->  H_2(Xhat; Q) ),

but Section 5.3 computes ker(B^T), the kernel of the matrix of pairings
with AMBIENT divisor classes, and uses only the containment
"ker(B^T) contains K(X)": a larger certified kernel makes every forced-zero
statement stronger, so non-smoothability is safe either way.

For the CONVERSE direction the containment runs the wrong way.  "No coordinate
is forced to zero on ker(B^T)" does not imply the same for K(X), so a
sufficiency conjecture phrased with the computable test would not be well posed.
This file removes the difficulty for the whole admissible framework at once.

    LEMMA.  Let Delta be in the admissible framework.  Then the correction term
    in Batyrev's formula for h^{1,1},

            sum over 2-faces F of Delta of  l*(F) * ( l(F*) - 1 ),

    vanishes.  Hence h^{1,1}(Xhat) is purely toric, the ambient divisor classes
    span H^{1,1}(Xhat; Q), and ker(B^T) = K(X).

    PROOF.  A singular 2-face has dual edge of lattice length one, by the
    framework condition, so its factor (l(F*) - 1) is zero.  A non-singular
    2-face is a unimodular triangle, so its factor l*(F) is zero.  Every term
    dies for one reason or the other.

The vanishing is NOT the same as the one that empties Mavlyutov's mechanism,
though the two are neighbours: his count is the mirror correction, summed over
2-faces of the DUAL, and it vanishes because every edge of Delta has lattice
length one.  Both hold on the framework, for related but distinct reasons, and
both are checked here.

Runs on X-circ and its mirror partner, on Delta_20, Delta_19 and Delta_7, and
over all 77 census polytopes of framework_77.json.

Run:  python3 kernel_equality.py     (from paper4/certificates/)
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
from hodge_numbers import (hodge_numbers, facets, two_faces, lattice_points,
                           face_interior_count, vgcd, vsub)
from examples import V_19, V_20, V_F1, polar

CH = [0]
def ok(label, cond):
    CH[0] += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {label}")
    assert cond, label


def corrections(V):
    facs = facets(V)
    Vd = [u for u, _, _ in facs]
    facsd = facets(Vd)
    out = []
    for (VV, ff) in ((V, facs), (Vd, facsd)):
        pts = lattice_points(VV, ff)
        s, bad = 0, []
        for I, fp in two_faces(VV, ff):
            tight = frozenset(i for i in range(len(ff)) if I <= ff[i][2])
            lint = face_interior_count(pts, ff, tight)
            dl = vgcd(vsub(ff[fp[0]][0], ff[fp[1]][0]))
            if lint * (dl - 1) != 0:
                bad.append((sorted(I), lint, dl))
            s += lint * (dl - 1)
        out.append((s, bad))
    return out


V9 = [tuple(map(int, r)) for r in
      json.load(open(os.path.join(HERE, "v09_candidate.json")))["V"]]
V7 = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (-6, -4, -1, 0),
      (0, 0, 0, 1), (-6, -4, 0, -1), (-3, -2, 1, -1)]

# Paper 4's own examples first: X-circ is the hypersurface of the polar of
# Delta_F1, X_Delta the hypersurface of Delta_F1 itself.  X_Delta carries an
# F_1 quadrilateral 2-face, so it is not in the admissible framework, but the
# hypothesis of the lemma -- every 2-face is a unimodular triangle or has dual
# edge of lattice length one -- holds for it too, and the conclusion follows.
CASES = (("X-circ (polar Delta_F1)", polar(V_F1)),
         ("X_Delta (Delta_F1)", V_F1),
         ("Delta_20", V_20),
         ("Delta_19", V_19),
         ("Delta_9", V9),
         ("Delta_7", V7))

print("== the two Batyrev correction terms, on the framework examples ==")
for nm, V in CASES:
    (s11, bad11), (s21, bad21) = corrections(V)
    h11, h21 = hodge_numbers(V)
    ok(f"{nm}: the h^11 correction is {s11} (witnesses {bad11}), so h^11 = "
       f"{h11} is purely toric and the ambient divisor classes span "
       "H^{1,1}(Xhat)", s11 == 0)
    ok(f"{nm}: the h^21 correction is {s21} (witnesses {bad21}), which is "
       f"Mavlyutov's count of non-polynomial deformations; h^21 = {h21}",
       s21 == 0)

print("\n== why each term vanishes, 2-face by 2-face ==")
for nm, V in CASES:
    facs = facets(V)
    pts = lattice_points(V, facs)
    r = {"singular, dual edge length 1": 0, "unimodular triangle": 0, "OTHER": 0}
    for I, fp in two_faces(V, facs):
        tight = frozenset(i for i in range(len(facs)) if I <= facs[i][2])
        lint = face_interior_count(pts, facs, tight)
        dl = vgcd(vsub(facs[fp[0]][0], facs[fp[1]][0]))
        if len(I) > 3 and dl == 1:
            r["singular, dual edge length 1"] += 1
        elif lint == 0:
            r["unimodular triangle"] += 1
        else:
            r["OTHER"] += 1
    print(f"    {nm}: " + ", ".join(f"{k} = {v}" for k, v in r.items()))
    ok(f"{nm}: every 2-face dies for one of the two stated reasons, with none "
       "left over", r["OTHER"] == 0)

print("\n== the other vanishing: every edge of Delta is unit ==")
for nm, V in CASES:
    facs = facets(V)
    # an edge of Delta is a pair of vertices whose common facets meet in
    # exactly that pair; taking all pairs inside a 2-face would pick up
    # diagonals, which are not edges and may well be long
    ln, nedges = [1], 0
    for i in range(len(V)):
        for j in range(i + 1, len(V)):
            common = [f for f in facs if i in f[2] and j in f[2]]
            if not common:
                continue
            inter = set(common[0][2])
            for f in common[1:]:
                inter &= set(f[2])
            if inter == {i, j}:
                nedges += 1
                ln.append(vgcd(vsub(V[i], V[j])))
    ok(f"{nm}: {nedges} edges, the longest of lattice length {max(ln)}, which "
       "is what "
       "empties Mavlyutov's mechanism, and is a different statement from the "
       "one above", max(ln) == 1)

print("\n== and over the whole census: all 77 framework polytopes ==")
fw = json.load(open(os.path.join(HERE, "framework_77.json")))
bad11 = bad21 = 0
for c in fw:
    V = [tuple(map(int, r)) for r in c["V"]]
    (s11, _), (s21, _) = corrections(V)
    bad11 += (s11 != 0)
    bad21 += (s21 != 0)
ok(f"the h^11 correction vanishes on all {len(fw)} of them "
   f"({bad11} exceptions), so the conclusion is a property of the framework "
   "and not of the examples", bad11 == 0)
ok(f"the h^21 correction vanishes on all {len(fw)} as well ({bad21} "
   "exceptions)", bad21 == 0)

print(f"\n{CH[0]} checks passed.")
print("""
CONCLUSION.  On the admissible framework the ambient Picard group detects all
of H^{1,1}(Xhat), so Paper 4's certified kernel is the topological kernel and
not merely an upper bound for it.  That is what makes a sufficiency conjecture
phrased in terms of the computable test well posed; with only the containment
it would not be.""")
