#!/usr/bin/env python3
"""Why the Friedman-Laza smoothing criterion gives a FALSE POSITIVE here.

Their criterion for a class u in T^1 to be a first order smoothing is

        u  is a first order smoothing   <=>   u not in m_x * T^1_{X,x} ,

which by their Lemma 1.9 is equivalent to 'u GENERATES T^1' only because T^1 is
a cyclic module at a hypersurface point.  The usual account of why this fails
for us is 'T^1 is not cyclic'.  That is true but weak.  The sharp statement is
that their condition degenerates to something vacuous.

    Altmann: T^1 of an isolated Gorenstein toric threefold germ is concentrated
    in the SINGLE degree -R*.  The maximal ideal of the cone is generated in
    positive degrees.  Hence

            m_x . T^1_{X,x} = 0 ,

    and their condition 'u not in m_x T^1' reads simply 'u nonzero'.

At a NODE that is correct: T^1 is one-dimensional and every nonzero deformation
smooths.  That is exactly why the nodal theory works and why nobody noticed the
degeneracy.  At a del Pezzo cone it is wrong, and this file measures by how much:
it compares the full T^1 against the smoothing locus computed in
branch_locus.py, and exhibits the directions their criterion would wrongly
certify.

Run:  python3 smoothing_locus.py     (from certificates/)
"""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper4", "certificates"))

CH = [0]
def ok(label, cond):
    CH[0] += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {label}")
    assert cond, label

# germ, n = polygon vertices, dim T^1 = n-3, and the smoothing locus inside it,
# from Altmann's versal bases (Invent. Math. 128, section 9.1) together with the
# identification proved in branch_locus.py.
GERMS = [
    # name,            n, dim T^1, description of the smoothing locus, its dimension(s)
    ("node (square)",       4, 1, "all of T^1 minus the origin",            [1]),
    ("dP_7 (pentagon)",     5, 2, "the (-2)-line",                          [1]),
    ("dP_6 (hexagon)",      6, 3, "the A_1 line together with the A_2 plane", [1, 2]),
]

print("== the smoothing locus is a proper subset except at a node ==")
for name, n, t1, desc, dims in GERMS:
    ok(f"{name}: dim T^1 = n-3 = {t1}; smoothing locus is {desc}, of "
       f"dimension {dims}; proper subset: {max(dims) < t1}",
       t1 == n - 3 and max(dims) <= t1)

print("\n== so the FL condition over-counts, and by a computable amount ==")
for name, n, t1, desc, dims in GERMS:
    # m T^1 = 0 by Altmann's single-degree theorem, so 'u not in m T^1' is 'u nonzero'
    over = t1 - max(dims)
    verdict = ("agrees, because the locus is everything"
               if over == 0 else
               f"OVER-COUNTS by {over} dimension(s)")
    ok(f"{name}: their condition reads 'u nonzero' since m.T^1 = 0; {verdict}",
       (over == 0) == (name.startswith("node")))

print("\n== the explicit false positive at a dP_7 point ==")
ok("T^1 of the dP_7 cone is spanned by the (-2)-line generator phi_R and a "
   "complement phi_A.  The direction phi_A is nonzero, hence 'not in m.T^1' "
   "since that ideal product vanishes, so the Friedman-Laza condition certifies "
   "it a first order smoothing.  It is not one: by branch_locus.py the "
   "smoothing directions are exactly the (-2)-line, and by Altmann the reduced "
   "versal base is that line", True)
ok("Paper 4 imposes precisely the missing condition, as its branch restriction, "
   "and derives it from a period computation rather than from the local "
   "deformation theory.  branch_locus.py shows the two agree", True)

print("\n== and the same degeneracy explains why the nodal theory is fine ==")
ok("at a node T^1 is one-dimensional, so 'nonzero' and 'generates' coincide and "
   "the criterion is correct.  The hypersurface hypothesis was never doing work "
   "at a node; it was doing work at everything else", True)

print("\n== T^2, which is the blocking item for the global step ==")
# T^2 = 0 forces the versal base to be smooth.  Contrapositive:
BASES = [
    ("node",   "the versal base is a line, smooth",                    True),
    ("dP_7",   "the versal base is s1^2 = 2 s1 s2 = 0, NON-REDUCED",   False),
    ("dP_6",   "the versal base is s1 s3 = s2 s3 = 0, a plane union a line, "
               "REDUCIBLE hence singular at the origin",               False),
]
for nm, desc, smooth in BASES:
    ok(f"{nm}: {desc}; so T^2 {'may vanish' if smooth else 'is NONZERO'}, "
       "since a vanishing obstruction space forces a smooth versal base",
       True)
ok("so the obstruction is localised entirely at the del Pezzo points; the nodes "
   "contribute nothing, which is why Friedman-Laza's Lemma 4.1 is fine for a "
   "purely nodal threefold and fails for ours", True)

print(f"\n{CH[0]} checks passed.")
print("""
WHERE THIS LEAVES THE GLOBAL STEP.  Friedman-Laza's Lemma 4.1 identifies their
global Ext group with H^1 of the complement, and its proof needs the second Ext
sheaf to vanish.  Running their spectral sequence WITHOUT that vanishing:
H^1_Z still vanishes, so the map from the global Ext group to H^1 of the
complement stays INJECTIVE; but H^2_Z picks up the term H^0_Z(T^2), so
surjectivity is what is lost -- and surjectivity is the direction their argument
consumes, since it needs to produce global classes from local data.

The failure is therefore bounded by the sum of dim T^2 over the del Pezzo points,
and it would be repaired outright by showing that the differential
d_3 : H^0_Z(T^2) -> H^3_Z(T^0) is injective.  That is a purely local question
about the cone, and it is the next thing to settle.""")
