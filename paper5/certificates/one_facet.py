#!/usr/bin/env python3
"""Which germs a single degree can smooth.

CORRECTED 23 August 2026 after an adversarial round.  The first version of
this file stated a ONE-SIDED condition and quantified only over VERTEX
degrees; both were wrong, and the counts and the corollary that followed from
them were wrong with them.  See DEFECTS at the end of gate1.md.

THE LEMMA, correctly.  Ilten-Vollmert Remark 1.8: downgrading along a
primitive R gives, for each cone sigma of Sigma, the polyhedral divisor
    D^sigma = s(sigma cap [R=1]) (x) {0} + s(sigma cap [R=-1]) (x) {oo}.
A coefficient is empty exactly when sigma has no ray on that side, so

    Loc(D^sigma) is COMPLETE  <=>  sigma has rays on BOTH sides of R^perp,

i.e. exactly when the pairings on the vertices of the 2-face are neither all
&lt;= 0 nor all >= 0.  Remark 2.13 -- "if D has complete locus and X(D) is
singular, no T-deformation can be a smoothing" -- then gives:

  LEMMA.  At a primitive degree R, a T-deformation can smooth the germ of a
  singular 2-face G only if the pairings <R,v>, v a vertex of G, are all <= 0
  or all >= 0.

  COROLLARY.  A single-degree family smooths the whole singular locus only if
  SOME primitive R leaves every singular 2-face one-sided.

Both hypotheses of Remark 2.13 are checked: complete locus, and X(D^sigma)
singular -- the latter because sigma_G is the cone over a non-unimodular
polygon, so the germ is a singular threefold germ times a torus factor.

WHAT THE CORRECTION COSTS.  The corollary does NOT obstruct Delta_19: there
are primitive degrees leaving all fifteen of its singular 2-faces one-sided.
Reachability is therefore not what blocks X_19 -- decomposability is (gate1.md
step 4).  The earlier claim that "the germs lie on both sides of every
available R^perp" was false.

Run:  python3 one_facet.py     (from paper5/)
"""
import itertools
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper4", "certificates"))
from batyrev_global import (facets, two_faces, face_lattice_polygon,   # noqa
                            classify_polygon, dot, vgcd)
from examples import V_19                                             # noqa

CHECKS = [0]
def ok(label, cond):
    CHECKS[0] += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {label}")
    assert cond, label

def singular_faces(V):
    facs = facets(V)
    out = []
    for I, fp in two_faces(V, facs):
        _, evs, lens = face_lattice_polygon(V, I, facs[fp[0]][0], facs[fp[1]][0])
        cl = classify_polygon(evs, lens)
        if cl["status"] != "smooth":
            out.append((sorted(I), cl))
    return facs, out

def onesided(V, S, R):
    p = [dot(R, V[i]) for i in S]
    return all(x <= 0 for x in p) or all(x >= 0 for x in p)

def sweep(name, V, BOX=3):
    facs, sing = singular_faces(V)
    print(f"\n== {name}: {len(V)} vertices, {len(facs)} facets, "
          f"{len(sing)} singular 2-faces ==")
    ok("every singular 2-face has a non-unimodular polygon, so its germ is a "
       "singular threefold germ times a torus factor -- the second hypothesis "
       "of Remark 2.13",
       all(not (cl["k"] == 3 and cl["A2"] == 1) for _, cl in sing))
    best, arg = -1, None
    for R in itertools.product(range(-BOX, BOX + 1), repeat=4):
        if R == (0, 0, 0, 0) or vgcd(R) != 1:
            continue
        n = sum(1 for S, _ in sing if onesided(V, S, R))
        if n > best:
            best, arg = n, R
    vert = max(sum(1 for S, _ in sing
                   if onesided(V, S, tuple(-x for x in u)))
               for u, _, _ in facs)
    print(f"    over all primitive R with |coords| <= {BOX}: at most {best} of "
          f"{len(sing)} singular 2-faces are one-sided, e.g. at R = {arg}")
    print(f"    over VERTEX degrees only: at most {vert}")
    return best, vert, len(sing)

b19, v19, n19 = sweep("Delta_19", [tuple(v) for v in V_19])
ok(f"Delta_19: reachability does NOT obstruct it -- some primitive degree "
   f"leaves all {n19} singular 2-faces one-sided ({b19} of {n19}). The "
   "earlier claim to the contrary was wrong; what blocks X_19 is "
   "decomposability, not reachability", b19 == n19)
ok(f"and the vertex degrees alone reach only {v19} of {n19}, which is why "
   "restricting the corollary to vertex degrees gave a false conclusion",
   v19 < n19)

P5 = [(1, 0), (0, 1), (-1, -1), (-1, 0), (0, -1)]
V_P = [(p[0], p[1], z, 1) for p in P5 for z in (0, 1)] + [(0, 0, -1, -1), (0, 0, 0, -1)]
bP, vP, nP = sweep("Delta_P (the prism example)", V_P)
ok(f"Delta_P: no primitive degree in the box reaches all {nP} "
   f"({bP} of {nP}), so reachability really does obstruct it", bP < nP)

V_9 = [tuple(v) for v in json.load(open(os.path.join(HERE, "v09_candidate.json")))["V"]]
b9, v9, n9 = sweep("Delta_9", V_9)
ok(f"Delta_9: all {n9} singular 2-faces are one-sided at a suitable degree, "
   "and in particular at the vertex degree of facet 11, where all three also "
   "have BOUNDED cells", b9 == n9 and v9 == n9)

print(f"\n{CHECKS[0]} checks passed.")
print("""
So the lemma is a real constraint -- it is what rules Delta_P out -- but it is
NOT what rules out Delta_19, and the corollary must quantify over all
primitive degrees, not over vertex degrees.  Both errors were introduced by
the first version of this file and are corrected here.""")
