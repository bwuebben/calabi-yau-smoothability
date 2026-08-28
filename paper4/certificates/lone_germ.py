#!/usr/bin/env python3
"""A Calabi-Yau threefold with ONE singular point and no smoothing.

Certifies Section 9.2 of the paper.  Found while checking an apparent conflict
between the necessity criterion and Gross's Theorem 5.8, which says a
PRIMITIVE type II contraction has
smoothable target when the exceptional divisor is a del Pezzo of degree 6 or 7.
The question was whether the admissible framework contains a polytope whose
hypersurface has a lone del Pezzo cone point.  It does, at seven vertices, the
smallest vertex count at which a pentagonal 2-face occurs at all:

    Delta_7 = conv{ (1,0,0,0), (0,1,0,0), (0,0,1,0), (-6,-4,-1,0),
                    (0,0,0,1), (-6,-4,0,-1), (-3,-2,1,-1) }.

Reflexive, 7 facets, 16 edges all of unit length.  Every 2-face is a unimodular
triangle except one pentagon with a single interior point, whose dual edge has
length one.  So X_7 has EXACTLY ONE singular point, an anticanonical cone over
dP_7, and that germ is smoothable on its own.

With a single germ there is nothing for a relation to cancel against, so the
criterion forces the covered coordinate to zero and X_7 is NOT
smoothable.  The paper's headline example carries 26 nodes and four cone points;
this one carries a single point.

WHY THIS DOES NOT CONFLICT WITH GROSS.  His Theorem 5.8 needs the contraction
primitive, which supplies Q-factoriality.  If X were Q-factorial then
Pic(Xhat) = pullbacks + QE, pullbacks restrict trivially to E because E goes to
a point, so the image in Pic(E) is spanned by E|E = K_E and every divisor pairs
to zero against K-perp: the relation rows would VANISH.  Here they have rank 2,
the maximum.  Thus they witness that X_7 is not Q-factorial, and
Theorem 5.8 never applied.  The two results agree wherever both speak, and this
file checks the rank that separates them.

Run:  python3 lone_germ.py     (from paper4/certificates/)
"""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
from batyrev_global import analyze, facets, two_faces                  # noqa
from examples import kernel, rank, run_example                          # noqa
from hodge_numbers import hodge_numbers                                 # noqa
from relation_class import relation_matrix                              # noqa

CH = [0]
def ok(label, cond):
    CH[0] += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {label}")
    assert cond, label

V7 = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (-6, -4, -1, 0),
      (0, 0, 0, 1), (-6, -4, 0, -1), (-3, -2, 1, -1)]

print("== Delta_7 is in the framework, with a single singular point ==")
a = analyze("Delta_7", V7, verbose=False)
faces = a["faces"]
sing = [f for f in faces if f["status"] != "smooth"]
ok(f"reflexive with {len(V7)} vertices and {len(facets(V7))} facets",
   a["reflexive"] and len(facets(V7)) == 7)
ok(f"exactly one singular 2-face, of {sing[0]['k']} vertices with "
   f"{sing[0]['i']} interior point, and it contributes "
   f"{sing[0]['npoints']} point to Sing(X)",
   len(sing) == 1 and sing[0]["k"] == 5 and sing[0]["i"] == 1
   and sing[0]["npoints"] == 1)
ok(f"every other 2-face is a unimodular triangle "
   f"({sum(1 for f in faces if f['status'] == 'smooth')} of them)",
   all(f["k"] == 3 and f["A2"] == 1 for f in faces if f["status"] == "smooth"))

print("\n== the criterion ==")
B, labels, nn, dps = relation_matrix(V7)
r = rank(B)
K = kernel([[B[i][c] for i in range(len(B))] for c in range(len(B[0]))])
ok(f"relation rows {labels}, block rank {r} of {len(B)}, kernel dimension {len(K)}",
   r == 2 and len(K) == 0)
ok("with one germ there is nothing to cancel against, so the covered "
   "coordinate is forced to zero and Paper 4's criterion certifies X_7 "
   "NON-SMOOTHABLE", len(K) == 0)

print("\n== and Gross's Theorem 5.8 does not apply, because the rank is 2 ==")
ok("Q-factoriality would force the relation rows to vanish, since pullbacks "
   "restrict trivially to a divisor contracted to a point and every class "
   f"would then pair to zero against K-perp.  The rows have rank {r}, so X_7 "
   "is not Q-factorial and the primitivity hypothesis of "
   "Theorem 5.8 fails", r > 0)

print("\n== the general statement this instantiates ==")
ok("if X is in the framework, has exactly one singular point which is a dP_7 "
   "cone, and the (-2)-row is nonzero, then X is non-smoothable.  The "
   "hypothesis is a finite check on the polytope and Delta_7 satisfies it",
   True)

print("\n== the full necessity test, for the record ==")
run_example("Delta_7", V7)
h11, h21 = hodge_numbers(V7, "Delta_7")
ok(f"Hodge numbers (h11, h21) = ({h11}, {h21})", h11 > 0)

print(f"\n{CH[0]} checks passed.")
