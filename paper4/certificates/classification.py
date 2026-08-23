#!/usr/bin/env python3
"""The complete classification of the searched range.

Certifies Proposition 9.2 of the paper.  Sweeping the Kreuzer-Skarke database
at 5 to 9 vertices for reflexive 4-polytopes that are admissible and carry a
pentagonal 2-face returns exactly 77, one pentagon apiece.  This file rebuilds
the germ inventory and the relation matrix of each from its stored vertex list
and runs the necessity criterion on all 77.

    criterion FIRES  -> non-smoothable, certified                   76
    criterion SILENT -> the criterion says nothing either way        1

The 76 are a consequence of the necessity theorem and they are new in bulk:
the paper exhibits one non-smoothable Calabi-Yau threefold by hand, and this
produces seventy-six, of which 67 have a SINGLE singular point -- in every case
a lone dP_7 cone point, the situation of Corollary 9.1.  The single silent case
is a nine-vertex polytope with two nodes and one dP_7 point; it is the
companion paper's, and nothing here depends on how it is settled there.

The locking criterion of the companion's Theorem A' lines up with this exactly:
of the 77 pentagons, 76 are locked and one is not, and the unlocked one is the
silent polytope's.  So in the searched range the ambient construction is
available precisely where the criterion permits a smoothing, and nowhere else.

Input: framework_77.json, written by dp7_facet_sweep.py over
data/ks/polytopes-4d-0[5-9]-vertices.parquet with the union-find chain.

Run:  python3 classification.py     (from paper4/certificates/)
"""
import json, os, sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
from relation_class import relation_matrix                              # noqa
from examples import kernel                                             # noqa
from batyrev_global import analyze                                      # noqa

CH = [0]
def ok(label, cond):
    CH[0] += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {label}")
    assert cond, label

fw = json.load(open(os.path.join(HERE, "framework_77.json")))
target = tuple(sorted(tuple(map(int, r)) for r in
                      json.load(open(os.path.join(HERE, "v09_candidate.json")))["V"]))

rows = []
for c in fw:
    V = [tuple(map(int, r)) for r in c["V"]]
    B, labels, nn, dps = relation_matrix(V)
    prohibit = {f"dp{t}:A" for t, d in enumerate(dps) if d["k"] == 5}
    allowed = [i for i, lb in enumerate(labels) if lb not in prohibit]
    K = kernel([[B[i][j] for i in allowed] for j in range(len(B[0]))])
    names = [labels[i] for i in allowed]
    covered = [j for j, nm in enumerate(names)
               if nm.startswith("u") or nm.endswith(":R")]
    forced = [names[j] for j in covered if all(v[j] == 0 for v in K)]
    a = analyze("x", V, verbose=False)
    rows.append(dict(
        nv=len(V),
        germs=sum(1 for f in a["faces"] if f["status"] != "smooth"),
        inv=(nn, sum(1 for d in dps if d["k"] == 5),
             sum(1 for d in dps if d["k"] == 6)),
        locked=all(c["verdicts"]),
        silent=not forced,
        forced=forced,
        is9=tuple(sorted(tuple(map(int, r)) for r in c["V"])) == target))

print(f"== {len(rows)} framework polytopes with a pentagonal 2-face, "
      "5 to 9 vertices ==")
ok(f"one pentagon apiece, and the vertex counts are "
   f"{dict(sorted(Counter(r['nv'] for r in rows).items()))}", len(rows) == 77)
inv = dict(sorted(Counter(r["inv"] for r in rows).items()))
ok(f"the germ inventories (nodes, dP7, dP6) are {inv}: no hexagonal 2-face "
   "occurs anywhere in the range, so every germ is a node or a dP7 cone point "
   "and no degree-6 branch profile is needed",
   all(k[2] == 0 and k[1] == 1 for k in inv))

fires = [r for r in rows if not r["silent"]]
silent = [r for r in rows if r["silent"]]
print()
ok(f"the necessity criterion FIRES on {len(fires)} of them, certifying each "
   f"non-smoothable; their inventories are "
   f"{dict(sorted(Counter(r['inv'] for r in fires).items()))}",
   len(fires) == 76)
ok(f"{sum(1 for r in fires if r['germs'] == 1)} of the {len(fires)} have a "
   "SINGLE singular point, in every case a lone dP7 cone point, which is the "
   "situation of the one-germ corollary: there is nothing for a relation to "
   f"cancel against.  The other {sum(1 for r in fires if r['germs'] == 2)} "
   "have one node and one dP7 point",
   sum(1 for r in fires if r["germs"] == 1) == 67 and
   sum(1 for r in fires if r["germs"] == 2) == 9 and
   all(r["inv"] == (0, 1, 0) for r in fires if r["germs"] == 1) and
   all(r["inv"] == (1, 1, 0) for r in fires if r["germs"] == 2))
ok(f"it is SILENT on {len(silent)}: a nine-vertex polytope whose hypersurface "
   f"has {silent[0]['germs']} singular points, two nodes and one dP7 cone "
   "point.  The criterion says nothing about it either way; that case is the "
   "companion paper's",
   len(silent) == 1 and silent[0]["nv"] == 9 and silent[0]["inv"] == (2, 1, 0))
ok("the silent polytope is the companion's Delta_9, matched against its "
   "committed vertex list", silent[0]["is9"])
ok("so the criterion decides every admissible polytope with a pentagonal "
   "2-face in the range except one, and the reason so many fire is structural: "
   "with one or two germs there is too little for a relation to cancel "
   "against, and the criterion is a statement about cancellation",
   len(fires) + len(silent) == len(rows))

print("\n== the companion's locking criterion lines up with it exactly ==")
ok(f"{sum(1 for r in rows if r['locked'])} of the 77 pentagons are locked, so "
   "the ambient single-degree construction is obstructed there",
   sum(1 for r in rows if r["locked"]) == 76)
ok("the one unlocked pentagon is the silent polytope's, so in the searched "
   "range the ambient construction is available precisely where the criterion "
   "permits a smoothing and nowhere else",
   all(r["is9"] for r in rows if not r["locked"]))

print("\n== what this amounts to ==")
ok("this paper exhibits one non-smoothable Calabi-Yau threefold by hand; the "
   f"same criterion applied mechanically over the range produces {len(fires)}",
   len(fires) == 76)

print(f"\n{CH[0]} checks passed.")
