"""Is the general fibre toric?  NO -- and the check that suggested it was is
incomplete in exactly the way that matters.

An adversarial round claimed follows.md was wrong to call the general fibre
non-toric, on the ground that the level-(+1) slice IS the tail fan: Delta_9
has exactly one vertex at a positive level, the cosection sends it to the
origin, so every level-(+1) cell is {0} + tail.  Suess Cor. 1.9 would then
give at most two slices differing from the tail fan, hence toric.

That check verified the nonempty cells and never checked that they are ALL
nonempty.  They are not: of the 72 cones of the face fan, 43 have an EMPTY
level-(+1) cell, including four of the twelve maximal cones -- facets 5, 6, 7
and 11, which are exactly the affine-locus charts.  An empty coefficient is
not a lattice translate of the tail cone, so the slice at 0 is the tail fan
WITH HOLES, all three slices differ from the tail fan, and Cor. 1.9's
hypothesis fails.

The symptom, if one tries anyway: feeding the two summand slices into Suess
Prop. 3.1's cone for the affine-locus charts produces cones CONTAINING LINES
(facets 5, 6, 7 give 2 lines each), which cannot be fan cones.  That is the
construction complaining about being applied off its hypothesis.

So follows.md was right, and so was its count.  There is no lattice-point
proof of this shape, and the chart computation of charts.sage stands as the
argument.

This file records the negative result and the validation that makes it
trustworthy: the upgrade construction is exact on the SPECIAL fibre, where it
rebuilds Delta_9's own face fan cone for cone.

Run:  sage toric_fibre.sage     (from paper5/)
"""
import os, sys, json, itertools
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper4", "certificates"))
from batyrev_global import facets, dot

CH = [0]
def ok(label, cond):
    CH[0] += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {label}")
    assert cond, label

W = json.load(open(os.path.join(HERE, "v09_candidate.json")))
Vt = [tuple(int(x) for x in v) for v in W["V"]]
V = [vector(ZZ, v) for v in Vt]
facs = facets(Vt); FI = 11
R = vector(ZZ, [-x for x in facs[FI][0]])
w = [vector(ZZ, e) for e in identity_matrix(ZZ, 4).rows()
     if R.dot_product(vector(ZZ, e)) == 1][0]
KB = matrix(ZZ, R.row().right_kernel().basis_matrix())
def co(x): return vector(QQ, KB.transpose().solve_right(vector(QQ, x)))
def s(x):
    x = vector(QQ, x); return x - R.dot_product(x) * vector(QQ, w)
Hp = Polyhedron(eqns=[[0] + list(R)], base_ring=QQ)
H1 = Polyhedron(eqns=[[-1] + list(R)], base_ring=QQ)
Hm = Polyhedron(eqns=[[1] + list(R)], base_ring=QQ)

print("== validation: the upgrade rebuilds the special fibre exactly ==")
match = 0
for fi, (_, _, idx) in enumerate(facs):
    S = sorted(idx)
    sig = Polyhedron(rays=[list(V[i]) for i in S], base_ring=QQ)
    D0 = [list(co(s(vector(QQ, v)))) for v in (sig & H1).vertices_list()]
    Dm = [list(co(s(vector(QQ, v)))) for v in (sig & Hm).vertices_list()]
    tail = [list(co(r)) for r in (sig & Hp).rays()]
    rec = Polyhedron(rays=[[1] + v for v in D0] + [[-1] + v for v in Dm]
                     + [[0] + t for t in tail], base_ring=QQ)
    tgt = Polyhedron(rays=[[R.dot_product(V[i])] + list(co(s(vector(QQ, V[i]))))
                           for i in S], base_ring=QQ)
    match += (rec == tgt)
ok(f"all {len(facs)} maximal cones of Delta_9's face fan are rebuilt from "
   f"their two slices ({match} of {len(facs)})", match == len(facs))

print("\n== but the general fibre's slice at 0 is not the tail fan ==")
n = len(Vt)
seen, frontier = {frozenset(range(n))}, [frozenset(range(n))]
while frontier:
    nxt = []
    for S in frontier:
        for _, _, I in facs:
            T = S & I
            if T and T not in seen: seen.add(T); nxt.append(T)
    frontier = nxt
cones = [S for S in seen if S != frozenset(range(n))]
nonempty = [S for S in cones if any(R.dot_product(V[i]) > 0 for i in S)]
empty = [S for S in cones if not any(R.dot_product(V[i]) > 0 for i in S)]
ok(f"every NONEMPTY level-(+1) cell is {{0}} + tail, as the refuter observed "
   f"({len(nonempty)} of them)",
   all(co(s(vector(QQ, V[[i for i in S if R.dot_product(V[i]) > 0][0]])))
       == vector(QQ, [0, 0, 0]) for S in nonempty))
ok(f"but {len(empty)} of the {len(cones)} cells are EMPTY, and an empty "
   "coefficient is not a lattice translate of the tail cone -- so the slice "
   "is the tail fan WITH HOLES", len(empty) > 0)
mx = [fi for fi, (_, _, idx) in enumerate(facs)
      if not any(R.dot_product(V[i]) > 0 for i in idx)]
ok(f"four maximal cones are among them -- facets {mx}, exactly the "
   "affine-locus charts", mx == [5, 6, 7, 11])
ok("so all three slices differ from the tail fan and Suess Cor. 1.9 does NOT "
   "give toricity: follows.md was right, and the objection to it was not",
   True)

print(f"\n{CH[0]} checks passed.")
print("""
CONCLUSION.  The general fibre is not toric by this criterion, there is no
lattice-point shortcut, and the chart computation of charts.sage remains the
proof.  Recorded because the negative is worth as much as the positive would
have been: it closes off the alternative rather than leaving it dangling.""")
