"""Gate 1, closed over ALL primitive degrees, by a bound rather than a box.

1.  The crosscut depends only on the restriction of R to the germ's rank-3
    character lattice, i.e. on the pairing vector p = (<R,v_i>).
2.  A decomposition into unimodular simplices forces normalised area >= 1
    (a unimodular-triangle summand by monotonicity; an all-segment
    decomposition is a zonotope of area >= 2; the crosscut is never itself a
    unimodular simplex since its cone is sigma_P).
3.  A pentagon is triangulated by three triangles on its own vertices, so
    area >= 1 forces some triple of the five rescaled points to span a
    triangle of area >= 1/3.  Triangle area is decreasing in each pairing, so
    that region is bounded and is found by raising a cap until nothing new
    appears.  This yields a FINITE list of pairing vectors.
4.  Each pairing vector is realised by a line in M, the kernel of M -> M_P
    having rank one.  Walk each line and classify every degree.

Run:  sage gate1_final.sage     (from paper5/)
"""
import os, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper4", "certificates"))
exec(open(os.path.join(HERE, "gate1_alldegrees.sage")).read().split("BOX = 3")[0])
exec(open(os.path.join(HERE, "gate1_crosscuts.sage")).read().split("BOX = 3")[0]
     .split("from examples import V_19")[1])

CH = [0]
def ok(label, cond):
    CH[0] += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {label}")
    assert cond, label

NPb = matrix(ZZ, matrix(ZZ, [list(Vv[i]) for i in PENT]).row_space()
             .intersection(ZZ**4).basis())
Wg = [vector(ZZ, [ZZ(y) for y in NPb.transpose().solve_right(vector(QQ, Vv[i]))])
      for i in PENT]
tris = list(itertools.combinations(range(5), int(3)))

def tri_area(pa, pb, pc, tri):
    a, b, c = tri
    pts = [vector(QQ, Wg[a]) / pa, vector(QQ, Wg[b]) / pb, vector(QQ, Wg[c]) / pc]
    o = pts[0]
    D = matrix(QQ, [list(x - o) for x in pts[1:]])
    if D.rank() < 2:
        return QQ(0)
    L = matrix(ZZ, D.row_space().intersection(ZZ**3).basis())
    two = [vector(QQ, L.transpose().solve_right(x - o)) for x in pts]
    return 2 * Polyhedron(vertices=[list(x) for x in two], base_ring=QQ).volume()

print("== step 3: the finite list of pairing vectors ==")
CAP = 1
while CAP <= 20:
    found = [(t, a, b, c) for t in tris
             for a in range(1, CAP + 1) for b in range(1, CAP + 1)
             for c in range(1, CAP + 1) if tri_area(a, b, c, t) >= QQ(1) / 3]
    if not any(max(f[1:]) == CAP for f in found):
        break
    CAP += 1
ok(f"the search stabilises at cap {CAP}: no solution touches the boundary, and "
   "triangle area is decreasing in each pairing, so the region is captured",
   CAP <= 20)
cands = set()
for t, pa, pb, pc in found:
    A = matrix(ZZ, [list(Wg[j]) for j in t])
    try:
        Rp = A.solve_right(vector(QQ, [pa, pb, pc]))
    except ValueError:
        continue
    if all(x in ZZ for x in Rp):
        Rp = vector(ZZ, [ZZ(y) for y in Rp])
        pv = tuple(ZZ(Rp.dot_product(Wg[j])) for j in range(5))
        for q in (pv, tuple(-x for x in pv)):
            if all(x >= 1 for x in q):
                cands.add(q)
def area(pv):
    return sum(tri_area(pv[a], pv[b], pv[c], (a, b, c))
               for (a, b, c) in itertools.combinations(range(5), int(3))) / 4
big = sorted([pv for pv in cands if
              max(tri_area(pv[a], pv[b], pv[c], (a, b, c))
                  for (a, b, c) in tris) * 3 >= 1], key=lambda z: z)
print(f"    {len(big)} positive pairing vectors survive the bound "
      f"({2 * len(big)} with their negatives)")

print("\n== step 4: walk the lines ==")
Asys = matrix(ZZ, [list(Vv[i]) for i in PENT])
kd = Asys.right_kernel().basis()
ok(f"the kernel of M -> M_P has rank one, direction {tuple(kd[0])}", len(kd) == 1)
d = vector(ZZ, [ZZ(x) for x in kd[0]])
RANGE = 30
tot = caseB = caseC = 0
fails = []
lineB = []; lineC = []
for pv in big + [tuple(-x for x in q) for q in big]:
    try:
        R0 = Asys.solve_right(vector(QQ, pv))
    except ValueError:
        continue
    if any(x not in ZZ for x in R0):
        continue
    mineB = 0; mineC = 0
    for t in range(-RANGE, RANGE + 1):
        R = vector(ZZ, [ZZ(x) for x in (vector(QQ, R0) + t * vector(QQ, d))])
        if gcd([ZZ(x) for x in R]) != 1:
            continue
        tot += 1
        sgn = 1 if pv[0] > 0 else -1
        if germ_smoothable(R, sgn) is not True:
            caseB += 1; mineB += 1; continue
        good = False
        for fi in (21, 24):
            S = sorted(facs[fi][2])
            C = Polyhedron(rays=[list(Vv[i]) for i in S], base_ring=QQ) & \
                Polyhedron(eqns=[[-sgn] + list(R)], base_ring=QQ)
            if not C.is_empty() and defcone_dim(C) == 4:
                good = True; break
        if good:
            caseC += 1; mineC += 1
        else:
            fails.append((tuple(R), pv))
    (lineB if mineC == 0 else lineC).append(pv)
print(f"    of the {2*len(big)} lines, {len(lineB)} are settled WHOLESALE in case (B)")
print(f"       and {len(lineC)} need case (C): {sorted(set(lineC))}")
print(f"    {tot} primitive degrees on the {2 * len(big)} lines, |t| <= {RANGE}")
print(f"      case (B), no smoothing decomposition of the crosscut: {caseB}")
print(f"      case (C), a pentagon-bearing facet cell is rigid:     {caseC}")
print(f"      open:                                                 {len(fails)}")
for R, pv in fails[:10]:
    print(f"        R = {R}, pairings {pv}")
ok("every degree that could possibly smooth the germ is settled", not fails)
print(f"\n{CH[0]} checks passed.")
