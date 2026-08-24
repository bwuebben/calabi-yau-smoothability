"""Ilten-Vollmert Definition 4.1 over a whole two-sided slice complex, at a
vertex degree, for a given polytope / facet / lattice dilation.

Generalises admissibility.sage (which did the prism example) so the same check
runs on the 9-vertex census candidate.

At R = -u_F one has <R,v> >= -1 with equality exactly on F, so the vertices of
negative pairing are exactly F's, every cell of the level-(-1) slice is
conv(G cap F) + tail, and the decomposition extends cellwise by
cell(G)^i = conv(phi_i(G cap F)) + tail.  Definition 4.1(ii) quantifies over
arbitrary collections of cells, not only pairs.  The exact finite check below
builds the intersection semilattice of the two summand subdivisions and tests
every one-cell extension; transitivity of the face relation then covers every
inclusion of collections.  The pairwise counts are retained as diagnostics.
"""
import os, sys, itertools, json
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper4", "certificates"))
from batyrev_global import facets, two_faces, face_lattice_polygon, classify_polygon

CH = [0]
def ok(label, cond):
    CH[0] += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {label}")
    assert cond, label

def run(name, Vt, FI, T):
    print(f"\n########## {name} ##########")
    V = [vector(ZZ, v) for v in Vt]
    facs = facets(Vt)
    FSET = frozenset(facs[FI][2])
    R = vector(ZZ, [-x for x in facs[FI][0]])
    KB = matrix(ZZ, R.row().right_kernel().basis_matrix())
    ok(f"R = {tuple(R)}: <R,v> = -1 exactly on facet {FI} = {sorted(FSET)}, "
       "and >= 0 elsewhere",
       all(R.dot_product(V[i]) == -1 for i in FSET) and
       all(R.dot_product(V[i]) >= 0 for i in range(len(V)) if i not in FSET))

    pf = [sorted(I) for I, fp in two_faces(Vt, facs) if FI in fp]
    edges = sorted({tuple(sorted(set(A) & set(B)))
                    for A, B in itertools.combinations(pf, 2) if len(set(A) & set(B)) == 2})
    eidx = {e: j for j, e in enumerate(edges)}
    ok(f"the facet has {len(edges)} edges and the given dilation has that length",
       len(T) == len(edges))

    adj = {}
    for (a, b) in edges:
        adj.setdefault(a, []).append(b); adj.setdefault(b, []).append(a)
    root = sorted(FSET)[0]
    phi0 = {root: vector(QQ, [0, 0, 0, 0])}
    st = [root]
    while st:
        a = st.pop()
        for b in adj[a]:
            w = phi0[a] + T[eidx[tuple(sorted((a, b)))]] * (V[b] - V[a])
            if b in phi0:
                assert phi0[b] == w, "not path independent"
            else:
                phi0[b] = w; st.append(b)
    phi1 = {a: (V[a] - V[root]) - phi0[a] for a in phi0}
    A = Polyhedron(vertices=[list(phi0[a]) for a in phi0], base_ring=QQ)
    B = Polyhedron(vertices=[list(phi1[a]) for a in phi1], base_ring=QQ)
    QF = Polyhedron(vertices=[list(V[a] - V[root]) for a in sorted(FSET)], base_ring=QQ)
    ok(f"A + B = F   (A: {A.n_vertices()} vertices, B: {B.n_vertices()})", A + B == QF)

    H = Polyhedron(eqns=[[1] + list(R)], base_ring=QQ)
    n = len(Vt)
    seen, frontier = {frozenset(range(n))}, [frozenset(range(n))]
    while frontier:
        nxt = []
        for S in frontier:
            for _, _, I in facs:
                Tt = S & I
                if Tt and Tt not in seen:
                    seen.add(Tt); nxt.append(Tt)
        frontier = nxt
    cells = {}
    for S in seen:
        if S == frozenset(range(n)):
            continue
        C = Polyhedron(rays=[list(V[i]) for i in sorted(S)], base_ring=QQ) & H
        if C.is_empty():
            continue
        onF = sorted(S & FSET)
        tail = [list(r) for r in C.rays()]
        cells[S] = (C,
                    Polyhedron(vertices=[list(phi0[a]) for a in onF], rays=tail, base_ring=QQ),
                    Polyhedron(vertices=[list(phi1[a]) for a in onF], rays=tail, base_ring=QQ))
    print(f"    {len(cells)} nonempty cells in the level-(-1) slice")

    shift = Polyhedron(vertices=[list(-V[root])], base_ring=QQ)
    rc = lambda P: Polyhedron(rays=[list(r) for r in P.rays()], base_ring=QQ)
    # NB the tail cones agree by CONSTRUCTION (both summands are built with
    # C.rays()), so testing them proves nothing -- defect D7.  What has content
    # is that the summands add up.
    bad = [S for S, (C, S0, S1) in cells.items() if (S0 + S1) != (C + shift)]
    ok(f"Definition 2.1: cell^0 + cell^1 = cell on all {len(cells)} cells "
       f"({len(bad)} failures); the tail cones agree by construction",
       not bad)
    # Definition 2.2 on EVERY cell, not just on the facet (defect D8)
    nonlat = []
    for S, (C, S0, S1) in cells.items():
        for P in (S0, S1):
            if any(x not in ZZ for v in P.vertices_list() for x in v):
                nonlat.append(sorted(S))
    ok(f"Definition 2.2 on all {len(cells)} cells: every summand has lattice "
       f"vertices, so no vertex is non-lattice anywhere ({len(nonlat)} "
       "failures)", not nonlat)

    bad_i = []
    for S, T2 in itertools.combinations(list(cells), int(2)):
        CS, S0, S1 = cells[S]; CT, T0, T1 = cells[T2]
        inter = CS & CT
        if inter.is_empty():
            continue
        key = S & T2
        if key not in cells or cells[key][0] != inter or \
           (S0 & T0) != cells[key][1] or (S1 & T1) != cells[key][2]:
            bad_i.append((sorted(S), sorted(T2)))
    ok(f"Definition 4.1(i): (D cap E)^i = D^i cap E^i ({len(bad_i)} failures)", not bad_i)

    # --- D1: Definition 4.1(i) is conditional on a nonempty intersection, so
    # the reduction of (ii) does not cover families whose cells are DISJOINT
    # -- and their summands can still meet.  Those cases are tested directly.
    def is_face(P, Q):
        # The empty polyhedron is allowed as a face in Definition 4.1(ii).
        if P.is_empty():
            return True
        if Q.is_empty():
            return False
        if P == Q:
            return True
        for d in range(Q.dim() + 1):
            for f in Q.faces(d):
                if f.as_polyhedron() == P:
                    return True
        return False

    def poly_key(P):
        if P.is_empty():
            return ("empty",)
        vec = lambda v: tuple(QQ(x) for x in v)
        return (tuple(sorted(vec(v) for v in P.vertices_list())),
                tuple(sorted(vec(v) for v in P.rays_list())),
                tuple(sorted(vec(v) for v in P.lines_list())))

    def exact_collection_check(label, data):
        """Check IV Definition 4.1(ii) for every collection of cells.

        A state is the pair of intersections of the 0- and 1-summands over a
        nonempty collection B of cells.  Adding one cell produces the state for
        A = B union {cell}.  For each of the three nonempty subsets of summand
        indices, Definition 4.1(ii) says that the corresponding sum for A is a
        face of the sum for B.  Every inclusion B subset A factors into such
        additions, and face relations are transitive.
        """
        summands = [(x[1], x[2]) for x in data.values()]
        states, queue = {}, []

        def add_state(pair):
            key = (poly_key(pair[0]), poly_key(pair[1]))
            if key not in states:
                states[key] = pair
                queue.append(key)

        def selected_sum(pair, mask):
            if mask == 1:
                return pair[0]
            if mask == 2:
                return pair[1]
            return pair[0] + pair[1]

        for pair in summands:
            add_state(pair)

        bad, transitions, head = [], 0, 0
        while head < len(queue):
            old_key = queue[head]
            head += 1
            old = states[old_key]
            for cell_no, cell in enumerate(summands):
                new = (old[0] & cell[0], old[1] & cell[1])
                transitions += 1
                for mask in (1, 2, 3):
                    if not is_face(selected_sum(new, mask),
                                   selected_sum(old, mask)):
                        bad.append((old_key, cell_no, mask,
                                    (poly_key(new[0]), poly_key(new[1]))))
                add_state(new)

        ok(f"{label} Definition 4.1(ii), arbitrary cell collections: "
           f"{len(states)} intersection states, {transitions} one-cell "
           f"extensions, {len(bad)} violations", not bad)
        return len(states), transitions
    bad_ii, pairs, notface = [], 0, []
    for S, T2 in itertools.permutations(list(cells), int(2)):
        if not (S < T2):
            continue
        CS, S0, S1 = cells[S]; CT, T0, T1 = cells[T2]
        if not is_face(CS, CT):
            # a contained-but-not-face pair would itself violate 4.1(ii) at
            # I = {0,1}; report it rather than skip it (defect D6)
            notface.append((sorted(S), sorted(T2)))
            continue
        pairs += 1
        if not is_face(S0, T0):
            bad_ii.append((0, sorted(S), sorted(T2)))
        if not is_face(S1, T1):
            bad_ii.append((1, sorted(S), sorted(T2)))
    ok(f"no cell is contained in another without being a face of it "
       f"({len(notface)} such pairs) -- such a pair would itself break "
       "4.1(ii) at I = {0,1}", not notface)
    ok(f"Definition 4.1(ii), I = {{0}} and I = {{1}}: kappa^i is a face of "
       f"lambda^i on all {pairs} face pairs ({len(bad_ii)} failures)",
       not bad_ii)
    keys = list(cells)
    dis_pairs = dis_hit = dis_bad = 0
    for a, b in itertools.combinations(keys, int(2)):
        Ca, A0, A1 = cells[a]; Cb, B0, B1 = cells[b]
        if not (Ca & Cb).is_empty():
            continue
        dis_pairs += 1
        I0, I1 = A0 & B0, A1 & B1
        if I0.is_empty() and I1.is_empty():
            continue
        dis_hit += 1
        for (P0, P1, Q0, Q1) in ((I0, I1, A0, A1), (I0, I1, B0, B1)):
            for pair in ((P0, Q0), (P1, Q1)):
                if not pair[0].is_empty() and not is_face(pair[0], pair[1]):
                    dis_bad += 1
            if not P0.is_empty() and not P1.is_empty():
                if not is_face(P0 + P1, Q0 + Q1):
                    dis_bad += 1
    ok(f"Definition 4.1(ii) on DISJOINT pairs, which (i) cannot reach: "
       f"{dis_pairs} disjoint pairs, {dis_hit} with a nonempty summand "
       f"intersection, {dis_bad} violations", dis_bad == 0)
    exact_collection_check("level-(-1)", cells)

    # The level-(+1) slice is deliberately left undecomposed.  Earlier
    # versions merely said this and never constructed that half of the slice
    # complex (defect D-G).  Build it and check the trivial decomposition
    # C = C + tail(C) against the same compatibility conditions.
    Hp = Polyhedron(eqns=[[-1] + list(R)], base_ring=QQ)
    pcells = {}
    for S in seen:
        if S == frozenset(range(n)):
            continue
        C = Polyhedron(rays=[list(V[i]) for i in sorted(S)], base_ring=QQ) & Hp
        if C.is_empty():
            continue
        Tcone = Polyhedron(vertices=[[0, 0, 0, 0]],
                           rays=[list(r) for r in C.rays()], base_ring=QQ)
        pcells[S] = (C, C, Tcone)
    print(f"    {len(pcells)} nonempty cells in the level-(+1) slice")

    pbad_sum = [S for S, (C, S0, S1) in pcells.items() if S0 + S1 != C]
    ok(f"level-(+1) Definition 2.1: the trivial decomposition adds up on "
       f"all {len(pcells)} cells ({len(pbad_sum)} failures)", not pbad_sum)

    # Definition 2.2 allows at most one non-lattice face in each evaluated
    # pair.  In C + tail(C), the second summand always contains the lattice
    # point 0 in every character-minimising face, so the condition is
    # automatic even when C has rational vertices.
    pbad_adm = [S for S, (_, _, S1) in pcells.items()
                if vector(QQ, [0, 0, 0, 0]) not in S1]
    ok(f"level-(+1) Definition 2.2: the tail summand contains the lattice "
       f"point 0 on all {len(pcells)} cells ({len(pbad_adm)} failures)",
       not pbad_adm)

    pbad_i = []
    for S, T2 in itertools.combinations(list(pcells), int(2)):
        CS, S0, S1 = pcells[S]; CT, T0, T1 = pcells[T2]
        inter = CS & CT
        if inter.is_empty():
            continue
        key = S & T2
        if key not in pcells or pcells[key][0] != inter or \
           (S0 & T0) != pcells[key][1] or (S1 & T1) != pcells[key][2]:
            pbad_i.append((sorted(S), sorted(T2)))
    ok(f"level-(+1) Definition 4.1(i) on every nonempty cell "
       f"intersection ({len(pbad_i)} failures)", not pbad_i)

    pbad_ii, ppairs, pnotface = [], 0, []
    for S, T2 in itertools.permutations(list(pcells), int(2)):
        if not (S < T2):
            continue
        CS, S0, S1 = pcells[S]; CT, T0, T1 = pcells[T2]
        if not is_face(CS, CT):
            pnotface.append((sorted(S), sorted(T2)))
            continue
        ppairs += 1
        if not is_face(S0, T0):
            pbad_ii.append((0, sorted(S), sorted(T2)))
        if not is_face(S1, T1):
            pbad_ii.append((1, sorted(S), sorted(T2)))
    ok(f"level-(+1): no contained cell fails to be a face "
       f"({len(pnotface)} failures)", not pnotface)
    ok(f"level-(+1) Definition 4.1(ii) on all {ppairs} face pairs "
       f"({len(pbad_ii)} failures)", not pbad_ii)

    pdis_pairs = pdis_hit = pdis_bad = 0
    for a, b in itertools.combinations(list(pcells), int(2)):
        Ca, A0, A1 = pcells[a]; Cb, B0, B1 = pcells[b]
        if not (Ca & Cb).is_empty():
            continue
        pdis_pairs += 1
        I0, I1 = A0 & B0, A1 & B1
        if I0.is_empty() and I1.is_empty():
            continue
        pdis_hit += 1
        for (P0, P1, Q0, Q1) in ((I0, I1, A0, A1), (I0, I1, B0, B1)):
            for pair in ((P0, Q0), (P1, Q1)):
                if not pair[0].is_empty() and not is_face(pair[0], pair[1]):
                    pdis_bad += 1
            if not P0.is_empty() and not P1.is_empty() and \
               not is_face(P0 + P1, Q0 + Q1):
                pdis_bad += 1
    ok(f"level-(+1) Definition 4.1(ii) on {pdis_pairs} disjoint pairs, "
       f"{pdis_hit} with a nonempty summand intersection "
       f"({pdis_bad} violations)", pdis_bad == 0)
    exact_collection_check("level-(+1)", pcells)

    print("\n  == Corollary 2.12 on the singular cells ==")
    def cone_smooth(P, KB):
        """Is TV(Cone(Q x {1})) smooth?  The cell lives in N' = N cap R^perp;
        express it in the BASIS KB of that lattice rather than by deleting a
        coordinate -- deleting one is a lattice isomorphism only when the
        degree has a unit coordinate, which is false for three of Delta_9's
        twelve vertex degrees (defect D5)."""
        def c(x):
            return list(KB.transpose().solve_right(vector(QQ, x)))
        gens = [vector(QQ, c(p) + [1]) for p in P.vertices_list()]
        gens += [vector(QQ, c(r) + [0]) for r in P.rays()]
        gens = [g * lcm([QQ(x).denominator() for x in g]) for g in gens]
        gens = [vector(ZZ, [ZZ(x) / gcd([ZZ(y) for y in g]) for x in g]) for g in gens]
        M = matrix(ZZ, gens); r = M.rank()
        if len(gens) != r:
            return False
        return gcd(M.minors(r)) == 1
    smoothed, left = [], []
    for I, fp in two_faces(Vt, facs):
        _, evs, lens = face_lattice_polygon(Vt, I, facs[fp[0]][0], facs[fp[1]][0])
        cl = classify_polygon(evs, lens)
        if cl["status"] == "smooth":
            continue
        S = frozenset(I)
        kind = {4: "node", 5: "dP7 ", 6: "dP6 "}.get(cl["k"], "?")
        if S not in cells:
            print(f"    {kind} {sorted(I)}: no cell at this level (unreachable)")
            left.append(sorted(I)); continue
        C, S0, S1 = cells[S]
        if not C.is_compact():
            print(f"    {kind} {sorted(I)}: unbounded cell")
        sm = cone_smooth(S0, KB) and cone_smooth(S1, KB)
        print(f"    {kind} {sorted(I)}: summands {S0.n_vertices()}+{S1.n_vertices()} "
              f"vertices -> {'SMOOTH' if sm else 'still singular'}")
        (smoothed if sm else left).append(sorted(I))
    ok(f"{len(smoothed)} of {len(smoothed) + len(left)} singular 2-faces are "
       f"smoothed by this one family", True)
    return len(left) == 0

# ---------------------------------------------------------------- the runs
P5 = [(1, 0), (0, 1), (-1, -1), (-1, 0), (0, -1)]
V_P = [tuple(int(x) for x in v) for v in
       [(p[0], p[1], z, 1) for p in P5 for z in (0, 1)] + [(0, 0, -1, -1), (0, 0, 0, -1)]]
run("Delta_P (prism example), facet 0", V_P, 0,
    (0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0))

W = json.load(open(os.path.join(HERE, "v09_candidate.json")))
V_C = [tuple(int(x) for x in v) for v in W["V"]]
full = run("the 9-vertex census candidate, facet 11", V_C, 11,
           (0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0))
print(f"\n{CH[0]} checks passed.")
if full:
    print("""
*** EVERY singular 2-face of the 9-vertex candidate is smoothed by ONE
single-degree Ilten-Vollmert family.  The ambient toric fourfold's whole
singular locus along the hypersurface's germs is removed in the general
fibre. ***""")
