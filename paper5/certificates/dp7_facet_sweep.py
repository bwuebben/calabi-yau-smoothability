#!/usr/bin/env python3
"""Is the gate-1 rigidity special to Delta_19, or does it always happen?

For a reflexive 4-polytope Delta with a pentagon 2-face P (five unit
edges, one interior point -> a dP7 cone point on the generic anticanonical
hypersurface), gate1.md shows for Delta_19 that the bounded cell of the
slice complex carrying P is Minkowski-rigid at EVERY degree, because a
chain of 2-faces forces all edge dilations of the pentagon-bearing facet
equal.  The chain is purely combinatorial, so the same test runs over the
Kreuzer-Skarke database.

Test per polytope: for each pentagon 2-face and each of the two facets
containing it, propagate

  * a triangular 2-face forces its three edge dilations equal;
  * a 2-face with exactly one unforced edge forces it;
  * a 2-face with exactly two unforced edges, ADJACENT to each other,
    forces both (adjacent edges of a polygon are never parallel).

If the chain reaches all edges of the facet, that facet -- and hence the
slice cell over it, which has the facet's face lattice whenever all its
rays pair negatively -- is indecomposable for every degree on the
pentagon's dual-edge line.  If the chain stalls, the polytope is a
CANDIDATE: the cell may decompose, and a single-degree Ilten-Vollmert
family may be able to reach its dP7 point.  Candidates are printed with
their vertices for follow-up; a stalled chain is not by itself a proof of
decomposability.

Run:  ./venv/bin/python dp7_facet_sweep.py ../data/ks/polytopes-4d-0{5,6,7}-vertices.parquet
      ./venv/bin/python dp7_facet_sweep.py ../data/ks/*.parquet --procs 8
"""
import argparse
import itertools
import json
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper4", "certificates"))

from ks_sweep import (facets_reflexive, two_faces_of, polygon_of_face,   # noqa
                      iter_vert_chunks, file_meta)
from batyrev_global import vgcd, vsub                                    # noqa
from examples import rank                                                # noqa


def pick_interior(evs, lens):
    """(2*area, boundary points, interior points) of the lattice polygon."""
    A2 = x = y = 0
    for e, l in zip(evs, lens):
        A2 += x * (l * e[1]) - y * (l * e[0])
        x += l * e[0]; y += l * e[1]
    A2 = abs(A2)
    b = sum(lens)
    return A2, b, (A2 - b + 2) // 2


def facet_complex(twofaces_in_facet):
    """Edges of a 3-polytope from its list of 2-faces (as vertex sets), plus
    each 2-face in cyclic order.  Returns None if the incidence data is not
    that of a 3-polytope (guards against degenerate input)."""
    edges = set()
    for A, B in itertools.combinations(twofaces_in_facet, 2):
        common = tuple(sorted(set(A) & set(B)))
        if len(common) == 2:
            edges.add(common)
    cycles = []
    for A in twofaces_in_facet:
        adj = {}
        for e in edges:
            if set(e) <= set(A):
                adj.setdefault(e[0], []).append(e[1])
                adj.setdefault(e[1], []).append(e[0])
        if sorted(adj) != sorted(A) or any(len(adj[v]) != 2 for v in A):
            return None, None
        cyc = [A[0]]; prev, cur = None, A[0]
        while len(cyc) < len(A):
            nxt = [x for x in adj[cur] if x != prev][0]
            cyc.append(nxt); prev, cur = cur, nxt
        if len(set(cyc)) != len(A):
            return None, None
        cycles.append(cyc)
    return sorted(edges), cycles


def chain_closes(twofaces_in_facet):
    """True iff the 2-face chain forces every edge dilation equal.

    UNION-FIND, so that the constant each edge is forced to is tracked.  An
    earlier version of this routine kept a flat set of forced edges and could
    not tell whether the already-forced edges of a 2-face carried the SAME
    constant; without that check both the single-edge and the adjacent-pair
    rule are unsound, and an adversarial round produced 3-polytopes where the
    flat version reported rigidity for a decomposable polytope.  The census
    numbers in the first version of this sweep were computed with the flat
    version and are superseded by the ones this produces.
    """
    edges, cycles = facet_complex(twofaces_in_facet)
    if edges is None:
        return None
    par = {e: e for e in edges}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: par[ra] = rb

    seed = None
    for A, cyc in zip(twofaces_in_facet, cycles):
        if len(A) == 3:
            es = [tuple(sorted((cyc[j], cyc[(j + 1) % 3]))) for j in range(3)]
            union(es[0], es[1]); union(es[1], es[2])
            if seed is None: seed = es[0]
    if seed is None:
        return False
    progress = True
    while progress:
        progress = False
        big = find(seed)
        for cyc in cycles:
            fe = [tuple(sorted((cyc[j], cyc[(j + 1) % len(cyc)])))
                  for j in range(len(cyc))]
            rest = [e for e in fe if find(e) != big]
            if not rest:
                continue
            if len(rest) == 1:
                union(rest[0], big); big = find(big); progress = True
            elif len(rest) == 2 and len(set(rest[0]) & set(rest[1])) == 1:
                union(rest[0], big); union(rest[1], big)
                big = find(big); progress = True
    return len({find(e) for e in edges}) == 1


def summand_dim(V, twofaces_in_facet):
    """dim of the space of edge dilations closing on every 2-face of the
    facet, computed on the lattice facet itself.  1 = indecomposable."""
    from fractions import Fraction as Fr
    edges, cycles = facet_complex(twofaces_in_facet)
    if edges is None:
        return None
    eidx = {e: k for k, e in enumerate(edges)}
    rows = []
    for cyc in cycles:
        for c in range(4):
            row = [Fr(0)] * len(edges)
            for j in range(len(cyc)):
                a, b = cyc[j], cyc[(j + 1) % len(cyc)]
                row[eidx[tuple(sorted((a, b)))]] += V[b][c] - V[a][c]
            rows.append(row)
    return len(edges) - rank(rows)


def analyse(V):
    """Per polytope: one record per pentagon 2-face, or [] if none."""
    V = [tuple(int(x) for x in v) for v in V]
    facs = facets_reflexive(V)
    tf = two_faces_of(V, facs)
    pentagons = []
    in_framework = True          # every 2-face inside the Paper 4 inventory?
    for I, (a, b) in tf.items():
        evs, lens = polygon_of_face(V, sorted(I), facs[a][0], facs[b][0])
        A2, _, i = pick_interior(evs, lens)
        k = len(evs)
        dl = vgcd(vsub(facs[a][0], facs[b][0]))
        # a 2-face is harmless only if it is a UNIMODULAR triangle; otherwise
        # it must be a unit-edge square / pentagon / hexagon carrying one point
        if not (k == 3 and A2 == 1):
            if any(l != 1 for l in lens) or dl != 1 or \
               not ((k == 4 and i == 0) or (k in (5, 6) and i == 1)):
                in_framework = False
        if k == 5 and i == 1 and all(l == 1 for l in lens):
            pentagons.append((sorted(I), a, b))
    if not pentagons:
        return []
    by_facet = {}
    for I, (a, b) in tf.items():
        by_facet.setdefault(a, []).append(sorted(I))
        by_facet.setdefault(b, []).append(sorted(I))
    out = []
    for I, a, b in pentagons:
        dl = vgcd(vsub(facs[a][0], facs[b][0]))
        verdicts, dims = [], []
        for f in (a, b):
            c = chain_closes(by_facet[f])
            verdicts.append(c)
            # the chain is only a sufficient test; decide the stalled ones
            dims.append(None if c else summand_dim(V, by_facet[f]))
        out.append(dict(pent=I, dual_len=int(dl), verdicts=verdicts, dims=dims,
                        framework=in_framework, _facets=[int(a), int(b)],
                        nfaces=[len(by_facet[a]), len(by_facet[b])]))
    return out


def _work(rows):
    agg = Counter()
    cands = []
    fws = []
    for V in rows:
        for r in analyse(V):
            agg["pentagon_faces"] += 1
            if r["framework"]:
                agg["IN_PAPER4_FRAMEWORK"] += 1
                rec = {k: v for k, v in r.items() if k != "_facets"}
                fws.append(dict(V=[[int(x) for x in y] for y in V],
                                facets=[int(x) for x in r["_facets"]], **rec))
            agg[f"dual_len_{min(r['dual_len'], 3)}"] += 1
            v = r["verdicts"]
            if None in v:
                agg["degenerate_incidence"] += 1
                if r["framework"]:
                    agg["DEGENERATE_AND_IN_FRAMEWORK"] += 1
            elif all(v):
                agg["rigid_both_facets"] += 1
                if r["framework"]:
                    agg["LOCKED_AND_IN_FRAMEWORK"] += 1
            elif all(d == 1 for d in r["dims"] if d is not None):
                agg["chain_stalled_but_rigid"] += 1
                if r["framework"]:
                    agg["STALLED_BUT_RIGID_AND_IN_FRAMEWORK"] += 1
            else:
                agg["DECOMPOSABLE_CELL"] += 1
                if r["framework"]:
                    agg["DECOMPOSABLE_AND_IN_FRAMEWORK"] += 1
                if len(cands) < 400:
                    facet = [f for f, d in zip(r["_facets"], r["dims"])
                             if d is not None and d >= 2][0]
                    rec = {k: v for k, v in r.items() if k != "_facets"}
                    cands.append(dict(V=[[int(x) for x in y] for y in V],
                                      facet=int(facet), **rec))
    if cands:
        agg["_cands"] = cands
    if fws:
        agg["_fw"] = fws
    return agg


def selftest():
    from examples import V_19
    r = analyse([tuple(v) for v in V_19])
    assert len(r) == 1, r
    assert r[0]["pent"] == [14, 15, 16, 17, 18], r
    assert r[0]["dual_len"] == 1, r
    assert r[0]["verdicts"] == [True, True], r
    assert r[0]["nfaces"] == [9, 9], r
    print("  selftest: Delta_19 reproduces gate1.md "
          "(one pentagon, dual length 1, chain closes on both facets)  [ok]")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--procs", type=int, default=1)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out", default="dp7_facet_sweep.json")
    args = ap.parse_args()

    print("== self-test ==")
    selftest()

    total = Counter()
    cands = []
    fwall = []
    for path in args.files:
        n_rows, n_verts = file_meta(path)
        agg = Counter()
        got = 0
        chunks = iter_vert_chunks(path, 2000, args.limit)
        if args.procs > 1:
            from multiprocessing import Pool
            with Pool(args.procs) as pool:
                for a in pool.imap_unordered(_work, chunks):
                    cands.extend(a.pop("_cands", []))
                    fwall.extend(a.pop("_fw", []))
                    agg.update(a); got += 1
        else:
            for rows in chunks:
                a = _work(rows)
                cands.extend(a.pop("_cands", []))
                fwall.extend(a.pop("_fw", []))
                agg.update(a); got += 1
        print(f"\n{os.path.basename(path)}: {n_rows} polytopes, {n_verts} vertices")
        for k in sorted(agg):
            print(f"    {k:<24} {agg[k]}")
        total.update(agg)

    print("\n== totals ==")
    for k in sorted(total):
        print(f"    {k:<24} {total[k]}")
    if cands:
        tmp = args.out + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(cands, fh, indent=1)
        os.replace(tmp, args.out)
        print(f"\n{len(cands)} candidate(s) written to {args.out}")
    if fwall:
        fwpath = (args.out or "fw.json").replace(".json", "_framework.json")
        with open(fwpath, "w") as fh:
            json.dump(fwall, fh, indent=1)
        print(f"{len(fwall)} framework pentagon(s) written to {fwpath}")
    elif not total["pentagon_faces"]:
        print("\nNo pentagon 2-face occurs in this range at all -- the sweep "
              "is vacuous here.")
    else:
        print(f"\nAll {total['pentagon_faces']} pentagon 2-faces in this range "
              "have BOTH their facets\nforced rigid by the chain: the gate-1 "
              "rigidity is not special to Delta_19.")


if __name__ == "__main__":
    main()
