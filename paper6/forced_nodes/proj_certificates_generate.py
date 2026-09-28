#!/usr/bin/env python3
"""Generate the data file list_polytopes.json.gz read by forced_nodes_check.py: the admissible polytopes of the list,
their two-faces, and (Proj) witnesses and Farkas certificates for every statement counted in Section 11.4.

    sage -python paper6/forced_nodes/proj_certificates_generate.py [--cap SEC] [--only ID ...]

Run from the repository root.  SageMath is used only for exact rational linear programming (the PPL backend of
MixedIntegerLinearProgram); faces, nodes, circuits and K are computed by the functions of forced_nodes_check.py.  The
checker does not need this program: it verifies everything this program writes.  How the vectors were found plays no
role in the check.

The linear system of a (partial) profile (Lemma "(Proj) as a linear program" of the paper): unknowns psi in Q^R,
equations <circ_q, psi> = 0 for the node parallelograms of the squares, the pentagons and the hexagons with a
prescribed subdivision, inequalities <I_eps, psi> >= 1 for the interior edges of the subdivisions of the pentagons
and of these hexagons.  It is infeasible if and only if some y >= 0 with sum(y) = 1 and some nu satisfy
sum y_eps I_eps = sum nu_q circ_q (Farkas), and then every profile extending the partial one is infeasible.

For a statement "(Proj) holds for some profile in a set S" (S = the profiles with every hexagon on A1, on A2, with
both branches, or with a branch pattern without forced node), a depth-first search over the subdivisions of the
hexagons either finds a profile in S with an integral witness psi, or exhausts S; each subtree is cut off by a
Farkas certificate of its partial profile, reduced greedily to as few hexagons as possible, and recorded as a nogood
(partial profile, integral y, integral nu).  The nogoods then cover S.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
from fractions import Fraction
from itertools import permutations, product
from math import lcm
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import forced_nodes_check as F  # noqa: E402
from sage.all import MixedIntegerLinearProgram, QQ, matrix  # noqa: E402
from sage.numerical.mip import MIPSolverException  # noqa: E402

T0 = time.time()
PAPER_NAMES = {"paper4/certificates/framework_77.json[50]": "Delta_9",
               "paper4/certificates/examples.py:V_19": "Delta_19",
               "paper4/certificates/examples.py:V_20": "Delta_20"}


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


# ------------------------------------------------------------------------------------------------------------------
# exact linear programming

def farkas(ineq, eq):
    """ineq, eq: lists of (id, row).  Returns (y, nu) as dicts of integers with sum y I = sum nu circ, y >= 0 not
    zero, or None if the system <I, psi> >= 1, <circ, psi> = 0 is feasible."""
    if not ineq:
        return None
    p = MixedIntegerLinearProgram(solver="PPL", maximization=False)
    y = p.new_variable(nonnegative=True)
    nu = p.new_variable(nonnegative=False)
    # Materialize every column before solve: indexing afterwards changes the LP.
    for j in range(len(eq)):
        _ = nu[j]
    n = len(ineq[0][1])
    p.add_constraint(sum(y[i] for i in range(len(ineq))) == 1)
    for c in range(n):
        ti = [i for i in range(len(ineq)) if ineq[i][1][c]]
        tj = [j for j in range(len(eq)) if eq[j][1][c]]
        if not ti and not tj:
            continue                                   # coordinate c appears in no row
        p.add_constraint(sum(ineq[i][1][c] * y[i] for i in ti) - sum(eq[j][1][c] * nu[j] for j in tj) == 0)
    p.set_objective(sum(y[i] for i in range(len(ineq))))
    try:
        p.solve()
    except MIPSolverException:
        return None
    yv = [Fraction(int(QQ(p.get_values(y[i])).numerator()), int(QQ(p.get_values(y[i])).denominator()))
          for i in range(len(ineq))]
    nv = [Fraction(int(QQ(p.get_values(nu[j])).numerator()), int(QQ(p.get_values(nu[j])).denominator()))
          for j in range(len(eq))]
    L = lcm(*[x.denominator for x in yv + nv]) if yv + nv else 1
    Y = {ineq[i][0]: int(yv[i] * L) for i in range(len(ineq)) if yv[i]}
    N = {eq[j][0]: int(nv[j] * L) for j in range(len(eq)) if nv[j]}
    return Y, N


def witness(ineq, eq, nR):
    p = MixedIntegerLinearProgram(solver="PPL", maximization=False)
    x = p.new_variable(nonnegative=False)
    for c in range(nR):
        _ = x[c]
    if not ineq:
        return [0] * nR
    for _, r in eq:
        if any(r):
            p.add_constraint(sum(r[c] * x[c] for c in range(nR) if r[c]) == 0)
    for _, r in ineq:
        if not any(r):
            return None                                # 0 >= 1 is infeasible
        p.add_constraint(sum(r[c] * x[c] for c in range(nR) if r[c]) >= 1)
    p.set_objective(0)
    try:
        p.solve()
    except MIPSolverException:
        return None
    v = [Fraction(int(QQ(p.get_values(x[c])).numerator()), int(QQ(p.get_values(x[c])).denominator()))
         for c in range(nR)]
    L = lcm(*[t.denominator for t in v]) if v else 1
    return [int(t * L) for t in v]


# ------------------------------------------------------------------------------------------------------------------
# searches

def system(P, prof):
    return P.edges(prof), [(q, P.circuit(abcd)) for q, abcd in P.nodes(prof)]


def reduce_nogood(P, prof, cert):
    """Greedily drop hexagons from a partial profile while the reduced system stays infeasible."""
    keep = dict(prof)
    sup = F.owners(list(cert[0]) + list(cert[1]))
    keep = {h: c for h, c in keep.items() if int(h) in sup}
    best = farkas(*system(P, keep))
    assert best is not None
    for h in sorted(keep, key=int):
        trial = {k: c for k, c in keep.items() if k != h}
        c2 = farkas(*system(P, trial))
        if c2 is not None:
            keep, best = trial, c2
    return keep, best


def order_hexes(P):
    hx = list(range(len(P.hexes)))
    verts = {t: set(P.hexes[t]["cycle"]) for t in hx}
    out = []
    left = list(hx)
    while left:
        nxt = left[0] if not out else max(left, key=lambda t: sum(len(verts[t] & verts[u]) for u in out))
        out.append(nxt)
        left.remove(nxt)
    return out


def search(P, allowed, options, cap):
    """Depth-first search for a profile whose pattern is allowed and which satisfies (Proj)."""
    order = order_hexes(P)
    nogoods = []
    st = dict(lp=0, found=None, capped=False)
    ts = time.time()

    def covered(assign):
        return any(all(assign.get(h) == c for h, c in ng["profile"].items()) for ng in nogoods)

    def rec(depth, assign):
        if st["found"] or st["capped"]:
            return
        if time.time() - ts > cap:
            st["capped"] = True
            return
        if covered(assign):
            return
        pat = tuple(assign[str(h)].split(":")[0] for h in order[:depth])
        full = depth == len(order)
        if not allowed(pat, full):
            return
        st["lp"] += 1
        cert = farkas(*system(P, assign))
        if cert is not None:
            prof, (y, nu) = reduce_nogood(P, assign, cert)
            nogoods.append(dict(profile=prof, y=y, nu=nu))
            return
        if full:
            psi = witness(*system(P, assign), len(P.R))
            assert psi is not None
            F.check_witness(P, assign, psi)
            st["found"] = dict(profile=dict(assign), psi=psi)
            return
        h = order[depth]
        for c in options:
            assign[str(h)] = c
            rec(depth + 1, assign)
            del assign[str(h)]

    rec(0, {})
    assert not st["capped"], "search capped"
    if st["found"]:
        return dict(witness=st["found"]), st["lp"]
    return dict(order=order, nogoods=nogoods), st["lp"]


def nodes_json(nd):
    return [[q, [list(p) for p in abcd]] for q, abcd in nd]


def good_patterns(P, mixed_only=False):
    hx = list(range(len(P.hexes)))
    assert len(hx) <= 12, "too many hexagons to enumerate branch patterns"
    good = []
    for pat in product(("A1", "A2"), repeat=len(hx)):
        if mixed_only and len(set(pat)) < 2:
            continue
        prof = {str(h): ("A1:0" if b == "A1" else "A2:0") for h, b in zip(hx, pat)}
        if not P.lattice(prof)[2]:
            good.append(pat)
    return good


def find_equivalence(A, B):
    """A matrix g in GL_4(Z) with g(A) = B (vertex sets), or None; search over images of a basis of vertices."""
    A = [tuple(v) for v in A]
    Bs = {tuple(v) for v in B}
    if len(A) != len(Bs):
        return None
    base = None
    for S in permutations(range(len(A)), 4):
        M = matrix(QQ, [A[i] for i in S])
        if M.det() != 0:
            base = S
            break
        if S[0] > 0:
            break
    Minv = matrix(QQ, [A[i] for i in base]).inverse()
    Bl = sorted(Bs)
    for T in permutations(range(len(Bl)), 4):
        G = (Minv * matrix(QQ, [Bl[j] for j in T])).transpose()   # g with g a_i = b_i (columns)
        if any(x.denominator() != 1 for x in G.list()) or abs(G.det()) != 1:
            continue
        img = {tuple(int(x) for x in (G * matrix(QQ, 4, 1, list(a))).list()) for a in A}
        if img == Bs:
            return [[int(x) for x in row] for row in G.rows()]
    return None


# ------------------------------------------------------------------------------------------------------------------

def short_id(label):
    if label in PAPER_NAMES:
        return PAPER_NAMES[label]
    if label.startswith("paper4/certificates/framework_77.json["):
        return "framework:" + label.split("[")[1].rstrip("]")
    if "list number" in label:
        n = int(label.split("list number ")[1].rstrip(")"))
        return {88: "Delta_88", 154: "Delta_154"}.get(n, f"list:{n}")
    if label.startswith("polar of"):
        return "X_circ"
    raise ValueError(label)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cap", type=float, default=900)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--prepare-only", action="store_true")
    ap.add_argument("--hints", type=Path)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out", default=str(HERE / "list_polytopes.json.gz"))
    A = ap.parse_args()
    hints = json.loads(A.hints.read_text())["witnesses"] if A.hints else {}
    src = F.list_from_sources()
    merged = {}
    for label, V in src:
        merged.setdefault(frozenset(V), (V, []))[1].append(label)
    log(f"{len(src)} entries, {len(merged)} vertex sets")
    prepared = Path(A.out).with_suffix(".prepared.json")
    checkpoint = Path(A.out).with_suffix(".checkpoint.json")
    completed = json.loads(checkpoint.read_text()) if checkpoint.exists() else {}
    old = {}
    if Path(A.out).exists():
        for r in json.load(gzip.open(A.out, "rt"))["polytopes"]:
            if "stored" in r:
                old[r["id"]] = r["stored"]
    polys = json.loads(prepared.read_text()) if prepared.exists() else []
    for key, (V, labels) in ([] if polys else merged.items()):
        fa = F.analyse_faces(V)
        if not fa["admissible"]:
            continue
        ids = [short_id(l) for l in labels]
        pid = next((i for i in ids if i.startswith(("Delta", "X_"))), ids[0])
        R = sorted(set(map(tuple, V)) | {f["centre"] for f in fa["pentagons"] + fa["hexagons"]})
        rec = dict(id=pid, sources=labels, vertices=[list(v) for v in V],
                   squares=[[list(p) for p in c] for c in fa["squares"]],
                   pentagons=[dict(cycle=[list(p) for p in f["cycle"]], centre=list(f["centre"]))
                              for f in fa["pentagons"]],
                   hexagons=[dict(cycle=[list(p) for p in f["cycle"]], centre=list(f["centre"]))
                             for f in fa["hexagons"]],
                   R=[list(p) for p in R])
        polys.append(rec)
    log(f"{len(polys)} admissible vertex sets")

    prepared.write_text(json.dumps(polys, separators=(",", ":")))
    if A.prepare_only:
        return

    # X° and its equivalent
    xc = next(r for r in polys if r["id"] == "X_circ")
    sig = lambda r: (len(r["vertices"]), len(r["squares"]), len(r["pentagons"]), len(r["hexagons"]))
    for r in polys:
        if r is not xc and sig(r) == sig(xc):
            g = find_equivalence(xc["vertices"], r["vertices"])
            if g:
                xc["equivalent_to"] = dict(id=r["id"], matrix=g)
                log(f"X_circ is equivalent to {r['id']}")
                break
    assert "equivalent_to" in xc

    processed = 0
    for n, rec in enumerate(polys):
        if rec["id"] in completed:
            saved = completed[rec["id"]]
            assert saved["vertices"] == rec["vertices"] and saved["sources"] == rec["sources"], "checkpoint input mismatch"
            rec.update(saved)
            continue
        if A.limit is not None and processed >= A.limit:
            break
        if A.only and rec["id"] not in A.only:
            continue
        P = F.Polytope(rec)
        t0 = time.time()
        if not P.hexes:
            nd, Kb, forced = P.lattice({})
            rec["nodes"] = nodes_json(nd)
            rec["results"] = dict(nodes=len(nd), rank_K=len(Kb), forced=forced)
            if nd and Kb:
                ineq, eq = system(P, {})
                cert = farkas(ineq, eq)
                if cert is None:
                    psi = witness(ineq, eq, len(P.R))
                    F.check_witness(P, {}, psi)
                    rec["proj"] = dict(witness=psi)
                else:
                    rec["proj"] = dict(farkas=dict(y=cert[0], nu=cert[1]))
                rec["results"]["proj"] = cert is None
        else:
            lp = 0
            rec["kinds"] = {}
            for kind in ("A1", "A2", "mixed"):
                opts = [c for c in F.HEX_TILINGS if kind == "mixed" or c.startswith(kind)]
                hint = hints.get(rec["id"], {}).get(kind)
                if hint:
                    F.check_witness(P, hint["profile"], hint["psi"])
                    assert set(hint["profile"]) == {str(i) for i in range(len(P.hexes))}
                    pattern = tuple(hint["profile"][str(i)].split(":")[0] for i in range(len(P.hexes)))
                    assert F.pattern_filter(kind)(pattern, True)
                    res, k = dict(witness=dict(hint)), 0
                else:
                    res, k = search(P, F.pattern_filter(kind), opts, A.cap)
                lp += k
                if "witness" in res:
                    res["witness"]["nodes"] = nodes_json(P.nodes(res["witness"]["profile"]))
                rec["kinds"][kind] = res
            feas = {k: rec["kinds"][k].get("witness") for k in rec["kinds"]}
            hx = list(range(len(P.hexes)))
            out = dict(feasible={k: bool(v) for k, v in feas.items()})
            if any(feas.values()) and rec["id"] != "X_circ":
                fz = {k: (P.lattice(v["profile"])[2] if v else None) for k, v in feas.items()}
                good_w = [k for k in ("A1", "A2", "mixed") if feas[k] and not fz[k]]
                if good_w:
                    rec["corollary_b"] = dict(witness=feas[good_w[0]])
                else:
                    good = good_patterns(P)
                    res, k = search(P, F.pattern_filter(None, set(good)) if False else
                                    _reorder_filter(P, None, good), F.HEX_TILINGS, A.cap)
                    lp += k
                    if "witness" in res:
                        res["witness"]["nodes"] = nodes_json(P.nodes(res["witness"]["profile"]))
                    else:
                        res["good_patterns"] = [list(g) for g in good]
                    rec["corollary_b"] = res
                out["corB_any"] = "witness" in rec["corollary_b"]
                for kind in ("A1", "A2"):
                    out["corB_" + kind] = bool(feas[kind]) and not fz[kind]
                if feas["mixed"]:
                    if not fz["mixed"]:
                        rec["corollary_b_mixed"] = dict(witness=feas["mixed"])
                    else:
                        good = good_patterns(P, mixed_only=True)
                        res, k = search(P, _reorder_filter(P, "mixed", good), F.HEX_TILINGS, A.cap)
                        lp += k
                        if "witness" in res:
                            res["witness"]["nodes"] = nodes_json(P.nodes(res["witness"]["profile"]))
                        else:
                            res["good_patterns"] = [list(g) for g in good]
                        rec["corollary_b_mixed"] = res
                    out["corB_mixed"] = "witness" in rec["corollary_b_mixed"]
                if "witness" in rec["corollary_b"]:
                    w = rec["corollary_b"]["witness"]
                    w.setdefault("nodes", nodes_json(P.nodes(w["profile"])))
            rec["results"] = out
            log(f"{rec['id']}: {len(P.hexes)} hexagons, {out}, {lp} LPs [{time.time() - t0:.1f}s]")
        processed += 1
        completed[rec["id"]] = rec
        temp = checkpoint.with_suffix(".tmp")
        temp.write_text(json.dumps(completed, separators=(",", ":")))
        temp.replace(checkpoint)
        if rec["id"] in old:
            rec["stored"] = old[rec["id"]]
        if (n + 1) % 50 == 0:
            log(f"{n + 1}/{len(polys)}")
    if len(completed) != len(polys):
        log(f"checkpoint retained: {len(completed)}/{len(polys)}; no complete data file written")
        return
    for rec in polys:
        result = rec["results"]
        if not rec["hexagons"]:
            rec["stored"] = dict(nodes=result["nodes"], rank_K=result["rank_K"],
                                 forced=len(result["forced"]), no_forced=not result["forced"],
                                 proj=result.get("proj"))
        else:
            rec["stored"] = dict(feasible=result["feasible"], corB={
                kind: bool(result.get("corB_" + kind, False)) for kind in ("any", "A1", "A2", "mixed")})
    polys.sort(key=_sort_key)
    doc = dict(description="Admissible polytopes of the list of Section 11.4, their two-faces, and (Proj) witnesses "
                           "and Farkas certificates; see forced_nodes_check.py and the README.",
               polytopes=polys)
    with gzip.open(A.out, "wt") as fo:
        json.dump(doc, fo, separators=(",", ":"))
    log(f"wrote {A.out}")


def _reorder_filter(P, kind, good):
    """The pattern filter for a set of good full patterns (indexed by hexagon), in the search order of order_hexes."""
    order = order_hexes(P)
    return F.pattern_filter(kind, {tuple(g[h] for h in order) for g in good})


def _sort_key(r):
    i = r["id"]
    if i.startswith("Delta_") or i == "X_circ":
        return (0, i)
    kind, n = i.split(":")
    return (1 if kind == "framework" else 2, int(n))


if __name__ == "__main__":
    main()
