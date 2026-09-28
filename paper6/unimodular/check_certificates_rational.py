#!/usr/bin/env python3
"""Second check of the unimodular-height certificates of Theorem F (Appendix A), in exact rational arithmetic.

Run from the repository root with SageMath's Python:

    sage -python paper6/unimodular/check_certificates_rational.py               # all 51,827 certificates
    sage -python paper6/unimodular/check_certificates_rational.py --procs 8
    sage -python paper6/unimodular/check_certificates_rational.py --sample 20   # 20 per file, plus the largest

Requirements: SageMath (polyhedra with the PPL backend, exact integer and rational matrices) and numpy.
This program shares no code with check_unimodular_certificates.py; it proves the conclusion of Theorem F from
the same heights by a different set of sufficient conditions.  For a certificate (V, pts, h), see the
description of the data in check_unimodular_certificates.py:

  (i')  Delta = conv(V) has exactly the vertices V; its polar Delta° is a lattice polytope; the nonzero lattice
        points of Delta°, enumerated independently (integral points of the polyhedron), are exactly the
        recorded points, and all lie on the boundary of Delta°;
  (ii') in each facet F_v = {<m,v> = -1} of Delta°, every lower facet of the convex hull of the lifted points
        (m, h(m)), m in F_v, computed exactly, contains exactly four lifted points, and these span a cone of
        determinant +-1.  So each facet is triangulated by the regular subdivision induced by h, and since the
        heights are global, the triangulations of two facets restrict to the same regular subdivision of
        their common face: the tetrahedra form a unimodular triangulation T of the boundary of Delta°
        refining its faces;
  (v')  for each tetrahedron t, lying in F_v, and each recorded point m not in t: with l_t the linear function
        equal to h on t, the bend m_t(m) = h(m) - l_t(m) and the bend e_t(m) = 1 + <m,v> of phi satisfy
        e >= 0, m_t(m) > 0 where e = 0; for the least integer lam >= 0 with m_t + lam e > 0 everywhere,
        g = h + lam phi and hcheck = g + phi exceed the linear functions agreeing with them on the cone over
        t at every recorded point outside it.  On the complete fan over T this is strict convexity of g and
        of hcheck; hcheck is integral since T is unimodular.
As a consistency check the number of tetrahedra is compared with the normalised volume of Delta°.
"""
import argparse
import glob
import gzip
import json
import os
import random
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CERT_DIR = os.path.join(HERE, "certificates")


def check(rec):
    from sage.all import Polyhedron, matrix, ZZ, QQ, vector
    assert all(len(v) == 4 and all(type(x) is int for x in v) for v in rec["V"]), "(i') vertices must be four-dimensional integer vectors"
    assert all(len(p) == 4 and all(type(x) is int for x in p) for p in rec["pts"]), "(i') points must be four-dimensional integer vectors"
    assert all(type(x) is int for x in rec["h"]), "(i') heights not integers"
    V = [tuple(int(x) for x in v) for v in rec["V"]]
    pts = [tuple(int(x) for x in p) for p in rec["pts"]]
    h = [int(x) for x in rec["h"]]
    H = dict(zip(pts, h))
    assert len(H) == len(pts) == len(h), "(i') point list"
    # (i')
    P = Polyhedron(vertices=V, backend="ppl")
    assert P.dim() == 4 and {tuple(int(x) for x in v) for v in P.vertices()} == set(V), "(i') V is not the vertex set"
    Q = P.polar()                                              # {m : <m,x> >= -1 on P}
    assert all(x in ZZ for v in Q.vertices() for x in v), "(i') polar polytope not a lattice polytope"
    lat = [tuple(int(x) for x in p) for p in Q.integral_points()]
    nz = [p for p in lat if any(p)]
    assert len(lat) == len(nz) + 1, "(i') origin"
    assert set(nz) == set(pts), "(i') recorded points differ from the nonzero lattice points of the polar polytope"
    Varr = np.array(V, dtype=np.int64)
    for p in nz:
        assert int((Varr @ np.array(p, dtype=np.int64)).min()) == -1, "(i') point not on the boundary"
    # (ii') the facets of Delta° are the F_v, v a vertex of Delta
    tets = {}
    for v in V:
        Fp = [p for p in pts if sum(a * b for a, b in zip(p, v)) == -1]
        assert len(Fp) >= 4, "(ii') facet with fewer than four points"
        if len(Fp) == 4:
            assert abs(matrix(ZZ, Fp).det()) == 1, "(ii') not unimodular"
            cells = [tuple(sorted(Fp))]
        else:
            lifted = [tuple(p) + (H[p],) for p in Fp]
            L = Polyhedron(vertices=lifted, backend="ppl")
            assert L.dim() == 4, "(ii') lifted facet degenerate"
            cells = []
            for ie in L.Hrepresentation():
                if ie.is_equation():
                    continue
                if list(ie.A())[4] <= 0:                       # upper or vertical facet
                    continue
                on = [q for q in lifted if ie.eval(vector(QQ, q)) == 0]
                assert len(on) == 4, "(ii') a lower cell is not a tetrahedron"
                T = [q[:4] for q in on]
                assert abs(matrix(ZZ, T).det()) == 1, "(ii') not unimodular"
                cells.append(tuple(sorted(T)))
        for t in cells:
            assert t not in tets, "(ii') tetrahedron in two facets"
            tets[t] = v
    volQ = Q.volume() * 24                                     # normalised volume of Delta°
    assert volQ == len(tets), "(ii') number of tetrahedra differs from the normalised volume"
    # (v')
    Parr = np.array(pts, dtype=object)
    harr = np.array(h, dtype=object)
    idx = {p: i for i, p in enumerate(pts)}
    data = []
    for t, v in tets.items():
        Mi = matrix(ZZ, t).inverse()
        assert all(x in ZZ for x in Mi.list())
        a = Mi * vector(ZZ, [H[p] for p in t])                  # l_t(x) = <a, x>
        a = np.array([int(x) for x in a], dtype=object)
        m = harr - Parr @ a
        e = 1 + Parr @ np.array(v, dtype=object)
        mask = np.ones(len(pts), dtype=bool)
        for p in t:
            mask[idx[p]] = False
        m, e = m[mask], e[mask]
        assert all(x >= 0 for x in e), "(v') phi not convex"
        z = e == 0
        assert all(x > 0 for x in m[z]), "(v') h not strictly convex inside a facet"
        data.append((m[~z], e[~z]))
    lam = 0
    for m, e in data:
        for mm, ee in zip(m, e):
            if mm + lam * ee <= 0:
                lam = max(lam, (-mm) // ee + 1)
    for m, e in data:
        assert all(mm + lam * ee > 0 for mm, ee in zip(m, e)), "(v') g not strictly convex"
        assert all(mm + (lam + 1) * ee > 0 for mm, ee in zip(m, e)), "(v') hcheck not strictly convex"
    result = {"row": rec["row"], "tets": len(tets), "pts": len(pts), "lam": int(lam)}
    if os.environ.get("P6_CERT_ARCHIVE"):
        # Point indices refer to the unchanged certificate point order.
        result["certificate"] = rec
        result["simplices"] = [[idx[p] for p in t] for t in tets]
        result["simplex_facets"] = [V.index(v) for v in tets.values()]
    return result


def safe(job):
    name, rec = job
    try:
        r = check(rec)
        r["ok"] = True
    except Exception as e:                                      # any failure is a rejection
        r = {"row": rec.get("row"), "ok": False, "err": f"{type(e).__name__}: {str(e)[:200]}"}
    r["file"] = name
    return r


def jobs(files, sample, seed):
    rng = random.Random(seed)
    for f in files:
        name = os.path.basename(f)
        with gzip.open(f, "rt") as fh:
            recs = [json.loads(l) for l in fh]
        if sample is not None and recs:
            pick = rng.sample(range(len(recs)), min(sample, len(recs)))
            big = max(range(len(recs)), key=lambda i: len(recs[i]["pts"]))
            if big not in pick:
                pick.append(big)
            recs = [recs[i] for i in sorted(pick)]
        for r in recs:
            yield (name, r)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("files", nargs="*", help="certificate files (default: all in paper6/unimodular/certificates)")
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--sample", type=int, default=None, metavar="K",
                    help="check K random certificates per file and the largest one of each file")
    ap.add_argument("--seed", type=int, default=1)
    A = ap.parse_args()
    files = A.files or sorted(glob.glob(os.path.join(CERT_DIR, "heights_v*.jsonl.gz")))
    if not files:
        ap.error("no certificate files found")
    t0 = time.time()
    n = ok = tt = 0
    per = {}
    from contextlib import nullcontext
    archive_path = os.environ.get("P6_CERT_ARCHIVE")
    archive_context = gzip.open(archive_path + ".partial", "wt", compresslevel=1) if archive_path else nullcontext(None)
    with archive_context as archive, Pool(A.procs) as pool:
        for r in pool.imap_unordered(safe, jobs(files, A.sample, A.seed), chunksize=2):
            if archive is not None:
                archive.write(json.dumps(r, separators=(",", ":")) + "\n")
                if n % 100 == 0:
                    archive.flush()
            n += 1
            f = per.setdefault(r["file"], [0, 0, 0])
            f[0] += 1
            if r["ok"]:
                ok += 1
                tt += r["tets"]
                f[1] += 1
                f[2] += r["tets"]
            else:
                print("REJECTED", r, flush=True)
            if n % 100 == 0:
                print(f"[{time.time() - t0:6.0f}s] {n} checked, {ok} pass", flush=True)
    if archive_path and n == ok:
        os.replace(archive_path + ".partial", archive_path)
    for name in sorted(per):
        print(f"{name}: checked {per[name][0]}, passing {per[name][1]}, tetrahedra {per[name][2]}")
    print(f"TOTAL: checked {n}, passing {ok}, rejected {n - ok}, tetrahedra {tt} [{time.time() - t0:.0f}s]",
          flush=True)
    if n == 0:
        ap.error("no certificate records checked")
    if not A.files and A.sample is None:
        assert (n, ok, tt) == (51827, 51827, 12001172), (n, ok, tt)
        print("agrees with Appendix A: all 51,827 certificates pass, 12,001,172 tetrahedra")
    sys.exit(0 if n == ok else 1)


if __name__ == "__main__":
    main()
