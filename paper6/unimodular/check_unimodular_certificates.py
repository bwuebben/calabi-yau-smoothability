#!/usr/bin/env python3
"""Integer check of the unimodular-height certificates of Theorem F (Appendix A).

Run from the repository root:

    python3 paper6/unimodular/check_unimodular_certificates.py            # all files
    python3 paper6/unimodular/check_unimodular_certificates.py --procs 8
    python3 paper6/unimodular/check_unimodular_certificates.py --quick 20 # first 20 of each file
    python3 paper6/unimodular/check_unimodular_certificates.py --table    # Table of admissible counts

Requirements: Python 3 standard library only.

Data.  paper6/unimodular/certificates/heights_vNN.jsonl.gz holds one JSON object per admissible reflexive
polytope Delta with NN vertices:
    row        0-based row index of Delta in the Kreuzer-Skarke file polytopes-4d-NN-vertices.parquet
               (pinned revision of manifests/ks_polytopes_4d_sha256.tsv);
    V          the vertices of Delta, as listed in that row;
    two_faces  the numbers of two-faces of Delta that are unimodular triangles, unit squares, dP7 pentagons
               and dP6 hexagons;
    pts        the nonzero lattice points m of the polar polytope Delta° = {m : <m,v> >= -1 for v in V};
    h          the integer height h(m) of each point of pts (same order): the certificate;
    quadratic_forms_tried   how many quadratic forms were tried before h was found (information only).

The certificate h determines a candidate triangulation: in each facet F_v = {<m,v> = -1} of Delta°, the lower
faces of the convex hull of the lifted points (m, h(m)), m in F_v.  This program computes these lower faces
by gift wrapping (a proposal only, see below), and then checks in exact integer arithmetic, with g = h + lam
on the boundary points (lam >= 0 an integer computed here), extended linearly on the cones over the candidate
tetrahedra, and hcheck = g + phi, where phi = 1 on the boundary of Delta°:

  (i)   every recorded point is a nonzero lattice point of the boundary of Delta° (<m,v> >= -1 for all
        vertices v, with equality for some v), and every height is an integer;
  (ii)  every candidate tetrahedron lies in a facet hyperplane of Delta° and has determinant +-1
        (determinants by fraction-free Bareiss elimination);
  (iii) every triangle of a candidate tetrahedron lies in exactly two candidate tetrahedra, and their fourth
        vertices lie on opposite sides of the hyperplane spanned by the triangle and the origin;
  (iv)  each of three fixed rays in general position lies in the interior of the cone over exactly one
        candidate tetrahedron (five rays are listed; a ray found on a wall is skipped, and at least three
        must be usable);
  (v)   across every wall, the bends of g and of hcheck are positive integers.

By the lemma of Appendix A on (i)-(v), these conditions imply that the candidate tetrahedra form a unimodular
triangulation T of the boundary of Delta° refining its faces, and that hcheck is an integral height on the
fan over T with hcheck and hcheck - phi strictly convex: the conclusion of Theorem F for Delta.

How the tetrahedra are proposed plays no role in this conclusion; only (i)-(v) are used.  The proposal is an
exact gift-wrapping computation of the lower hull in each facet: a first lower cell is found by tilting a
supporting hyperplane from the lowest lifted point until it touches four points, and each further cell is
found by rotating the hyperplane of a known cell about one of its triangles until it meets the next lifted
point.  If the lower hull has a cell that is not a tetrahedron, the certificate is rejected.
"""
import argparse
import glob
import gzip
import itertools
import json
import os
import sys
import time
from fractions import Fraction
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
CERT_DIR = os.path.join(HERE, "certificates")
REPO = os.path.dirname(os.path.dirname(HERE))


class Reject(Exception):
    pass


def need(cond, msg):
    if not cond:
        raise Reject(msg)


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3]


# ------------------------------------------------------------------ exact integer linear algebra

def det4(M):
    """Determinant of a 4x4 integer matrix by Bareiss fraction-free elimination."""
    a = [list(r) for r in M]
    n, sign, prev = 4, 1, 1
    for k in range(n - 1):
        if a[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if a[i][k] != 0), None)
            if sw is None:
                return 0
            a[k], a[sw] = a[sw], a[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] = (a[i][j] * a[k][k] - a[i][k] * a[k][j]) // prev
        prev = a[k][k]
    return sign * a[n - 1][n - 1]


def det3(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def inverse_unimodular(M, D):
    """Integer inverse of a 4x4 integer matrix M (rows) with det M = D = +-1, by Cramer's rule (adjugate).
    Returns Minv with Minv[i][j] such that x = sum_j c_j M[j] for c_j = sum_i x_i Minv[i][j]."""
    inv = [[0] * 4 for _ in range(4)]
    for i in range(4):
        for j in range(4):
            rows = [[M[r][s] for s in range(4) if s != j] for r in range(4) if r != i]
            cof = (-1) ** (i + j) * det3(*rows)
            inv[j][i] = cof * D                      # (adj M)_{ji} / D, and 1/D = D
    for i in range(4):                               # exactness: M * inv = identity
        for j in range(4):
            need(sum(M[i][k] * inv[k][j] for k in range(4)) == (1 if i == j else 0), "(ii) inverse")
    return inv


def coords(inv, x):
    """Coefficients c with x = sum c_j M_j."""
    return [x[0] * inv[0][j] + x[1] * inv[1][j] + x[2] * inv[2][j] + x[3] * inv[3][j] for j in range(4)]


def int_kernel(rows):
    """A nonzero integer vector orthogonal to the given integer 4-vectors (rank < 4), or None."""
    A = [[Fraction(x) for x in r] for r in rows]
    piv, r = [], 0
    for c in range(4):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        A[r] = [x / A[r][c] for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c)
        r += 1
    free = [c for c in range(4) if c not in piv]
    if not free:
        return []
    out = []
    for fc in free:
        w = [Fraction(0)] * 4
        w[fc] = Fraction(1)
        for i, pc in enumerate(piv):
            w[pc] = -A[i][fc]
        den = 1
        for x in w:
            den = den * x.denominator // _gcd(den, x.denominator)
        out.append([int(x * den) for x in w])
    return out


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return abs(a)


def rank(rows):
    return 4 - len(int_kernel(rows)) if rows else 0


# ------------------------------------------------------------------ proposal: lower hull by gift wrapping

def first_cell(P, H):
    """Four indices spanning a lower cell of the lifted points (P[i], H[i]) of one facet, found by tilting a
    supporting hyperplane.  P: points of the facet (all with <p,v> = -1).  Linear functionals a on Z^4 restrict
    to affine functions on the facet hyperplane."""
    n = len(P)
    hmin = min(H)
    contact = [i for i in range(n) if H[i] == hmin]
    # The supporting functional a starts as the constant hmin on the facet and is tilted by functionals w
    # vanishing on the current contact set; we carry the slacks H[p] - a(p) >= 0, zero exactly on the contacts.
    slack = [Fraction(H[i] - hmin) for i in range(n)]
    while rank([P[i] for i in contact]) < 4:
        ker = int_kernel([P[i] for i in contact])
        w = None
        for cand in ker:
            vals = [dot(cand, P[i]) for i in range(n)]
            if any(x > 0 for x in vals):
                w, wv = cand, vals
                break
            if any(x < 0 for x in vals):
                w, wv = [-x for x in cand], [-x for x in vals]
                break
        need(w is not None, "facet degenerate")
        s = min(slack[i] / wv[i] for i in range(n) if wv[i] > 0)
        slack = [slack[i] - s * wv[i] for i in range(n)]
        need(min(slack) >= 0, "tilting")
        contact = [i for i in range(n) if slack[i] == 0]
    need(len(contact) == 4, "(ii) a lower face of the lifted points is not a tetrahedron")
    return tuple(contact)


def lower_cells(P, H):
    """All lower cells (as 4-tuples of local indices) of the lifted points of one facet, by gift wrapping."""
    n = len(P)
    if n == 4:
        return [(0, 1, 2, 3)]
    start = tuple(sorted(first_cell(P, H)))
    cells = {start}
    tri = {}
    queue = [start]
    while queue:
        t = queue.pop()
        M = [P[i] for i in t]
        D = det4(M)
        need(abs(D) == 1, "(ii) candidate tetrahedron not unimodular")
        inv = inverse_unimodular(M, D)
        lam = [coords(inv, P[i]) for i in range(n)]     # barycentric coordinates w.r.t. t (integers)
        ht = [H[i] for i in t]
        for k in range(4):                               # triangle opposite t[k]
            tau = tuple(t[j] for j in range(4) if j != k)
            lst = tri.setdefault(tau, [])
            if t not in lst:
                lst.append(t)
            if len(lst) >= 2:
                continue
            best, tie = None, False
            for i in range(n):
                ld = lam[i][k]
                if ld >= 0:
                    continue
                m = H[i] - (lam[i][0] * ht[0] + lam[i][1] * ht[1] + lam[i][2] * ht[2] + lam[i][3] * ht[3])
                # rotation parameter m / (-ld); minimal ratio wins (compare exactly)
                if best is None or m * best[2] < best[1] * (-ld):
                    best, tie = (i, m, -ld), False
                elif m * best[2] == best[1] * (-ld):
                    tie = True
            need(not tie, "(ii) a lower face of the lifted points is not a tetrahedron")
            if best is None:
                continue                                 # triangle on the boundary of the facet
            new = tuple(sorted(tau + (best[0],)))
            lst.append(new)
            if new not in cells:
                cells.add(new)
                queue.append(new)
    return sorted(cells)


# ------------------------------------------------------------------ the checks (i)-(v)

RAYS = [(1000003, 999983, -1000037, 1000039), (-7919, 104729, 1299709, -15485863), (3, -5, 7, 11),
        (-1000033, 1000081, 999979, -999961), (104723, -104717, 104711, 104707)]


def check(rec):
    V = [tuple(v) for v in rec["V"]]
    pts = [tuple(p) for p in rec["pts"]]
    h = rec["h"]
    # (i)
    need(all(len(v) == 4 for v in V) and all(len(p) == 4 for p in pts),
         "(i) coordinates must be four-dimensional")
    need(len(h) == len(pts) and len(set(pts)) == len(pts), "(i) point list")
    need(all(type(x) is int for x in h), "(i) heights not integers")
    need(all(type(x) is int for p in pts for x in p) and all(type(x) is int for v in V for x in v),
         "(i) coordinates not integers")
    for p in pts:
        vals = [dot(p, v) for v in V]
        need(any(p) and min(vals) == -1, "(i) not a nonzero boundary lattice point of the polar polytope")
    # proposal, facet by facet
    tets = set()
    for v in V:
        idx = [i for i, p in enumerate(pts) if dot(p, v) == -1]
        need(len(idx) >= 4, "facet with fewer than four points")
        for c in lower_cells([pts[i] for i in idx], [h[i] for i in idx]):
            tets.add(tuple(sorted(idx[j] for j in c)))
    tets = sorted(tets)
    # (ii)
    inv = {}
    for t in tets:
        need(any(all(dot(pts[i], v) == -1 for i in t) for v in V), "(ii) tetrahedron not in a facet")
        M = [pts[i] for i in t]
        D = det4(M)
        need(abs(D) == 1, "(ii) tetrahedron not unimodular")
        inv[t] = inverse_unimodular(M, D)
    # (iii)
    walls = {}
    for t in tets:
        for tau in itertools.combinations(t, 3):
            walls.setdefault(tau, []).append((t, next(i for i in t if i not in tau)))
    for tau, ap in walls.items():
        need(len(ap) == 2, "(iii) triangle not in exactly two tetrahedra")
        s1 = det4([pts[i] for i in tau] + [pts[ap[0][1]]])
        s2 = det4([pts[i] for i in tau] + [pts[ap[1][1]]])
        need(s1 * s2 < 0, "(iii) the two tetrahedra lie on the same side")
    # (iv)
    good = 0
    for r in RAYS:
        hits, onwall = 0, False
        for t in tets:
            c = coords(inv[t], r)
            mn = min(c)
            if mn > 0:
                hits += 1
            elif mn == 0:
                onwall = True
                break
        if onwall:
            continue
        need(hits == 1, "(iv) covering degree is not one")
        good += 1
    need(good >= 3, "(iv) fewer than three usable rays")
    # (v) bends: the wall tau is a facet of t = tau + {a} and of t' = tau + {b}; write pts[b] = sum c_j t_j.
    data = []
    for tau, ((t, a), (_, b)) in walls.items():
        c = coords(inv[t], pts[b])
        m = h[b] - sum(ci * h[u] for ci, u in zip(c, t))
        data.append((m, 1 - sum(c)))                     # bend of h, and bend of phi (0 inside a facet)
    for m, e in data:
        need(e >= 0, "(v) phi not convex")
        if e == 0:
            need(m > 0, "(v) h not strictly convex inside a facet")
    lam = max([0] + [(-m) // e + 1 for m, e in data if e > 0 and m <= 0])
    for m, e in data:
        need(m + lam * e > 0, "(v) g not strictly convex")
        need(m + lam * e + e > 0, "(v) hcheck not strictly convex")
    return {"row": rec["row"], "ok": True, "tets": len(tets), "walls": len(walls), "lam": lam}


def safe(rec):
    try:
        return check(rec)
    except Reject as e:
        return {"row": rec.get("row"), "ok": False, "err": str(e)}
    except Exception as e:                               # malformed data counts as rejection
        return {"row": rec.get("row"), "ok": False, "err": f"{type(e).__name__}: {e}"}


def records(path, quick=None):
    with gzip.open(path, "rt") as f:
        for k, line in enumerate(f):
            if quick is not None and k >= quick:
                return
            yield json.loads(line)


def run(files, procs, quick):
    t0 = time.time()
    total_n = total_ok = total_t = 0
    failed = []
    with Pool(procs) as pool:
        for path in files:
            n = ok = tt = 0
            t1 = time.time()
            for res in pool.imap(safe, records(path, quick), chunksize=4):
                n += 1
                if res["ok"]:
                    ok += 1
                    tt += res["tets"]
                else:
                    failed.append((os.path.basename(path), res))
                    print("REJECTED", os.path.basename(path), res, flush=True)
                if n % 1000 == 0:
                    print(f"  [{time.time() - t0:7.0f}s] {os.path.basename(path)}: {n} checked", flush=True)
            print(f"{os.path.basename(path)}: checked {n}, passing {ok}, tetrahedra {tt} "
                  f"[{time.time() - t1:.0f}s]", flush=True)
            total_n += n
            total_ok += ok
            total_t += tt
    print(f"TOTAL: checked {total_n}, passing {total_ok}, rejected {len(failed)}, tetrahedra {total_t} "
          f"[{time.time() - t0:.0f}s]", flush=True)
    return total_n, total_ok, total_t


def table(files):
    """Admissible polytopes by number of vertices, from the certificate files; the numbers of polytopes per file
    are read from the required bundled provenance/classification_inventory.json."""
    rows = {}
    ref = os.path.join(os.path.dirname(HERE), "provenance", "classification_inventory.json")
    if os.path.exists(ref):
        for f in json.load(open(ref))["files"]:
            rows[f["vertices"]] = f["rows"]
    if not rows:
        raise ValueError("classification inventory missing or empty; counts unavailable")
    print(f"{'vertices':>8} {'polytopes':>12} {'admissible':>10} {'with hexagon':>12} {'no singular 2-face':>18}")
    T = [0, 0, 0, 0]
    for path in files:
        nv = int(os.path.basename(path)[len("heights_v"):len("heights_v") + 2])
        adm = hexa = smooth = 0
        for r in records(path):
            assert len(r["V"]) == nv
            tf = r["two_faces"]
            adm += 1
            hexa += tf["hexagon"] > 0
            smooth += tf["square"] == tf["pentagon"] == tf["hexagon"] == 0
        if nv not in rows:
            raise ValueError(f"classification inventory lacks vertex count {nv}")
        k = rows[nv]
        print(f"{nv:>8} {k:>12,} {adm:>10,} {hexa:>12,} {smooth:>18,}")
        for i, x in enumerate((k, adm, hexa, smooth)):
            T[i] += x
    print(f"{'total':>8} {T[0]:>12,} {T[1]:>10,} {T[2]:>12,} {T[3]:>18,}")
    print(f"without hexagon: {T[1] - T[2]:,}")
    return T


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("files", nargs="*", help="certificate files (default: all in paper6/unimodular/certificates)")
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--quick", type=int, default=None, metavar="K", help="check only the first K of each file")
    ap.add_argument("--table", action="store_true", help="print the counts of admissible polytopes and stop")
    A = ap.parse_args()
    files = A.files or sorted(glob.glob(os.path.join(CERT_DIR, "heights_v*.jsonl.gz")))
    if not files:
        ap.error("no certificate files found")
    if A.table:
        T = table(files)
        if T[1] == 0:
            ap.error("no certificate records found")
        if not A.files:
            assert T[1:] == [51827, 11424, 741], T
            print("agrees with Appendix A: 51,827 admissible, 11,424 with a hexagon, 741 with no singular 2-face")
        return
    n, ok, tt = run(files, A.procs, A.quick)
    if n == 0:
        ap.error("no certificate records checked")
    if not A.files and A.quick is None:
        assert (n, ok, tt) == (51827, 51827, 12001172), (n, ok, tt)
        print("agrees with Appendix A: all 51,827 certificates pass, 12,001,172 tetrahedra")
    sys.exit(0 if n == ok else 1)


if __name__ == "__main__":
    main()
