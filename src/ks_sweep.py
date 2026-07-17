#!/usr/bin/env python3
"""
Kreuzer-Skarke sweep (paper Sec. 6, Q1 + Q3): scan the KS database of
reflexive 4-polytopes for non-smoothable 2-faces.

Data: https://huggingface.co/datasets/calabi-yau-data/polytopes-4d
(parquet, one file per vertex count, columns: vertices [normal form],
vertex_count, facet_count, point_count, dual_point_count, h11, h12,
euler_characteristic).  Download files into data/ks/ first.

Per polytope Delta we compute all 2-faces with their induced-lattice
polygons and classify each:
  * smooth (standard triangle),
  * unit-edge rigid / def-only / smoothable (the census trichotomy; for
    rigid or def-only faces the generic anticanonical X is NON-SMOOTHABLE
    by the paper's Corollary — theorem-grade),
  * non-unit-edge and with no decomposition into unit segments + standard
    triangles: recorded as CANDIDATE non-smoothable only (the germ is
    non-isolated; the correct criterion there is Filip's 0-mutability,
    arXiv:2504.04486, which we have not implemented) — counted separately,
    never claimed.

For each unit-edge non-smoothable face we decide which of the 87 census
classes it is, by exact GL2(Z) matching: a lattice equivalence maps edges to
edges, so it is determined by the images of two independent edge vectors;
trying all ordered image pairs decides equivalence.  Faces equivalent to no
census class are non-smoothable local models OUTSIDE the |.|<=2 census box
and are recorded.

Cross-validation per file: on a random sample of rows our Batyrev Hodge
numbers (src/hodge_numbers.py) must equal the dataset's h11/h12 columns.

Also hunted: a polytope, all edges unit, whose unique non-smooth 2-face is
def-only with dual edge of lattice length 1 -> a compact CY threefold with a
SINGLE deformable-but-never-smoothable point (the paper's sharpened Q2).

Run:  ./venv/bin/python src/ks_sweep.py data/ks/polytopes-4d-06-vertices.parquet
      ./venv/bin/python src/ks_sweep.py data/ks/*.parquet --procs 8
"""
from fractions import Fraction
from itertools import combinations
from multiprocessing import Pool
import argparse
import json
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from toric_census import rigid, smoothing_components                       # noqa: E402
from batyrev_global import (dot, vsub, vgcd, hyperplane_normal,            # noqa: E402
                            int_kernel, solve_int_coords)
from plant_search import nonsmoothable_classes                             # noqa: E402
from hodge_numbers import hodge_numbers                                    # noqa: E402

from math import atan2, gcd                                                # noqa: E402


# ---------------- exact facet/face machinery (n up to ~14 vertices) -------
def facets_reflexive(V):
    """Facets of reflexive conv(V): [(u primitive, idx frozenset)].
    Levels are asserted to be exactly 1 (the dataset is reflexive)."""
    out = {}
    for S in combinations(range(len(V)), 4):
        u = hyperplane_normal(*[V[i] for i in S])
        if u == (0, 0, 0, 0):
            continue
        c = dot(u, V[S[0]])
        vals = [dot(u, v) for v in V]
        if all(x <= c for x in vals):
            pass
        elif all(x >= c for x in vals):
            u = tuple(-x for x in u); c = -c; vals = [-x for x in vals]
        else:
            continue
        g = vgcd(u)
        assert c == g, "polytope not reflexive?"
        u = tuple(x // g for x in u)
        if u not in out:
            out[u] = frozenset(i for i, x in enumerate(vals) if x == g)
    return list(out.items())


def polygon_of_face(V, idx, u1, u2):
    """Induced-lattice polygon of the 2-face: (primitive edge vectors ccw,
    edge lattice lengths)."""
    pts = [V[i] for i in idx]
    ker = int_kernel([list(u1), list(u2)])
    p0 = pts[0]
    coords = [solve_int_coords(ker, vsub(p, p0)) for p in pts]
    cx = Fraction(sum(c[0] for c in coords), len(coords))
    cy = Fraction(sum(c[1] for c in coords), len(coords))
    coords.sort(key=lambda c: atan2(c[1] - cy, c[0] - cx))
    k = len(coords)
    evs, lens = [], []
    for i in range(k):
        d = vsub(coords[(i + 1) % k], coords[i])
        g = vgcd(d)
        evs.append((d[0] // g, d[1] // g))
        lens.append(g)
    return evs, lens


def two_faces_of(V, facs):
    found = {}
    for (a, b) in combinations(range(len(facs)), 2):
        I = facs[a][1] & facs[b][1]
        if len(I) >= 3:
            key = I
            if key not in found:
                found[key] = (a, b)
    # keep only genuine 2-faces: intersection of exactly 2 facets is 2-dim
    # for a 4-polytope whenever it has >= 3 affinely independent points; the
    # >= 3 shared vertices of two distinct facets of a 4-polytope span the
    # 2-face (ridges of facets are 2-faces).
    return found


# ================== FAST ENGINE (integer-exact, numpy-batched) ==================
# Reliability contract:
#   * all arithmetic is int64; a load-time assert (max|coord| <= 5000) plus the
#     Hadamard bound (normal entries <= (sqrt3*2c)^3, support values times 4c)
#     keeps every intermediate < 2^62 — numpy wraps silently on overflow, so
#     the assert is the load-bearing guard;
#   * no floats anywhere: integer cofactor determinants (not np.linalg.det),
#     np.gcd, and an exact integer half-plane comparator for the ccw sort;
#   * the pure-python engine above stays as the reference: --selftest compares
#     both engines on random rows of every file before sweeping it.

COORD_BOUND = 5000

def ccw_sorted_int(coords):
    """Exact ccw sort of distinct integer points around their centroid
    (quadrant + cross product; no floats)."""
    from functools import cmp_to_key
    k = len(coords)
    cx = sum(c[0] for c in coords)
    cy = sum(c[1] for c in coords)
    def vec(p):
        return (k * p[0] - cx, k * p[1] - cy)
    def half(v):
        return 0 if (v[1] > 0 or (v[1] == 0 and v[0] > 0)) else 1
    def cmp(p, q):
        vp, vq = vec(p), vec(q)
        hp, hq = half(vp), half(vq)
        if hp != hq:
            return hp - hq
        cr = vp[0] * vq[1] - vp[1] * vq[0]
        return -1 if cr > 0 else (1 if cr < 0 else 0)
    return sorted(coords, key=cmp_to_key(cmp))


_KERNEL_CACHE = {}

def _face_kernel(u1, u2):
    """Saturated kernel basis of [u1, u2] + a pivot 2x2 minor, cached: facet
    normal pairs recur heavily across polytopes."""
    key = (u1, u2)
    hit = _KERNEL_CACHE.get(key)
    if hit is not None:
        return hit
    k1, k2 = int_kernel([list(u1), list(u2)])
    piv = None
    for i in range(4):
        for j in range(i + 1, 4):
            d = k1[i] * k2[j] - k1[j] * k2[i]
            if d:
                piv = (i, j, d)
                break
        if piv:
            break
    res = (k1, k2, *piv)
    _KERNEL_CACHE[key] = res
    return res


def polygon_of_face_int(V, idx, u1, u2):
    """Integer-only version of polygon_of_face (no Fractions): induced-lattice
    coordinates via Cramer on a pivot 2x2 minor of the saturated kernel basis;
    exact-division asserts replace the rational solve."""
    pts = [V[i] for i in idx]
    k1, k2, i, j, d = _face_kernel(u1, u2)
    k1i, k1j, k2i, k2j = k1[i], k1[j], k2[i], k2[j]
    p0 = pts[0]
    p0i, p0j = p0[i], p0[j]
    coords = []
    for p in pts:
        ti, tj = p[i] - p0i, p[j] - p0j
        xn = ti * k2j - tj * k2i
        yn = k1i * tj - k1j * ti
        x, y = xn // d, yn // d
        assert x * d == xn and y * d == yn, "non-integral face coordinate"
        # (no 4-row re-verification needed: a face point lies on both facets,
        # hence p - p0 is in the saturated kernel span by construction; the
        # per-file engine selftest guards the implementation end to end)
        coords.append((x, y))
    coords = ccw_sorted_int(coords)
    k = len(coords)
    evs, lens = [], []
    for a in range(k):
        dxy = vsub(coords[(a + 1) % k], coords[a])
        g = vgcd(dxy)
        evs.append((dxy[0] // g, dxy[1] // g))
        lens.append(g)
    return evs, lens


def facets_np(verts, subsets):
    """Batched facet enumeration.  verts: (B, n, 4) int64; subsets: (C, 4)
    index array.  Returns per-polytope list of (u tuple, idx frozenset).
    Raises AssertionError if any supporting functional is not at level 1."""
    import numpy as np
    B, n, _ = verts.shape
    pts = verts[:, subsets]                     # (B, C, 4, 4)
    d = pts[:, :, 1:] - pts[:, :, :1]           # (B, C, 3, 4)
    cols = [[c for c in range(4) if c != i] for i in range(4)]
    u = np.empty(pts.shape[:2] + (4,), dtype=np.int64)   # (B, C, 4)
    for i in range(4):
        m = d[..., cols[i]]                     # (B, C, 3, 3)
        det = (m[..., 0, 0] * (m[..., 1, 1] * m[..., 2, 2] - m[..., 1, 2] * m[..., 2, 1])
               - m[..., 0, 1] * (m[..., 1, 0] * m[..., 2, 2] - m[..., 1, 2] * m[..., 2, 0])
               + m[..., 0, 2] * (m[..., 1, 0] * m[..., 2, 1] - m[..., 1, 1] * m[..., 2, 0]))
        u[..., i] = det if i % 2 == 0 else -det
    nonzero = (u != 0).any(-1)                  # (B, C)
    vals = np.einsum('bcx,bnx->bcn', u, verts)  # (B, C, n)
    v0 = np.take_along_axis(
        vals, subsets[None, :, 0][..., None], axis=2)[..., 0]   # (B, C)
    mx, mn = vals.max(-1), vals.min(-1)
    upper = (mx == v0) & nonzero
    lower = (mn == v0) & nonzero & ~upper
    supp = upper | lower
    sign = np.where(upper, 1, -1)
    g = np.gcd.reduce(np.abs(u), axis=-1)       # (B, C)
    lev_ok = (sign * v0 == g)
    assert bool(lev_ok[supp].all()), "supporting functional not at level 1 (not reflexive?)"
    out = []
    for b in range(B):
        idxs = np.flatnonzero(supp[b])
        seen = {}
        for c in idxs:
            gg = int(g[b, c])
            uu = tuple(int(x) * int(sign[b, c]) // gg for x in u[b, c])
            if uu not in seen:
                onfac = np.flatnonzero(vals[b, c] * int(sign[b, c]) == gg)
                seen[uu] = frozenset(int(x) for x in onfac)
        out.append(list(seen.items()))
    return out


_TRICHO_CACHE = {}

def classify_face_cached(evs, lens, matcher):
    """Face verdict: 'long' (non-unit edge), 'tri' (standard triangle),
    'smoothable', or (status, census_id, edges) for a non-smoothable
    unit-edge face.  Cached by the sorted edge multiset — exact, since every
    verdict we cache is a function of the multiset alone."""
    if any(l >= 2 for l in lens):
        return "long"
    k = len(evs)
    if k == 3:
        a, b = evs[0], evs[1]
        if abs(a[0] * b[1] - a[1] * b[0]) == 1:
            return "tri"                         # standard triangle: smooth
    key = tuple(sorted(evs))
    hit = _TRICHO_CACHE.get(key)
    if hit is not None:
        return hit
    if smoothing_components(list(evs)) > 0:
        _TRICHO_CACHE[key] = "smoothable"
        return "smoothable"
    A2, x, y = 0, 0, 0
    for e in evs:
        A2 += x * e[1] - y * e[0]
        x += e[0]; y += e[1]
    A2 = abs(A2)
    i = (A2 - k + 2) // 2
    st = "rigid" if rigid(list(evs)) else "def-only"
    cid = matcher.match(list(evs), k, A2, i)
    res = (st, cid, [list(e) for e in sorted(evs)])
    _TRICHO_CACHE[key] = res
    return res


def classify_polytope_fast(V, facs, matcher):
    """Same result dict as classify_polytope, but from precomputed facets and
    with the integer polygon path + trichotomy cache."""
    tf = {}
    for (a, b) in combinations(range(len(facs)), 2):
        I = facs[a][1] & facs[b][1]
        if len(I) >= 3 and I not in tf:
            tf[I] = (a, b)
    res = dict(nbad_unit=0, nbad_long=0, statuses=[], census_ids=[],
               noncensus=[], all_unit=True, single_defonly=False,
               n_smooth=0, nfaces=len(tf), max_dual=1, nbad_nontri=0,
               both_sides_hit=False, both_sides_nontri=False)
    bad_faces = []
    for I, (a, b) in tf.items():
        u1, u2 = facs[a][0], facs[b][0]
        dl = vgcd(vsub(u1, u2))                    # l(F*) = edge length of Delta*
        if dl > res["max_dual"]:
            res["max_dual"] = dl
        evs, lens = polygon_of_face_int(V, sorted(I), u1, u2)
        cl = classify_face_cached(evs, lens, matcher)
        if cl == "long":
            res["all_unit"] = False
            res["nbad_long"] += 1
        elif cl == "tri":
            res["n_smooth"] += 1
        elif cl == "smoothable":
            pass
        else:
            st, cid, evlist = cl
            res["nbad_unit"] += 1
            if len(evlist) >= 4:                    # non-triangle => non-quotient
                res["nbad_nontri"] += 1
            res["statuses"].append(st)
            res["census_ids"].append(cid)
            if cid is None:
                res["noncensus"].append(evlist)
            npts = vgcd(vsub(u1, u2))
            bad_faces.append((st, npts))
    if (res["all_unit"] and res["nbad_unit"] == 1 and res["nbad_long"] == 0
            and bad_faces and bad_faces[0][0] == "def-only"
            and bad_faces[0][1] == 1):
        res["single_defonly"] = True
    # both-sides-unit = Delta AND Delta* both unit-edged (max l(F*) == 1);
    # a "hit" also carries a non-smoothable face; non-triangle => non-quotient.
    res["both_sides_hit"] = (res["all_unit"] and res["max_dual"] == 1
                             and res["nbad_unit"] > 0)
    res["both_sides_nontri"] = res["both_sides_hit"] and res["nbad_nontri"] > 0
    return res


# ---------------- census matching (exact GL2(Z) equivalence) ----------------
def _match_gl2(E_face, E_class):
    """Is there M in GL2(Z) with M . E_face = E_class (as multisets)?
    M maps edges to edges, so it is determined by the images of two
    independent edges of the face."""
    tgt = tuple(sorted(E_class))
    a = E_face[0]
    b = next(e for e in E_face[1:] if a[0] * e[1] - a[1] * e[0] != 0)
    det_ab = a[0] * b[1] - a[1] * b[0]
    for ap in E_class:
        for bp in E_class:
            det_ap = ap[0] * bp[1] - ap[1] * bp[0]
            if abs(det_ap) != abs(det_ab):
                continue
            # solve M a = ap, M b = bp exactly
            den = det_ab
            m00 = ap[0] * b[1] - bp[0] * a[1]
            m01 = bp[0] * a[0] - ap[0] * b[0]
            m10 = ap[1] * b[1] - bp[1] * a[1]
            m11 = bp[1] * a[0] - ap[1] * b[0]
            if any(x % den for x in (m00, m01, m10, m11)):
                continue
            M = ((m00 // den, m01 // den), (m10 // den, m11 // den))
            if abs(M[0][0] * M[1][1] - M[0][1] * M[1][0]) != 1:
                continue
            img = tuple(sorted((M[0][0] * e[0] + M[0][1] * e[1],
                                M[1][0] * e[0] + M[1][1] * e[1])
                               for e in E_face))
            if img == tgt:
                return True
    return False


class CensusMatcher:
    def __init__(self):
        self.classes = [(tuple(e), m) for e, m in nonsmoothable_classes()]
        self.by_inv = {}
        for idx, (edges, m) in enumerate(self.classes):
            self.by_inv.setdefault((m["k"], m["A2"], m["i"]), []).append(idx)

    def match(self, evs, k, A2, i):
        for idx in self.by_inv.get((k, A2, i), []):
            if _match_gl2(list(evs), list(self.classes[idx][0])):
                return idx
        return None


# ---------------- per-polytope classification ----------------
def classify_polytope(V, matcher):
    """Returns a compact result dict for one polytope."""
    facs = facets_reflexive(V)
    tf = two_faces_of(V, facs)
    res = dict(nbad_unit=0, nbad_long=0, statuses=[], census_ids=[],
               noncensus=[], all_unit=True, single_defonly=False,
               n_smooth=0, nfaces=len(tf), max_dual=1, nbad_nontri=0,
               both_sides_hit=False, both_sides_nontri=False)
    bad_faces = []
    for I, (a, b) in tf.items():
        u1, u2 = facs[a][0], facs[b][0]
        dl = vgcd(vsub(u1, u2))                    # l(F*) = edge length of Delta*
        if dl > res["max_dual"]:
            res["max_dual"] = dl
        evs, lens = polygon_of_face(V, I, u1, u2)
        k = len(evs)
        A2, x, y = 0, 0, 0
        for e, l in zip(evs, lens):
            A2 += x * (l * e[1]) - y * (l * e[0])
            x += l * e[0]; y += l * e[1]
        A2 = abs(A2)
        bnd = sum(lens)
        i = (A2 - bnd + 2) // 2
        unit = all(l == 1 for l in lens)
        if not unit:
            res["all_unit"] = False
        if k == 3 and A2 == 1:
            res["n_smooth"] += 1
            continue
        sc = smoothing_components(evs) if unit else None
        if unit:
            if sc == 0:
                st = "rigid" if rigid(evs) else "def-only"
                res["nbad_unit"] += 1
                if k >= 4:                          # non-triangle => non-quotient
                    res["nbad_nontri"] += 1
                cid = matcher.match(evs, k, A2, i)
                res["census_ids"].append(cid)
                if cid is None:
                    res["noncensus"].append([list(e) for e in sorted(evs)])
                npts = vgcd(vsub(u1, u2))
                bad_faces.append((st, npts))
                res["statuses"].append(st)
        else:
            # non-isolated germ; candidate only (Filip criterion not implemented)
            # cheap NECESSARY check: no partition of the edge multiset WITH
            # multiplicity into segments/triangles => no 0-mutable candidate
            # of the obvious product form; we only count it.
            res["nbad_long"] += 1
    if (res["all_unit"] and res["nbad_unit"] == 1 and res["nbad_long"] == 0
            and bad_faces and bad_faces[0][0] == "def-only"
            and bad_faces[0][1] == 1):
        res["single_defonly"] = True
    # both-sides-unit = Delta AND Delta* both unit-edged (max l(F*) == 1);
    # a "hit" also carries a non-smoothable face; non-triangle => non-quotient.
    res["both_sides_hit"] = (res["all_unit"] and res["max_dual"] == 1
                             and res["nbad_unit"] > 0)
    res["both_sides_nontri"] = res["both_sides_hit"] and res["nbad_nontri"] > 0
    return res


# ---------------- workers ----------------
_MATCHER = None

def _init_worker():
    global _MATCHER
    _MATCHER = CensusMatcher()


def _accumulate(agg, r, V):
    agg["n"] += 1
    if r["nbad_unit"] > 0:
        agg["nonsmoothable"] += 1
        agg["bad_unit_faces"] += r["nbad_unit"]
        for st in r["statuses"]:
            agg["statuses"][st] += 1
        for cid in r["census_ids"]:
            if cid is not None:
                agg["census_hits"][cid] = agg["census_hits"].get(cid, 0) + 1
        for nc in r["noncensus"]:
            if len(agg["noncensus"]) < 50:
                agg["noncensus"].append(nc)
    if r["nbad_long"] > 0:
        agg["bad_long"] += 1
    if r["single_defonly"]:
        agg["single_defonly"].append([list(map(int, v)) for v in V])
    if r["both_sides_hit"]:                          # Delta & Delta* both unit-edged
        agg["both_sides"] += 1
        if r["both_sides_nontri"]:
            agg["both_sides_nontri"] += 1
        if len(agg["both_sides_polys"]) < 50:
            agg["both_sides_polys"].append(
                dict(V=[list(map(int, v)) for v in V],
                     nontri=bool(r["both_sides_nontri"]),
                     statuses=list(r["statuses"])))


def _empty_agg():
    return dict(n=0, nonsmoothable=0, bad_unit_faces=0, bad_long=0,
                census_hits={}, noncensus=[], single_defonly=[],
                both_sides=0, both_sides_nontri=0, both_sides_polys=[],
                statuses={"rigid": 0, "def-only": 0})


def _work_chunk(rows):
    """Reference engine worker (pure python)."""
    global _MATCHER
    agg = _empty_agg()
    for V in rows:
        r = classify_polytope([tuple(v) for v in V], _MATCHER)
        _accumulate(agg, r, V)
    return agg


def _work_chunk_fast(args):
    """Fast engine worker: numpy facet batch + integer polygon path."""
    global _MATCHER
    import numpy as np
    verts, subsets = args
    agg = _empty_agg()
    B = verts.shape[0]
    # sub-batch to bound memory: B*C*n*8 bytes for vals
    C, n = subsets.shape[0], verts.shape[1]
    step = max(16, min(B, 40_000_000 // max(1, C * n * 8)))
    for s in range(0, B, step):
        vb = verts[s:s + step]
        all_facs = facets_np(vb, subsets)
        for bi in range(vb.shape[0]):
            V = [tuple(int(x) for x in row) for row in vb[bi]]
            r = classify_polytope_fast(V, all_facs[bi], _MATCHER)
            _accumulate(agg, r, V)
    return agg


def merge(a, b):
    a["n"] += b["n"]; a["nonsmoothable"] += b["nonsmoothable"]
    a["bad_unit_faces"] += b["bad_unit_faces"]; a["bad_long"] += b["bad_long"]
    for k, v in b["census_hits"].items():
        a["census_hits"][k] = a["census_hits"].get(k, 0) + v
    for st in a["statuses"]:
        a["statuses"][st] += b["statuses"][st]
    a["noncensus"].extend(b["noncensus"][:max(0, 50 - len(a["noncensus"]))])
    a["single_defonly"].extend(b["single_defonly"])
    a["both_sides"] += b["both_sides"]
    a["both_sides_nontri"] += b["both_sides_nontri"]
    a["both_sides_polys"].extend(
        b["both_sides_polys"][:max(0, 50 - len(a["both_sides_polys"]))])
    return a


# ---------------- validation ----------------
def validate_file(path, nsample=25, seed=20260710, max_dual_points=120):
    """Our Batyrev Hodge numbers must reproduce the dataset's h11/h12 on a
    random sample.  Rows are sampled among those with a small dual polytope:
    hodge_numbers enumerates the dual's facets with the slow exact path, and
    high-vertex-count files contain rows whose duals make that take minutes;
    the convention check is equally decisive on small-dual rows."""
    import pyarrow.parquet as pq
    t = pq.read_table(path, columns=["h11", "h12", "dual_point_count"])
    dpc = t.column("dual_point_count").to_pylist()
    pool = [i for i, d in enumerate(dpc) if d <= max_dual_points]
    if len(pool) < nsample:                      # fall back to smallest duals
        pool = sorted(range(len(dpc)), key=lambda i: dpc[i])[:max(nsample, 50)]
    rng = random.Random(seed)
    idxs = list({pool[rng.randrange(len(pool))] for _ in range(min(nsample, len(pool)))})
    h11s, h12s = t.column("h11").to_pylist(), t.column("h12").to_pylist()
    rows = sample_vertex_rows(path, idxs)
    for i in idxs:
        V = [tuple(v) for v in rows[i]]
        got = hodge_numbers(V, "")
        # the dataset lists polytopes in the dual (Newton-polytope) convention:
        # its h11/h12 are our h21/h11 for the face-fan reading of the vertices.
        # (The KS list is closed under polarity, so sweeping all rows in our
        # convention still covers every reflexive polytope exactly once.)
        assert got == (h12s[i], h11s[i]), \
            f"{path} row {i}: hodge mismatch {got} vs {(h12s[i], h11s[i])}"
    return len(idxs)


# ---------------- driver (streaming I/O: bounded memory at any file size) ----
def file_meta(path):
    """(num_rows, vertex_count) from parquet metadata + one row."""
    import pyarrow.parquet as pq
    pf = pq.ParquetFile(path)
    nrows = pf.metadata.num_rows
    first = next(pf.iter_batches(batch_size=1, columns=["vertex_count"]))
    return nrows, first.column(0)[0].as_py()


def _batch_to_np(batch, n):
    """RecordBatch['vertices'] -> (B, n, 4) int64 with the coordinate-bound
    assert (numpy wraps silently on overflow; this is the guard)."""
    import numpy as np
    flat = batch.column(0).flatten().flatten().to_numpy()
    verts = np.ascontiguousarray(flat.astype(np.int64).reshape(-1, n, 4))
    mx = int(np.abs(verts).max())
    assert mx <= COORD_BOUND, \
        f"coordinate bound violated ({mx} > {COORD_BOUND}): int64 margin not proven"
    return verts


def iter_vert_chunks(path, chunk, limit=None):
    """Stream (B, n, 4) int64 chunks without materializing the file."""
    import pyarrow.parquet as pq
    pf = pq.ParquetFile(path)
    _, n = file_meta(path)
    done = 0
    for batch in pf.iter_batches(batch_size=chunk, columns=["vertices"]):
        verts = _batch_to_np(batch, n)
        if limit is not None:
            if done >= limit:
                return
            verts = verts[:limit - done]
        done += verts.shape[0]
        yield verts


def sample_vertex_rows(path, idxs):
    """Fetch specific rows' vertex lists by streaming batches (never
    materializes the vertices column)."""
    import pyarrow.parquet as pq
    want = set(idxs)
    out = {}
    pf = pq.ParquetFile(path)
    off = 0
    for batch in pf.iter_batches(batch_size=65536, columns=["vertices"]):
        hit = [i for i in want if off <= i < off + batch.num_rows]
        for i in hit:
            out[i] = batch.column(0)[i - off].as_py()
            want.discard(i)
        off += batch.num_rows
        if not want:
            break
    return out


def selftest_engines(path, nsample, seed=20260710):
    """Reference engine == fast engine on random rows (full result dicts,
    order-insensitively).  Streams only the sampled rows."""
    import numpy as np
    rng = random.Random(seed)
    m = CensusMatcher()
    nrows, n = file_meta(path)
    subsets = np.array(list(combinations(range(n), 4)), dtype=np.int64)
    idxs = list({rng.randrange(nrows) for _ in range(min(nsample, nrows))})
    rows = sample_vertex_rows(path, idxs)
    for i in idxs:
        V = [tuple(v) for v in rows[i]]
        verts1 = np.array([rows[i]], dtype=np.int64)
        ref = classify_polytope(V, m)
        fac = facets_np(verts1, subsets)[0]
        fast = classify_polytope_fast(V, fac, m)
        for kk in ("nbad_unit", "nbad_long", "all_unit", "single_defonly",
                   "n_smooth", "nfaces", "max_dual", "both_sides_hit",
                   "both_sides_nontri"):
            assert ref[kk] == fast[kk], (i, kk, ref[kk], fast[kk])
        assert sorted(ref["statuses"]) == sorted(fast["statuses"]), i
        assert sorted(map(str, ref["census_ids"])) == sorted(map(str, fast["census_ids"])), i
        assert sorted(map(str, ref["noncensus"])) == sorted(map(str, fast["noncensus"])), i
    return len(idxs)


def sweep_file(path, procs, chunk=2000, limit=None, engine="fast"):
    import numpy as np
    t0 = time.time()
    agg = _empty_agg()
    if engine == "fast":
        nrows, n = file_meta(path)
        subsets = np.array(list(combinations(range(n), 4)), dtype=np.int64)
        C = subsets.shape[0]
        chunk = max(200, min(20000, 60_000_000 // max(1, C * n),
                             (nrows + 4 * procs - 1) // (4 * procs)))
        # Pool.imap_unordered has NO backpressure (its feeder thread drains the
        # iterable eagerly), so gate the generator with a semaphore released
        # per consumed result: at most ~3 chunks per worker are in flight.
        import threading
        sem = threading.BoundedSemaphore(procs * 3)

        def jobs():
            for verts in iter_vert_chunks(path, chunk, limit):
                sem.acquire()
                yield (verts, subsets)

        with Pool(processes=procs, initializer=_init_worker) as pool:
            for part in pool.imap_unordered(_work_chunk_fast, jobs()):
                sem.release()
                agg = merge(agg, part)
    else:
        import pyarrow.parquet as pq
        t = pq.read_table(path, columns=["vertices"])
        rows = t.column("vertices").to_pylist()
        if limit:
            rows = rows[:limit]
        chunks = [rows[i:i + chunk] for i in range(0, len(rows), chunk)]
        with Pool(processes=procs, initializer=_init_worker) as pool:
            for part in pool.imap_unordered(_work_chunk, chunks):
                agg = merge(agg, part)
    agg["secs"] = round(time.time() - t0, 1)
    return agg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--procs", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument("--limit", type=int, default=None, help="rows per file cap")
    ap.add_argument("--json", default=None)
    ap.add_argument("--engine", choices=["fast", "ref"], default="fast")
    ap.add_argument("--selftest", type=int, default=60,
                    help="rows per file compared between engines (fast mode)")
    ap.add_argument("--no-validate", action="store_true")
    args = ap.parse_args()

    matcher = CensusMatcher()
    print(f"census matcher: {len(matcher.classes)} non-smoothable classes  "
          f"engine={args.engine}", flush=True)

    out = {}
    for path in args.files:
        name = os.path.basename(path)
        if not args.no_validate:
            ns = validate_file(path)
            print(f"{name}: hodge cross-check OK on {ns} random rows", flush=True)
            if args.engine == "fast" and args.selftest:
                nv = selftest_engines(path, args.selftest)
                print(f"{name}: engine selftest (ref == fast) OK on {nv} random rows",
                      flush=True)
        agg = sweep_file(path, args.procs, limit=args.limit, engine=args.engine)
        frac = agg["nonsmoothable"] / agg["n"] if agg["n"] else 0.0
        print(f"{name}: n={agg['n']}  NON-SMOOTHABLE (unit-edge bad face): "
              f"{agg['nonsmoothable']} ({100*frac:.2f}%)  "
              f"[rigid faces {agg['statuses']['rigid']}, def-only faces "
              f"{agg['statuses']['def-only']}]  with non-unit-edge 2-faces: "
              f"{agg['bad_long']}  census classes seen: "
              f"{len(agg['census_hits'])}  outside-box models: "
              f"{len(agg['noncensus'])}  SINGLE-DEFONLY hits: "
              f"{len(agg['single_defonly'])}  [{agg['secs']}s]", flush=True)
        print(f"{name}: both-sides-unit (Delta & Delta* unit-edged) + "
              f"non-smoothable: {agg['both_sides']}  "
              f"** NON-TRIANGLE/non-quotient: {agg['both_sides_nontri']} **",
              flush=True)
        # decode census ids -> class metadata
        agg["census_hits"] = {
            repr(list(matcher.classes[int(k)][0])): v
            for k, v in agg["census_hits"].items()}
        out[name] = agg

    seen = set()
    for a in out.values():
        seen.update(a["census_hits"].keys())
    print(f"\nTOTAL distinct census classes occurring: {len(seen)} / 87", flush=True)
    if args.json:
        with open(args.json, "w") as f:
            json.dump(out, f, indent=1, default=str)
        print(f"results -> {args.json}")


if __name__ == "__main__":
    main()
