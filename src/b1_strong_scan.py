#!/usr/bin/env python3
"""
Test conjecture B1-strong at scale: is EVERY 2-face of a both-sides-unit reflexive
4-polytope either a triangle or a CENTRALLY-SYMMETRIC polygon (zonotope)?

For every both-sides-unit polytope (Delta unit-edge AND max l(F*)==1), inspect its
non-triangle 2-faces; flag any that is NOT centrally symmetric (edge multiset !=
its own negation).  0 flags across the scan => B1-strong holds on the population,
which implies B1 (mirror non-isolation).

Run (multiprocessed): ./venv/bin/python src/b1_strong_scan.py data/ks/polytopes-4d-0{5,6,7,8,9}-vertices.parquet
"""
from collections import Counter
from multiprocessing import Pool
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from toric_census import rigid, smoothing_components
from ks_sweep import facets_reflexive, two_faces_of, polygon_of_face
from batyrev_global import vgcd, vsub


def centrally_symmetric(evs):
    c = Counter(tuple(e) for e in evs)
    neg = Counter((-e[0], -e[1]) for e in evs)
    return c == neg


def classify(V):
    facs = facets_reflexive(V)
    tf = two_faces_of(V, facs)
    all_unit = True
    max_dual = 1
    faces = []                                    # (k, evs, status_is_nonsmooth)
    for I, (a, b) in tf.items():
        u1, u2 = facs[a][0], facs[b][0]
        evs, lens = polygon_of_face(V, I, u1, u2)
        dl = vgcd(vsub(u1, u2))
        max_dual = max(max_dual, dl)
        if any(l >= 2 for l in lens):
            all_unit = False
        faces.append((len(evs), tuple(evs)))
    return all_unit, max_dual, faces


def _work(rows):
    n = both = 0
    viol = []                       # (k, evs) non-triangle, non-centrally-symmetric on a both-sides poly
    viol_nonsmooth = 0
    for V in rows:
        V = [tuple(int(x) for x in v) for v in V]
        all_unit, max_dual, faces = classify(V)
        n += 1
        if not (all_unit and max_dual == 1):
            continue
        both += 1
        for k, evs in faces:
            if k >= 4 and not centrally_symmetric(evs):
                if len(viol) < 20:
                    viol.append((k, list(evs), V))
                if smoothing_components(list(evs)) == 0:
                    viol_nonsmooth += 1
    return n, both, viol, viol_nonsmooth


def main():
    import pyarrow.parquet as pq
    procs = max(1, (os.cpu_count() or 4) - 2)
    N = B = VN = 0
    VIOL = []
    for path in sys.argv[1:]:
        t0 = time.time()
        pf = pq.ParquetFile(path)
        with Pool(procs) as pool:
            pend = []
            for batch in pf.iter_batches(batch_size=4000, columns=["vertices"]):
                pend.append(pool.apply_async(_work, (batch.column(0).to_pylist(),)))
                if len(pend) >= procs * 3:
                    n, both, viol, vn = pend.pop(0).get()
                    N += n; B += both; VN += vn
                    VIOL.extend(viol[:max(0, 20 - len(VIOL))])
            for p in pend:
                n, both, viol, vn = p.get()
                N += n; B += both; VN += vn
                VIOL.extend(viol[:max(0, 20 - len(VIOL))])
        print(f"{os.path.basename(path)}: cumulative n={N} both-sides={B} "
              f"non-tri-non-symmetric-faces={len(VIOL)} (of which non-smoothable={VN})"
              f"  [{time.time()-t0:.0f}s]", flush=True)
    print("\n" + "=" * 60)
    print(f"TOTAL: n={N}  both-sides-unit={B}")
    print(f"non-triangle NON-centrally-symmetric both-sides faces: {len(VIOL)} "
          f"(non-smoothable among them: {VN})")
    if not VIOL:
        print("==> B1-STRONG HOLDS on this population: every both-sides 2-face is a\n"
              "    triangle or a centrally-symmetric zonotope.  (Implies B1.)")
    else:
        print("==> B1-STRONG FALSE: found non-triangle non-symmetric both-sides faces:")
        for k, evs, V in VIOL[:10]:
            print(f"    k={k} smoothable={'yes' if smoothing_components(list(evs)) else 'NO'} "
                  f"edges={evs}  V={V}")


if __name__ == "__main__":
    main()
