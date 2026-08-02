#!/usr/bin/env python3
"""
Find every BOTH-SIDES-UNIT reflexive polytope in the given KS files and dump the
full 2-face inventory of each: for each 2-face, (k, i, l(F*), status).  Goal:
see an explicit SMOOTHABLE non-triangle both-sides face (self-duality says one
must exist), to understand what non-smoothability rules out.

Run: ./venv/bin/python src/b1_dump.py data/ks/polytopes-4d-0{5,6,7}-vertices.parquet
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from toric_census import rigid, smoothing_components
from ks_sweep import facets_reflexive, two_faces_of, polygon_of_face
from batyrev_global import vgcd, vsub
import pyarrow.parquet as pq


def faces_full(V):
    facs = facets_reflexive(V)
    tf = two_faces_of(V, facs)
    all_unit = True
    max_dual = 1
    rows = []
    for I, (a, b) in tf.items():
        u1, u2 = facs[a][0], facs[b][0]
        evs, lens = polygon_of_face(V, I, u1, u2)
        dl = vgcd(vsub(u1, u2))
        max_dual = max(max_dual, dl)
        if any(l >= 2 for l in lens):
            all_unit = False
        k = len(evs)
        A2 = 0; x = 0; y = 0
        for e in evs:
            A2 += x * e[1] - y * e[0]; x += e[0]; y += e[1]
        A2 = abs(A2); i = (A2 - k + 2) // 2
        if k == 3 and A2 == 1:
            st = "smooth"
        elif rigid(list(evs)):
            st = "RIGID"
        elif smoothing_components(list(evs)) == 0:
            st = "DEF-ONLY"
        else:
            st = "smoothable"
        rows.append((k, i, dl, st, tuple(sorted(evs))))
    return all_unit, max_dual, rows


def main():
    total = 0
    for path in sys.argv[1:]:
        pf = pq.ParquetFile(path)
        nv = os.path.basename(path)
        for batch in pf.iter_batches(batch_size=5000, columns=["vertices"]):
            for V in batch.column(0).to_pylist():
                V = [tuple(int(x) for x in v) for v in V]
                all_unit, max_dual, rows = faces_full(V)
                if all_unit and max_dual == 1:           # BOTH-SIDES-UNIT
                    total += 1
                    nontri = [r for r in rows if r[0] >= 4]
                    print(f"\n=== BOTH-SIDES #{total} [{nv}]  V={V}")
                    print(f"    #2-faces={len(rows)}  non-triangle 2-faces={len(nontri)}")
                    from collections import Counter
                    kinds = Counter((k, st) for k, i, dl, st, _ in rows)
                    print(f"    face (k,status) histogram: {dict(kinds)}")
                    for (k, i, dl, st, evs) in rows:
                        if k >= 4 or st != "smooth":
                            print(f"      k={k} i={i} l(F*)={dl} {st:11s} edges={list(evs)}")
    print(f"\nTOTAL both-sides-unit found: {total}")


if __name__ == "__main__":
    main()
