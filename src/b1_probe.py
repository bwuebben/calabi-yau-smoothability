#!/usr/bin/env python3
"""
B1 probe (mirror-non-isolation, Track B1).  For a reflexive 4-polytope, does a
unit-edge 2-face F with k>=4 vertices and i(F)>=1 interior points FORCE some
2-face G with l(G*)>=2 (i.e. Delta* non-unit-edge, mirror non-isolated)?

We test the exact trigger on the 6-vertex KS file (smallest with non-triangle
faces): tabulate, per polytope, whether it has a unit non-triangle face of each
kind {square(k=4,i=0), (k>=4,i>=1)-smoothable, non-smoothable-nontriangle} and
its max l(F*).  A both-sides candidate = all Delta edges unit AND max l(F*)==1.

Run:  ./venv/bin/python src/b1_probe.py data/ks/polytopes-4d-06-vertices.parquet
"""
from collections import Counter
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from toric_census import rigid, smoothing_components
from ks_sweep import facets_reflexive, two_faces_of, polygon_of_face
from batyrev_global import vgcd, vsub


def face_stats(V):
    facs = facets_reflexive(V)
    tf = two_faces_of(V, facs)
    all_unit = True
    max_dual = 1
    kinds = set()                      # which non-triangle unit face kinds occur
    # for the implication: does a (k>=4,i>=1) unit face coexist with max l == 1?
    has_k4_i1 = False
    has_square = False
    has_nontri_nonsmooth = False
    for I, (a, b) in tf.items():
        u1, u2 = facs[a][0], facs[b][0]
        evs, lens = polygon_of_face(V, I, u1, u2)
        dl = vgcd(vsub(u1, u2))
        if dl > max_dual:
            max_dual = dl
        if any(l >= 2 for l in lens):
            all_unit = False
            continue
        k = len(evs)
        if k == 3:
            continue
        # unit non-triangle face: compute i via Pick (2A = 2i + k - 2)
        A2 = 0; x = 0; y = 0
        for e in evs:
            A2 += x * e[1] - y * e[0]; x += e[0]; y += e[1]
        A2 = abs(A2)
        i = (A2 - k + 2) // 2
        sc = smoothing_components(list(evs))
        if k == 4 and i == 0:
            has_square = True
        if i >= 1:
            has_k4_i1 = True
            if sc == 0:
                has_nontri_nonsmooth = True
    return dict(all_unit=all_unit, max_dual=max_dual,
                has_k4_i1=has_k4_i1, has_square=has_square,
                has_nontri_nonsmooth=has_nontri_nonsmooth)


def main():
    import pyarrow.parquet as pq
    path = sys.argv[1]
    t0 = time.time()
    # counters for the implications
    n = 0
    # A) (k>=4,i>=1) face AND max l == 1  (would REFUTE the i>=1 trigger)
    refute_i1 = 0
    # B) square face, Delta unit-edge, max l == 1  (square on a both-sides poly)
    square_both = 0
    # C) non-triangle non-smoothable face AND max l == 1 (REFUTE the headline B1)
    refute_headline = 0
    examples = []
    pf = pq.ParquetFile(path)
    for batch in pf.iter_batches(batch_size=5000, columns=["vertices"]):
        for V in batch.column(0).to_pylist():
            V = [tuple(int(x) for x in v) for v in V]
            r = face_stats(V)
            n += 1
            if r["has_k4_i1"] and r["all_unit"] and r["max_dual"] == 1:
                refute_i1 += 1
                if len(examples) < 5:
                    examples.append(("i1_maxl1", V))
            if r["has_square"] and r["all_unit"] and r["max_dual"] == 1:
                square_both += 1
            if r["has_nontri_nonsmooth"] and r["all_unit"] and r["max_dual"] == 1:
                refute_headline += 1
    print(f"n={n}  [{time.time()-t0:.0f}s]")
    print(f"(k>=4,i>=1) unit face WITH Delta unit-edge & max l(F*)=1  : {refute_i1}"
          f"   <- 0 => 'i>=1 & k>=4 => mirror non-isolated' holds on v06")
    print(f"square (k=4,i=0) face WITH Delta unit-edge & max l(F*)=1   : {square_both}"
          f"   <- >0 => squares CAN be both-sides (trigger is i>=1, not k>=4)")
    print(f"non-triangle NON-SMOOTHABLE face & Delta unit & max l=1    : {refute_headline}"
          f"   <- 0 => headline B1 holds on v06")
    for tag, V in examples:
        print(f"    {tag}: {V}")


if __name__ == "__main__":
    main()
