#!/usr/bin/env python3
"""Both certificate checks on 135 corrupted certificates (Appendix A, "What the proof rests on", item 3).

Run from the repository root with SageMath's Python (the rational check needs Sage):

    sage -python paper6/unimodular/corrupted_certificates_test.py

Fifteen certificates (three chosen at random from each of the files with 7, 12, 18, 24 and 29 vertices) are
modified in nine ways: all heights zero; heights negated; random heights; heights given by a linear function;
a point removed; the point of largest norm removed; three random transpositions of heights; the origin added
as a point; heights divided by 10^5 (rounded down).  This gives 135 modified certificates.  Each is given to
the integer check (check_unimodular_certificates.py) and to the rational check (check_certificates_rational.py).
A modified certificate may still be a valid certificate (for instance when the transposed heights still induce
regular unimodular triangulations); the test asserts that on every one of them
the two programs agree (both accept or both reject), and that both accept the fifteen unmodified certificates.
"""
import copy
import gzip
import importlib.util
import json
import os
import random
import sys

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


INT = load("check_unimodular_certificates")
RAT = load("check_certificates_rational")

rng = random.Random(7)
recs = []
for n in ["07", "12", "18", "24", "29"]:
    with gzip.open(os.path.join(HERE, "certificates", f"heights_v{n}.jsonl.gz"), "rt") as f:
        R = [json.loads(line) for line in f]
    recs += rng.sample(R, 3)


def corrupt(r, how):
    r = copy.deepcopy(r)
    k = len(r["pts"])
    if how == "zero heights":
        r["h"] = [0] * k
    elif how == "negated heights":
        r["h"] = [-x for x in r["h"]]
    elif how == "random heights":
        r["h"] = [rng.randrange(10 ** 9) for _ in range(k)]
    elif how == "linear heights":
        a = [rng.randrange(-50, 50) for _ in range(4)]
        r["h"] = [sum(x * y for x, y in zip(a, p)) for p in r["pts"]]
    elif how == "a point removed":
        i = rng.randrange(k)
        del r["pts"][i]
        del r["h"][i]
    elif how == "point of largest norm removed":
        i = max(range(k), key=lambda j: sum(x * x for x in r["pts"][j]))
        del r["pts"][i]
        del r["h"][i]
    elif how == "heights transposed":
        for _ in range(3):
            i, j = rng.sample(range(k), 2)
            r["h"][i], r["h"][j] = r["h"][j], r["h"][i]
    elif how == "origin added":
        r["pts"].append([0, 0, 0, 0])
        r["h"].append(0)
    elif how == "heights rescaled":
        r["h"] = [x // 10 ** 5 for x in r["h"]]
    return r


def accepts_int(r):
    return INT.safe(r)["ok"]


def accepts_rat(r):
    return RAT.safe(("test", r))["ok"]


HOWS = ["zero heights", "negated heights", "random heights", "linear heights", "a point removed",
        "point of largest norm removed", "heights transposed", "origin added", "heights rescaled"]

if __name__ == "__main__":
    a0 = sum(accepts_int(r) for r in recs)
    b0 = sum(accepts_rat(r) for r in recs)
    print(f"unmodified: accepted by the integer check {a0}/{len(recs)}, by the rational check {b0}/{len(recs)}")
    assert a0 == b0 == len(recs)
    total = disagree = 0
    for how in HOWS:
        a = b = 0
        for r in recs:
            c = corrupt(r, how)
            va, vb = accepts_int(c), accepts_rat(c)
            a += va
            b += vb
            total += 1
            if va != vb:
                disagree += 1
                print("  DIFFERENT ANSWERS", how, r["row"], va, vb)
        print(f"{how:30s}: accepted by the integer check {a}/{len(recs)}, by the rational check {b}/{len(recs)}",
              flush=True)
    print(f"{total} modified certificates, different answers on {disagree}")
    assert total == 135 and disagree == 0
