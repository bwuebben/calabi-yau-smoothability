#!/usr/bin/env python3
"""Require the example checker to reject altered coordinates, cycles and counts."""
import contextlib, io, json, tempfile
from pathlib import Path
import examples_check as check

HERE=Path(__file__).resolve().parent

def main():
    cases=[
        ('X9','incorrect basis coordinate',lambda d:d['cycle_coordinates_printed'][0].__setitem__(0,2)),
        ('X19','incorrect Euler characteristic',lambda d:d['cycles'][1]['printed'].__setitem__('euler_characteristic',6)),
        ('X19','incorrect intersection count',lambda d:d['cycles'][2]['printed'].__setitem__('interior_gamma_points',40)),
        ('X20','altered boundary coefficient',lambda d:d['cycles'][5]['boundary_legs'][0].__setitem__(2,7)),
        ('X20','altered triangle field',lambda d:d['cycles'][5]['triangles'][0]['field'].__setitem__(0,999)),
    ]
    with tempfile.TemporaryDirectory(prefix='example-controls-') as tmp:
        for key,label,change in cases:
            d=json.loads((HERE/f'{key.lower()}.json').read_text());change(d)
            p=Path(tmp)/'altered.json';p.write_text(json.dumps(d))
            check.FAILURES.clear()
            try:
                with contextlib.redirect_stdout(io.StringIO()):check.run_example(key,p)
            except AssertionError:
                # An invalid coefficient may prevent the later basis solve;
                # require an earlier explicit failed mathematical check.
                if not check.FAILURES:raise
            if not check.FAILURES:raise AssertionError(f'{key}: alteration accepted: {label}')
            print(f'Rejected {key}: {label}')
    print('All five altered examples rejected.')

if __name__=='__main__':main()
