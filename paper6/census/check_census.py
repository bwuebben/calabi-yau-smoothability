#!/usr/bin/env python3
"""Tabulate the retained tropical-cycle computation results of Section 11.

This checks record consistency and recomputes the reported totals. It does not
reconstruct hundreds of bases or certify the underlying cycle computations.
The file distinguishes per-case records from the two aggregate/report-only
items. The README explains those limits.
"""
import argparse
from collections import Counter
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def counts(d):
    free = d['hexagon_free']
    assert len({x['name'] for x in free}) == len(free)
    inside = [x for x in free if not x['status'].startswith('out of scope')]
    realised = [x for x in inside if x['realised_equals_K']]
    partial = [x for x in inside if not x['realised_equals_K']]
    for x in realised:
        assert x['realised_rank'] == x['rank_K'] and x['realised_index'] == 1
    for x in partial:
        assert 0 <= x['realised_rank'] < x['rank_K']
    assert all(x['boundary_shaped_span_K'] and x['boundary_index'] == 1 for x in inside)
    coranks = [x['rank_K']-x['realised_rank'] for x in partial]
    negative = Counter(x['min_negative_vertices_needed'] for x in inside)
    belt = d['belt_summary']
    assert belt['hexagon_free']['belts'] == belt['hexagon_free']['disc'] == belt['hexagon_free']['one_field']
    assert len(belt['hexagon_free_nonconforming']) == belt['hexagon_free_gamma']['belts_meeting_Gamma_inside']
    relative = d['relative_bases']
    assert len({x['name'] for x in relative}) == len(relative)
    assert all(x['balanced_coefficients_span_K'] and x['all_balanced_networks_bound'] and x['index_one'] for x in relative)
    hx = d['relative_hexagonal_bases']
    assert len({(x['name'], x['profile']) for x in hx}) == len(hx)
    assert all(x['balanced_coefficients_span_K'] and x['all_balanced_networks_bound'] for x in hx)
    hexes = d['hexagonal']
    assert len({(x['name'],x['group']) for x in hexes}) == len(hexes)
    for x in hexes:
        if x['realised_equals_K']:
            assert x['realised_rank'] == x['rank_K'] and x['index'] == 1
    a1 = [x for x in hexes if x['group'] == 'A1']
    a2 = [x for x in hexes if x['group'] == 'A2_or_mixed']
    enum = d['enumerated_boundary_summary']['polytopes']
    b = belt['hexagon_free_belt_spans']
    return dict(NList=len(free), NPfail=len(free)-len(inside), NScope=len(inside),
        NBdryEnum=enum, NBdryIP=len(inside)-enum, NMinNegZero=negative[0],
        NMinNegPos=len(inside)-negative[0], NMinNegTwo=negative[2], NMinNegFour=negative[4],
        NMinNegSix=negative[6], NMinNegEight=negative[8], NBelts=b['belts_generate_K'],
        NRealised=len(realised), NRemaining=len(partial), NCorankOne=coranks.count(1),
        NCorankMax=max(coranks), NFailSheet=sum('sheet' in x['failure_types'] for x in partial),
        NFailFill=sum('filling' in x['failure_types'] for x in partial),
        NFailLocal=sum('local_replacement' in x['failure_types'] for x in partial),
        NThmCBases=len(relative), NThmCGlobal=sum(not x['facet_networks_span_K'] for x in relative),
        NThmCHexBases=len(hx), NHexRealAone=sum(x['realised_equals_K'] for x in a1),
        NHexRealAtwo=sum(x['realised_equals_K'] for x in a2), NHexProfAtwo=len(a2),
        NBeltsStrict=b['conforming_belts_generate_K'],
        NBeltsCrossPoly=b['belts_generate_K']-b['conforming_belts_generate_K'],
        NBeltsChecked=belt['hexagon_free']['belts'],
        NBeltsCross=len(belt['hexagon_free_nonconforming']))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--json', type=Path, help='write recomputed totals')
    a = ap.parse_args()
    values = counts(json.loads((HERE / 'results.json').read_text()))
    for key, value in values.items():
        print(f'{key}: {value}')
    print('Recorded summaries are internally consistent; no geometric replay is claimed.')
    if a.json:
        a.json.write_text(json.dumps(values, indent=2) + '\n')


if __name__ == '__main__':
    main()
