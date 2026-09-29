#!/usr/bin/env python3
"""Tabulate the height-condition census of Remark 3.27(2) from height_condition.json.

For every hexagon-free polytope of the list that satisfies (Proj), the data file records, at the pairs (P, facet) of
a node parallelogram P and a facet containing it, whether the height condition holds for the subdivision Sigma' of
the (Proj) linear program and for its pulling refinement F_C (the file's description gives the definitions). This
program uses only the standard library. It checks the records against each other, against results.json and against
the admissible polytopes of forced_nodes/list_polytopes.json.gz, and recomputes the totals quoted in the paper:
Sigma' violates the height condition on 7 of the 364 polytopes, with largest apex height 2 or 3, and F_C satisfies
it on all of them. It also checks that Sigma' satisfies the condition whenever every boundary lattice point is a ray
of Sigma', as Remark 3.27(1) proves. The subdivisions themselves are not stored, so this is a check of recorded
results, not a recomputation of the cells.
"""
import argparse
import gzip
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXPECTED = dict(polytopes=364, sigma_violating_polytopes=7, sigma_max_apex_heights=[2, 3],
                pulling_satisfying_polytopes=364)


def check(d, census, admissible):
    rows = d['polytopes']
    names = [r['name'] for r in rows]
    assert len(set(names)) == len(names), 'duplicate polytope'
    inside = [x for x in census['hexagon_free'] if not x['status'].startswith('out of scope')]
    assert names == [x['name'] for x in inside], 'polytopes differ from the (Proj) cases of results.json'
    assert all(r['nodes'] == x['nodes'] for r, x in zip(rows, inside)), 'node numbers differ from results.json'
    lst = {p['id']: p for p in admissible['polytopes']}
    scope = {p['id'] for p in lst.values() if not p['hexagons'] and not p.get('equivalent_to')
             and 'witness' in p.get('proj', {})}
    ids = [r['list_id'] for r in rows]
    assert len(set(ids)) == len(ids) and set(ids) == scope, \
        'list ids are not the hexagon-free polytopes of the list with a (Proj) witness'
    for r in rows:
        s, f, name = r['sigma'], r['pulling'], r['name']
        assert lst[r['list_id']]['stored']['nodes'] == r['nodes'], name
        assert r['pairs'] == 2 * r['nodes'], name
        assert 0 <= r['boundary_points_not_rays'] < r['boundary_points'], name
        heights = {int(k): v for k, v in s['vertex_heights'].items()}
        assert min(heights) >= 1 and all(v > 0 for v in heights.values()), name
        assert sum(heights.values()) >= r['pairs'], name
        assert s['max_apex_height'] == max(heights), name
        assert 0 <= s['violating_pairs'] <= r['pairs'], name
        assert (s['violating_pairs'] > 0) == (s['max_apex_height'] >= 2), name
        if r['boundary_points_not_rays'] == 0:
            assert s['violating_pairs'] == 0, f'{name}: Remark 3.27(1) contradicted'
        assert f['violating_pairs'] == 0 and f['max_apex_height'] == 1, name
        assert all(f[k] is True for k in ('facet_volumes', 'every_boundary_point_a_vertex',
                                          'no_other_lattice_points_in_cells', 'two_skeleton_unchanged', 'regular')), name
    bad = [r for r in rows if r['sigma']['violating_pairs']]
    return dict(
        polytopes=len(rows),
        pairs=sum(r['pairs'] for r in rows),
        all_boundary_points_rays=sum(r['boundary_points_not_rays'] == 0 for r in rows),
        sigma_violating_polytopes=len(bad),
        sigma_violating_pairs=sum(r['sigma']['violating_pairs'] for r in rows),
        sigma_max_apex_heights=sorted({r['sigma']['max_apex_height'] for r in bad}),
        sigma_violations=[(r['name'], r['list_id'], r['sigma']['violating_pairs'], r['sigma']['max_apex_height'])
                          for r in bad],
        pulling_satisfying_polytopes=sum(r['pulling']['violating_pairs'] == 0 for r in rows),
        vertex_heights_sigma=dict(sorted(sum((Counter({int(k): v for k, v in r['sigma']['vertex_heights'].items()})
                                              for r in rows), Counter()).items())))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--json', type=Path, help='write recomputed totals')
    a = ap.parse_args()
    d = json.loads((HERE / 'height_condition.json').read_text())
    census = json.loads((HERE / 'results.json').read_text())
    admissible = json.loads(gzip.open(HERE.parent / 'forced_nodes/list_polytopes.json.gz').read())
    v = check(d, census, admissible)
    print(f"{v['polytopes']} hexagon-free polytopes satisfying (Proj); {v['pairs']} pairs (node parallelogram, facet); "
          f"every boundary lattice point is a ray of Sigma' on {v['all_boundary_points_rays']}")
    print(f"Sigma' violates the height condition on {v['sigma_violating_polytopes']} polytopes "
          f"({v['sigma_violating_pairs']} pairs), largest apex heights {v['sigma_max_apex_heights']}:")
    for name, lid, n, h in v['sigma_violations']:
        print(f'  {name} ({lid}): {n} pairs, largest apex height {h}')
    print(f"The pulling refinement F_C satisfies it on {v['pulling_satisfying_polytopes']} of {v['polytopes']}.")
    bad = {k: (v[k], e) for k, e in EXPECTED.items() if v[k] != e}
    if bad:
        raise SystemExit(f'MISMATCH with the numbers printed in the paper: {bad}')
    print('Recorded results are consistent and agree with the paper; the cells are not recomputed here.')
    if a.json:
        a.json.write_text(json.dumps(v, indent=2) + '\n')


if __name__ == '__main__':
    main()
