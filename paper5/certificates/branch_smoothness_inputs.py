#!/usr/bin/env python3
"""Exact inputs to analytic globalization and branch smoothness.

Run: sage -python paper5/certificates/branch_smoothness_inputs.py
This verifies the counterexample's Hodge numbers, active branch projections,
and the isomorphisms between the original and nodal relation kernels.
Smoothness of deformation bases is proved in the accompanying paper,
not certified by these finite calculations.
"""
from itertools import product
import json
from pathlib import Path
from sage.all import Polyhedron, QQ, matrix, vector
from a2_branch_selection_counterexample import V


def serial(v):
    return [str(x) for x in v]


def main():
    checks = []
    def check(name,condition):
        assert condition,name
        checks.append(name)
    data = json.loads(Path(__file__).with_name('a2_branch_selection_counterexample.json').read_text())
    check('same explicit polytope',data['vertices'] == [list(v) for v in V])
    P = Polyhedron(vertices=V,base_ring=QQ)
    polar = Polyhedron(vertices=[list(-vector(h[1:])/h[0]) for h in P.inequality_generator()],base_ring=QQ)
    pts = list(polar.integral_points())
    facet_interiors = []
    for face in polar.faces(3):
        F = face.as_polyhedron()
        facet_interiors += [v for v in F.integral_points() if F.relative_interior_contains(v)]
    check('polar lattice count',len(pts) == 26)
    check('no polar facet interior points',not facet_interiors)
    # The mirror correction vanishes because every original edge is unit.
    check('all original edges have two lattice points',all(
          len(face.as_polyhedron().integral_points()) == 2 for face in P.faces(1)))
    h21 = len(pts)-5-len(facet_interiors)
    check('resolution h21',h21 == 21)
    h,e1,e2,e3 = matrix.identity(QQ,4).rows()
    J = matrix.diagonal(QQ,[1,-1,-1,-1])
    C = [e1,h-e1-e2,e2,h-e2-e3,e3,h-e1-e3]
    canonical = -3*h+e1+e2+e3
    for i,curve in enumerate(C):
        check(f'boundary curve {i} square -1',curve*J*curve == -1)
        check(f'boundary curve {i} anticanonical degree one',-canonical*J*curve == 1)
        check(f'boundary plus anticanonical nef {i}',all(
              (curve-canonical)*J*other >= 0 for other in C))
    reports = []
    for rec in data['profiles']:
        A = matrix(QQ,rec['matrix'])
        K = A.left_kernel().basis_matrix()
        active_ranks = []
        for i in range(2):
            ids = [j for j,label in enumerate(rec['labels']) if label.startswith(f'dp6_{i}:')]
            active_ranks.append(int(K.matrix_from_columns(ids).rank()))
        check(rec['profile']+' both selected branches active',all(active_ranks))
        rows = [vector(QQ,list(r)+[0,0]) for r in A.rows()[:30]]
        labels = [f'old_node_{i}' for i in range(30)]
        entries = [(j,j,1) for j in range(30)]
        original = 30
        for i,ch in enumerate(rec['profile']):
            start = len(rows)
            contracted = [e1,h-e2-e3] if ch == 'L' else [e1,e2,e3]
            for k,curve in enumerate(contracted):
                r = vector(QQ,27)
                for vertex,boundary in zip(data['dp6_faces'][i],C):
                    r[vertex] = curve*J*boundary
                r[25+i] = curve*J*canonical
                rows.append(r);labels.append(f'new_node_{i}:{k}')
            if ch == 'L':
                # ell=e1-(h-e2-e3).
                entries += [(original,start,1),(original,start+1,-1)]
                original += 1
            else:
                # a*A+b*B=-a*e1+b*e2+(a-b)*e3.
                entries += [(original,start,-1),(original,start+2,1),
                            (original+1,start+1,1),(original+1,start+2,-1)]
                original += 2
        N = matrix(QQ,rows)
        U = matrix(QQ,A.nrows(),N.nrows())
        for i,j,x in entries:
            U[i,j] = x
        KN = N.left_kernel().basis_matrix()
        check(rec['profile']+' exact curve substitution',U*N == A.augment(matrix(QQ,A.nrows(),2)))
        check(rec['profile']+' substitution injective',U.rank() == A.nrows())
        check(rec['profile']+' kernel maps to kernel',K*U*N == 0)
        check(rec['profile']+' equal kernel dimensions',K.nrows() == KN.nrows())
        check(rec['profile']+' kernel isomorphism',(K*U).row_space() == KN.row_space())
        forced = [j for j in range(N.nrows()) if KN.column(j) == 0]
        check(rec['profile']+' exactly one forced nodal coefficient',len(forced) == 1)
        reports.append(dict(profile=rec['profile'],active_projection_ranks=active_ranks,
            nodal_count=N.nrows(),nodal_relation_rank=int(N.rank()),
            kernel_dimension=K.nrows(),locally_trivial_dimension=h21,
            predicted_branch_dimension=h21+K.nrows(),
            forced_node=labels[forced[0]],
            node_labels=labels,nodal_matrix=[serial(r) for r in N],
            substitution=[serial(r) for r in U]))
    check('predicted dimensions',[r['predicted_branch_dimension'] for r in reports] == [34,34,34,35])
    target = Path(__file__).with_suffix('.json')
    target.write_text(json.dumps(dict(polar_lattice_points=[serial(v) for v in pts],h21=h21,
        reports=reports,checks=checks,passed=len(checks),
        scope='finite inputs; analytic smoothness and blow-down are proved separately'),indent=2)+'\n')
    print(f'PASS: {len(checks)} exact branch-smoothness input checks')
    for r in reports:
        print(r['profile'],'active ranks',r['active_projection_ranks'],
              'nodes',r['nodal_count'],'kernel',r['kernel_dimension'],
              'predicted base dimension',r['predicted_branch_dimension'])


if __name__ == '__main__':
    main()
