#!/usr/bin/env python3
"""Independent exact certificate for a counterexample to the weak criterion.

Run: sage -python paper5/certificates/a2_branch_selection_counterexample.py
The proof of nonsmoothability additionally uses the simultaneous-resolution
necessity argument. This script certifies its polyhedral and linear inputs;
it does not certify an analytic theorem by numerical experimentation.
No import of the search, its data, or its face/circuit routines is used.
"""
from itertools import combinations, product
import json
from pathlib import Path
from random import Random
from sage.all import Polyhedron, QQ, ZZ, gcd, lcm, matrix, vector

V = [(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),
     (1,-1,1,-1),(0,0,1,-1),(1,-1,0,0),(-1,1,0,0),
     (0,0,-1,1),(-1,1,-1,1),(0,0,0,-1),(0,0,-1,0),
     (-1,1,-1,0),(0,-1,1,0),(-1,0,1,0),(0,-1,0,1),
     (-1,0,0,1),(-1,0,1,-1),(0,-1,1,-1),(-1,0,-1,1),
     (0,-1,-1,1),(0,-1,-1,0),(0,-1,0,-1),(-1,0,0,-1),
     (-1,0,-1,0)]


def index3(rows):
    B = matrix(ZZ, rows)
    return gcd([abs(B.matrix_from_columns(cols).det())
                for cols in combinations(range(4), 3)])


def serial(v):
    return [str(x) for x in v]


def main():
    checks = []
    def check(name, condition):
        assert condition, name
        checks.append(name)
    P = Polyhedron(vertices=V, base_ring=QQ)
    check('dimension and vertices', P.dim() == 4 and P.n_vertices() == 25)
    check('reflexivity', all(h[0] > 0 and all(x/h[0] in ZZ for x in h[1:])
                             for h in P.inequality_generator()))
    normals = [vector(QQ, h[1:])/h[0] for h in P.inequality_generator()]
    pos = {tuple(v): i for i, v in enumerate(V)}
    edges = {frozenset(pos[tuple(v)] for v in e.vertices()) for e in P.faces(1)}
    check('all 71 edges primitive', len(edges) == 71 and all(
        gcd(list(vector(V[i])-vector(V[j]))) == 1 for i,j in map(tuple, edges)))
    nodes, sixes, centers = [], [], []
    smooth = 0
    for face in P.faces(2):
        ids = sorted(pos[tuple(v)] for v in face.vertices())
        if len(ids) == 3:
            check(f'unimodular triangle {ids}', index3([V[i] for i in ids]) == 1)
            smooth += 1
            continue
        check(f'allowed polygon size {ids}', len(ids) in (4,6))
        supporting = [n for n in normals if all(n*vector(V[i]) == -1 for i in ids)]
        check(f'dual edge length one {ids}', len(supporting) == 2 and
              gcd([ZZ(x) for x in supporting[0]-supporting[1]]) == 1)
        adjacent = {i: sorted(j for j in ids if frozenset((i,j)) in edges) for i in ids}
        check(f'boundary graph {ids}', all(len(a) == 2 for a in adjacent.values()))
        cyc = [ids[0], adjacent[ids[0]][0]]
        while len(cyc) < len(ids):
            cyc.append(next(j for j in adjacent[cyc[-1]] if j != cyc[-2]))
        check(f'boundary cycle {ids}', cyc[0] in adjacent[cyc[-1]] and
              len(set(cyc)) == len(ids))
        if len(ids) == 4:
            check(f'unit node polygon {ids}', vector(V[cyc[0]])+vector(V[cyc[2]]) ==
                  vector(V[cyc[1]])+vector(V[cyc[3]]) and
                  index3([V[i] for i in cyc[:3]]) == 1)
            nodes.append(cyc)
        else:
            center = sum((vector(QQ,V[i]) for i in ids), vector(QQ,4))/6
            check(f'integral center {ids}', all(x in ZZ for x in center))
            check(f'unimodular hexagon star {ids}', all(
                index3([center,V[cyc[i]],V[cyc[(i+1)%6]]]) == 1 for i in range(6)))
            check(f'all six boundary curves have square -1 {ids}', all(
                vector(V[cyc[i-1]])+vector(V[cyc[(i+1)%6]]) ==
                vector(V[cyc[i]])+center for i in range(6)))
            sixes.append(cyc)
            centers.append(tuple(center))
    nodes.sort(key=lambda c: sorted(c))
    sixes.sort(key=lambda c: sorted(c))
    check('complete singularity inventory', len(nodes) == 30 and len(sixes) == 2 and smooth == 37)
    lattice_points = set(map(tuple, P.integral_points()))
    check('all lattice points listed', lattice_points ==
          set(V) | set(centers) | {(0,0,0,0)} and len(lattice_points) == 28)
    # There are no facet-interior points. Each two-face with an interior
    # point has dual edge length one: Batyrev's non-toric correction is zero.
    # Thus h11=28-5=23, and the 27 boundary divisors span H2 dually.
    # The two added center columns in every root/circuit row are zero.
    h,e1,e2,e3 = matrix.identity(QQ,4).rows()
    J = matrix.diagonal(QQ,[1,-1,-1,-1])
    curves = [e1,h-e1-e2,e2,h-e2-e3,e3,h-e1-e3]
    local_classes = [-(h-e1-e2-e3),-e1+e3,e2-e3]
    canonical = -3*h+e1+e2+e3
    check('root classes pair to zero with exceptional divisor',
          all(alpha*J*canonical == 0 for alpha in local_classes))
    check('A1 and A2 root systems',
          matrix(QQ,[[a*J*b for b in local_classes] for a in local_classes]) ==
          matrix(QQ,[[-2,0,0],[0,-2,1],[0,1,-2]]))
    def embed(cyc, coeff):
        r = vector(QQ,25)
        for i,x in zip(cyc,coeff):
            r[i] = x
        check('curve functional annihilates principal divisors',
              all(sum(r[i]*V[i][j] for i in range(25)) == 0 for j in range(4)))
        return r
    fixed = [embed(c,[1,-1,1,-1]) for c in nodes]
    blocks = [[embed(c,[alpha*J*D for D in curves]) for alpha in local_classes]
              for c in sixes]
    reports = []
    for choice in product('LP',repeat=2):
        rows = list(fixed)
        labels = [f'node_{i}' for i in range(30)]
        pairs = []
        for i,ch in enumerate(choice):
            j = len(rows)
            if ch == 'L':
                rows.append(blocks[i][0]); labels.append(f'dp6_{i}:L')
                pairs.append((j,None))
            else:
                rows += blocks[i][1:]; labels += [f'dp6_{i}:a',f'dp6_{i}:b']
                pairs.append((j,j+1))
        A = matrix(QQ,rows)
        K = A.left_kernel().basis_matrix()
        forms = matrix.identity(QQ,A.nrows()).rows()
        flabels = list(labels)
        for i,(j,k) in enumerate(pairs):
            if k is not None:
                forms.append(forms[j]-forms[k]); flabels.append(f'dp6_{i}:a-b')
        forced = [(label,f) for label,f in zip(flabels,forms) if K*f == 0]
        weak = all(K.column(j) != 0 for j in range(30)) and all(
            K.column(j) != 0 or (k is not None and K.column(k) != 0) for j,k in pairs)
        check(f'{choice} exactly one obstruction', len(forced) == 1)
        name,f = forced[0]
        D = A.solve_right(f)
        scale = lcm(x.denominator() for x in D)
        D *= scale
        check(f'{choice} integral dual obstruction', A*D == scale*f and all(x in ZZ for x in D))
        rec = dict(profile=''.join(choice), labels=labels,
                   matrix=[serial(r) for r in A], rank=int(A.rank()),
                   kernel_dimension=int(K.nrows()), weak_pass=bool(weak),
                   forced=name, forced_functional=serial(f),
                   dual_divisor=serial(D), dual_scale=str(scale))
        # Independently re-express necessity on the simultaneous partial
        # model. A2 is replaced by three actual (-1,-1) curves e1,e2,e3;
        # their intersection with the appropriate center divisor is -1.
        # A forced node here proves nonsmoothability directly, without
        # invoking the new sufficiency theorem or a tangent approximation.
        nodal_rows = [vector(QQ,list(r)+[0,0]) for r in fixed]
        nodal_labels = [f'old_node_{i}' for i in range(30)]
        node_indices = list(range(30))
        for i,ch in enumerate(choice):
            if ch == 'L':
                nodal_rows.append(vector(QQ,list(blocks[i][0])+[0,0]))
                nodal_labels.append(f'dp6_{i}:L')
            else:
                for k,exceptional_curve in enumerate((e1,e2,e3)):
                    r = vector(QQ,27)
                    for vertex,Dcurve in zip(sixes[i],curves):
                        r[vertex] = exceptional_curve*J*Dcurve
                    r[25+i] = exceptional_curve*J*canonical
                    node_indices.append(len(nodal_rows))
                    nodal_rows.append(r)
                    nodal_labels.append(f'new_node_{i}:e{k+1}')
        N = matrix(QQ,nodal_rows)
        KN = N.left_kernel().basis_matrix()
        forced_nodes = [j for j in node_indices if KN.column(j) == 0]
        check(f'{choice} a node is forced on simultaneous partial model',bool(forced_nodes))
        j = forced_nodes[0]
        nf = vector(QQ,[int(k == j) for k in range(N.nrows())])
        ND = N.solve_right(nf)
        nscale = lcm(x.denominator() for x in ND)
        ND *= nscale
        check(f'{choice} direct nodal dual obstruction',N*ND == nscale*nf)
        rec['partial_model'] = dict(labels=nodal_labels,matrix=[serial(r) for r in N],
            forced_nodes=[nodal_labels[j] for j in forced_nodes],
            certified_node=nodal_labels[j],functional=serial(nf),
            dual_divisor=serial(ND),dual_scale=str(nscale),
            kernel_dimension=int(KN.nrows()))
        if weak:
            required = [f for f in forms if K*f != 0]
            best = None
            rng = Random(20260905)
            for _ in range(4000):
                w = vector(QQ,[rng.randint(-6,6) for _ in range(K.nrows())])*K
                if all(w*f for f in required):
                    score = (max(abs(x) for x in w),sum(abs(x) for x in w))
                    if best is None or score < best[0]:
                        best = (score,w)
            if best is None:
                # Deterministic finite polynomial-avoidance bound.
                for t in range(1,len(required)*max(1,K.nrows()-1)+2):
                    w = vector(QQ,[t**i for i in range(K.nrows())])*K
                    if all(w*f for f in required):
                        best = (None,w); break
            check(f'{choice} explicit weak relation found', best is not None)
            w = best[1]
            check(f'{choice} exact relation and nonvanishing', w*A == 0 and
                  all(w*f for f in required) and all(w[j] or (k is not None and w[k])
                                                       for j,k in pairs))
            rec['weak_witness'] = serial(w)
        reports.append(rec)
    check('four profiles exhausted', [r['profile'] for r in reports] == ['LL','LP','PL','PP'])
    check('weak premise holds exactly on mixed profiles',
          [r['profile'] for r in reports if r['weak_pass']] == ['LP','PL'])
    check('uniform profiles fail at nodes', all(reports[i]['forced'].startswith('node_') for i in (0,3)))
    # The independent cyclic marking can permute the three discriminant
    # lines relative to the search's marking.
    check('mixed profiles lie in A2 discriminant', all(
          reports[i]['forced'].startswith('dp6_') and
          reports[i]['forced'].split(':')[1] in ('a','b','a-b') for i in (1,2)))
    check('ranks and kernel dimensions',
          [(r['rank'],r['kernel_dimension']) for r in reports] == [(19,13),(20,13),(20,13),(20,14)])
    full_labels = [f'n{i}' for i in range(30)]+['lambda0','a0','b0','lambda1','a1','b1']
    full = matrix(QQ,fixed+[r for block in blocks for r in block])
    universal = []
    for record_index,expected in (
            (0,{'n13':1,'a0':-1,'b0':1,'a1':-1,'b1':1}),
            (1,{'b0':1,'b1':1}),
            (3,{'n16':1,'lambda0':1,'lambda1':1})):
        D = vector(QQ,reports[record_index]['dual_divisor'])
        values = vector(QQ,[expected.get(label,0) for label in full_labels])
        check('universal divisor identity '+str(expected), full*D == values)
        universal.append(dict(divisor=serial(D),nonzero_pairings=expected))
    out = dict(vertices=V,facet_normals=[serial(n) for n in normals],
               lattice_points=[serial(p) for p in sorted(lattice_points)],
               f_vector=[int(x) for x in P.f_vector()],nodes=nodes,dp6_faces=sixes,
               smooth_two_faces=smooth,local_intersection_classes=[serial(a) for a in local_classes],
               h11=23,non_toric_correction=0,profiles=reports,
               universal_divisor_identities=universal,checks=checks,passed=len(checks),
               scope='exact inputs to the counterexample proof; analytic necessity is proved separately')
    target = Path(__file__).with_suffix('.json')
    target.write_text(json.dumps(out,indent=2)+'\n')
    print(f'PASS: {len(checks)} exact checks; 30 nodes, 2 dP6 points, h11=23')
    print('dp6 faces:',sixes)
    for r in reports:
        print(r['profile'],'rank',r['rank'],'kernel',r['kernel_dimension'],
              'weak',r['weak_pass'],'forced',r['forced'])
        print('  D:',r['dual_divisor'],'scale',r['dual_scale'])
        if 'weak_witness' in r:
            print('  witness:',r['weak_witness'])


if __name__ == '__main__':
    main()
