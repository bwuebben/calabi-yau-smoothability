#!/usr/bin/env python3
"""Standard-library rational replay of the counterexample's linear inputs.

The polytope geometry is checked separately by the independent Sage script.
Here all matrix products and ranks use fractions.Fraction, without Sage.
"""
from fractions import Fraction as Q
import json
from pathlib import Path


def dot(a,b):
    return sum((x*y for x,y in zip(a,b)), Q(0))


def matvec(A,v):
    return [dot(r,v) for r in A]


def leftvec(v,A):
    return [dot(v,c) for c in zip(*A)]


def rank(A):
    A = [r[:] for r in A]
    k = 0
    for j in range(len(A[0])):
        pivot = next((i for i in range(k,len(A)) if A[i][j]),None)
        if pivot is None:
            continue
        A[k],A[pivot] = A[pivot],A[k]
        scale = A[k][j]
        A[k] = [x/scale for x in A[k]]
        for i in range(k+1,len(A)):
            scale = A[i][j]
            if scale:
                A[i] = [x-scale*y for x,y in zip(A[i],A[k])]
        k += 1
    return k


def main():
    data = json.loads(Path(__file__).with_name('a2_branch_selection_counterexample.json').read_text())
    passed = 0
    def check(condition):
        nonlocal passed
        assert condition
        passed += 1
    matrices = {}
    for rec in data['profiles']:
        A = [[Q(x) for x in r] for r in rec['matrix']]
        matrices[rec['profile']] = A
        check(rank(A) == rec['rank'])
        check(len(A)-rank(A) == rec['kernel_dimension'])
        check(matvec(A,list(map(Q,rec['dual_divisor']))) ==
              [Q(rec['dual_scale'])*Q(x) for x in rec['forced_functional']])
        check(all(dot(r,list(col)) == 0 for r in A for col in zip(*data['vertices'])))
        if rec['weak_pass']:
            w = list(map(Q,rec['weak_witness']))
            check(not any(leftvec(w,A)))
            check(all(w[:30]))
            for i,ch in enumerate(rec['profile']):
                labels = rec['labels']
                if ch == 'L':
                    check(w[labels.index(f'dp6_{i}:L')] != 0)
                else:
                    a = w[labels.index(f'dp6_{i}:a')]
                    b = w[labels.index(f'dp6_{i}:b')]
                    check(bool(a or b) and not a*b*(a-b))
        part = rec['partial_model']
        N = [[Q(x) for x in r] for r in part['matrix']]
        check(matvec(N,list(map(Q,part['dual_divisor']))) ==
              [Q(part['dual_scale'])*Q(x) for x in part['functional']])
        check(len(N)-rank(N) == part['kernel_dimension'])
        check('node_' in part['certified_node'])
    full = (matrices['LL'][:30]+[matrices['LL'][30]]+
            matrices['PL'][30:32]+[matrices['LL'][31]]+matrices['LP'][31:33])
    labels = [f'n{i}' for i in range(30)]+['lambda0','a0','b0','lambda1','a1','b1']
    # Unrestricted tangent calculation: the cone image has dimension five,
    # and the separately checked H_M identity puts it in b0+b1=0.
    full_rank = rank(full)
    node_rank = rank(full[:30])
    check(full_rank == 20)
    check(len(full)-full_rank == 16)
    check(node_rank == 19)
    check(30-node_rank == 11)
    check((len(full)-full_rank)-(30-node_rank) == 5)
    for profile,A in matrices.items():
        ids = list(range(30))
        for i,branch in enumerate(profile):
            ids += [30+3*i] if branch == 'L' else [31+3*i,32+3*i]
        check([full[j] for j in ids] == A)
    for identity in data['universal_divisor_identities']:
        check(matvec(full,list(map(Q,identity['divisor']))) ==
              [Q(identity['nonzero_pairings'].get(label,0)) for label in labels])
    # Independent monomial intersection: lcm is union of supports for these
    # squarefree monomials. Discard nonminimal supports after each product.
    primes = [set('yvb'),set('ybu'),set('xvb'),set('xu')]
    generators = {frozenset()}
    for prime in primes:
        candidates = {support | {variable} for support in generators for variable in prime}
        generators = {support for support in candidates
                      if not any(other < support for other in candidates)}
    check(generators == {frozenset(s) for s in ['xy','xb','uv','ub']})
    pair_dimensions = [37-len(primes[i] | primes[j])
                       for i in range(4) for j in range(i+1,4)]
    check(pair_dimensions == [33,33,32,32,33,33])
    print(f'PASS: {passed} exact rational replay checks, no Sage')


if __name__ == '__main__':
    main()
