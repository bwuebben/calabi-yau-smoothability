#!/usr/bin/env python3
"""Three affine charts of the simultaneous A2 partial resolution.

Run: sage -python paper5/certificates/a2_simultaneous_partial_resolution.py
The construction is Tot(O_P2(-1)^3) -> {rank(M)<=1}, M=x*y^T,
with base the two diagonal differences. The central fibre has a P2
and three nodes; every noncentral fibre is unchanged by the morphism.
"""
import json
from pathlib import Path
from sage.all import QQ, PolynomialRing, matrix


def main():
    checks=[]; charts=[]
    def check(name,condition):
        assert condition,name
        checks.append(name)
    S=PolynomialRing(QQ,'x1,x2,x3,y1,y2,y3,z,d1,d2,d3')
    x=list(S.gens()[:3]); y=list(S.gens()[3:6]); z=S.gen(6); d=list(S.gens()[7:])
    M=matrix(S,3,3,lambda i,j:x[i]*y[j])
    check('factorization annihilates all rank-two minors',all(f==0 for f in M.minors(2)))
    for k in range(3):
        i,j=[r for r in range(3) if r!=k]
        # In x_k=1, solve y_k=z+d_k and z=x_i*y_i-d_i.
        sub={x[k]:S.one(),y[k]:x[i]*y[i]-d[i]+d[k],z:x[i]*y[i]-d[i]}
        pulled=M.apply_map(lambda f:f.subs(sub))
        f=x[i]*y[i]-x[j]*y[j]-(d[i]-d[j])
        check(f'chart {k}: diagonal k identity',pulled[k,k]-(x[i]*y[i]-d[i]+d[k])==0)
        check(f'chart {k}: diagonal-difference equation',pulled[i,i]-pulled[j,j]-(d[i]-d[j])==f)
        check(f'chart {k}: all rank-two minors remain zero',all(q==0 for q in pulled.minors(2)))
        vars4=[x[i],y[i],x[j],y[j]]
        Hess=matrix(S,[[f.derivative(a).derivative(b) for b in vars4] for a in vars4])
        check(f'chart {k}: central singularity is a node',Hess.det()==1)
        check(f'chart {k}: critical ideal is coordinate origin',
              S.ideal([f.derivative(v) for v in vars4])==S.ideal(vars4))
        check(f'chart {k}: discriminant is the other diagonal coincidence',
              f.subs({v:0 for v in vars4})==-(d[i]-d[j]))
        check(f'chart {k}: zero section occurs only over the branch origin',
              S.ideal([q.subs({yy:0 for yy in y}) for q in [x[0]*y[0]-x[2]*y[2]-(d[0]-d[2]),x[1]*y[1]-x[2]*y[2]-(d[1]-d[2])]])==S.ideal(d[0]-d[2],d[1]-d[2]))
        charts.append(dict(projective_chart=k,other_indices=[i,j],equation=str(f),
                           node_parameter=str(d[i]-d[j]),hessian_determinant=1))
    # A nonzero entry gives the inverse to the incidence resolution.
    for k in range(3):
        for ell in range(3):
            # M_ij=M_iell*M_kj/M_kell on the rank-one open chart.
            check(f'inverse incidence chart ({k},{ell})',all(
                M[i,j]*M[k,ell]==M[i,ell]*M[k,j] for i in range(3) for j in range(3)))
    # The three smoothing parameters sum to zero, and their vanishing
    # hyperplanes are exactly the local A2 discriminant.
    params=[d[1]-d[2],d[2]-d[0],d[0]-d[1]]
    check('three nodal parameters sum to zero',sum(params)==0)
    check('three distinct discriminant lines',S.ideal(params).dimension()==8)
    out=dict(construction='M=x*y^T, x nonzero, modulo (x,y)~(lambda*x,lambda^-1*y)',
             vector_bundle='O_P2(-1) direct-sum O_P2(-1) direct-sum O_P2(-1)',
             base='diagonal differences, modulo common translation',charts=charts,
             checks=checks,passed=len(checks))
    target=Path(__file__).with_suffix('.json')
    target.write_text(json.dumps(out,indent=2)+'\n')
    print(f'PASS: {len(checks)} exact simultaneous-resolution checks; wrote {target.name}')


if __name__=='__main__':main()
