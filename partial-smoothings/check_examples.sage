"""Independent exact reconstruction for the examples in Partial smoothings.

Run from this directory with SageMath 10.9: sage -c "load('check_examples.sage')".
The input is the ordered vertices and cyclic hexagon markings. The relation
and divisor tests use local face circuits, with no stored charge matrices.
"""
from sage.all import Polyhedron, QQ, ZZ, gcd, matrix, vector

DELTA = [
    (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1), (1,-1,1,-1),
    (0,0,1,-1), (1,-1,0,0), (-1,1,0,0), (0,0,-1,1), (-1,1,-1,1),
    (0,0,0,-1), (0,0,-1,0), (-1,1,-1,0), (0,-1,1,0), (-1,0,1,0),
    (0,-1,0,1), (-1,0,0,1), (-1,0,1,-1), (0,-1,1,-1), (-1,0,-1,1),
    (0,-1,-1,1), (0,-1,-1,0), (0,-1,0,-1), (-1,0,0,-1), (-1,0,-1,0),
]
DELTA_MARKS = [(13,15,20,21,22,18), (14,16,19,24,23,17)]
CONTROL = [
    (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1), (-1,1,1,0),
    (1,-1,0,0), (1,-1,-1,1), (0,0,1,-1), (0,0,-1,1), (-1,1,1,-1),
    (-1,1,0,0), (1,-1,-1,0), (0,0,0,-1), (0,0,-1,0), (-1,1,0,-1),
    (0,-1,0,1), (-1,0,1,0), (1,-2,-1,1), (0,-1,-1,1), (-1,0,1,-1),
    (1,-2,-1,0), (-2,1,1,-1), (0,-1,0,-1), (-1,0,-1,0), (-2,1,0,-1),
    (0,-1,-1,-1), (-1,0,-1,-1),
]
CONTROL_MARKS = [(19,22,20,17,15,16), (18,15,16,21,24,23),
                 (25,26,23,18,17,20), (25,22,19,21,24,26)]


def generate(vertices, cycles):
    """Form the intersection vectors in the prime divisors of the refined fan."""
    points = [vector(ZZ, p) for p in vertices]
    centers = [vector(ZZ, sum((points[i] for i in face), vector(ZZ, 4))/6) for face in cycles]
    rays = points + centers
    poly = Polyhedron(vertices=vertices, base_ring=QQ)
    lookup = {tuple(p): i for i,p in enumerate(points)}
    all_faces = sorted(tuple(sorted(lookup[tuple(v)] for v in f.vertices())) for f in poly.faces(2))
    normals = [vector(QQ, h[1:])/h[0] for h in poly.inequalities()]
    assert all(h[0] > 0 for h in poly.inequalities()) and all(x in ZZ for v in normals for x in v)
    assert len(poly.vertices()) == len(vertices) and set(map(tuple, poly.integral_points())) == (
        set(map(tuple, rays)) | {(0,0,0,0)})
    assert set(f for f in all_faces if len(f) == 6) == {tuple(sorted(c)) for c in cycles}
    axes = matrix.identity(ZZ, len(rays))

    def unit_triangle(i,j,k):
        return matrix(ZZ,[rays[i],rays[j],rays[k]]).elementary_divisors() == [1,1,1]

    nodes = []
    for face in all_faces:
        if len(face) == 3:
            assert unit_triangle(*face)
            continue
        pair = [u for u in normals if all(u*points[i] == -1 for i in face)]
        assert len(pair) == 2 and gcd(vector(ZZ,pair[0]-pair[1])) == 1
        if len(face) == 4:
            i, *rest = face
            positive = [j for j in rest if 2*(points[i]+points[j]) == sum((points[k] for k in face),vector(ZZ,4))]
            j, = positive
            k,l = sorted(set(rest)-{j})
            assert unit_triangle(i,j,k)
            nodes.append(axes[i]+axes[j]-axes[k]-axes[l])
        else:
            assert len(face) == 6
    assert len(nodes) == 30
    boundaries = []
    for c,face in enumerate(cycles,len(vertices)):
        wall = []
        for j,x in enumerate(face):
            n,m = face[j-1],face[(j+1)%6]
            assert unit_triangle(c,x,m)
            wall.append(axes[n]+axes[m]-axes[x]-axes[c])
        boundaries.append(matrix(ZZ,wall))
    N = matrix(ZZ,nodes)
    assert all(B * matrix(ZZ,rays) == 0 for B in [N,*boundaries])
    return N, boundaries[:2]


def roots(B):
    return [B[0]-B[3], B[4]-B[0], B[2]-B[4]]


def divisor(width, indices, coefficient=1):
    return vector(ZZ,[coefficient*int(i in indices) for i in range(width)])


def verify_profile(N, surfaces, profile, witness, forced, rank, active_rank):
    C = N
    for B, mark in zip(surfaces,profile):
        C = C.stack(B.matrix_from_rows([0,3] if mark == "L" else [0,2,4]))
    K = C.left_kernel().basis_matrix()
    assert C.nrows()-K.nrows() == rank
    assert [i for i,column in enumerate(K.columns()) if not any(column)] == forced
    assert vector(ZZ,witness)*C == 0
    assert [i for i,z in enumerate(witness) if not z] == forced
    active = C.matrix_from_rows([i for i,z in enumerate(witness) if z])
    assert [e for e in active.elementary_divisors() if e] == [1]*active_rank
    print(profile, "columns",C.nrows(),"rank",rank,"relation dimension",K.nrows(),
          "forced",forced,"active Smith factors",[1]*active_rank)
    return C,active


d = [1,-2,1,1,1,-2,-1,1,-1,2,-2,2,1,1,-1,-1,-1,1,-1,1,1,1,-1,-2,2,-1,2,1,-2,2]
p = [1,2,1,-2,-1,-2,2,-1,2,-1,1,-1,1,1,-1,2,-2,0,-1,1,2,-1,2,-1,1,-1,1,-1,-1,1]
q = [-2,-1,1,1,-1,-1,2,1,-1,-1,-2,-2,2,-1,-1,1,1,-2,-1,1,1,-3,-1,2,-1,2,1,2,2,-2]

N, blocks = generate(DELTA, DELTA_MARKS)
C, _ = verify_profile(N,blocks,"LP",d+[1,-1,-1,0,1],[33],22,21)
assert C*divisor(C.ncols(), [2,3,13,14,15,16]) == divisor(35,[33])

N, blocks = generate(CONTROL, CONTROL_MARKS)
Nroots = matrix(ZZ,roots(blocks[0])+roots(blocks[1]))
C, _ = verify_profile(N,blocks,"LP",[x+3*y for x,y in zip(p,q)]+[2,-2,6,-9,3],[],22,22)
verify_profile(N,blocks,"LL",p+[2,-2,0,0],[17,32,33],23,20)
A = N.stack(Nroots.matrix_from_rows([0,3]))
assert A*divisor(A.ncols(),[2,3,4,15,16]) == divisor(32,[17])
assert A*divisor(A.ncols(),[3,6,8,15,17,18,23]) == divisor(32,[31])

C, active = verify_profile(N,blocks,"PP",q+[0,0,0,2,-3,1],[30,31,32],24,21)
generators = matrix(ZZ, [[0,1,0,0,0,0],[0,0,1,1,0,0]])
assert (active.right_kernel().basis_matrix()*Nroots.T).row_module().basis_matrix() == generators
D1 = divisor(C.ncols(),[11,13,20,22,25,26],-1)
D2 = divisor(C.ncols(),[3,6,8,15,17,18,23,28])
assert active*D1 == 0 and active*D2 == 0
assert Nroots*D1 == generators[0] and Nroots*D2 == generators[1]
assert active.rank() == N.rank()+1 and N.stack(Nroots).rank() == N.rank()+2
print("PP survivors on the six roots:", list(generators.rows()))

# Volumes determining the integral functional dual to H on the kernel.
fan_poly = Polyhedron(vertices=CONTROL, base_ring=QQ)
polar_vertices = [list(vector(QQ,h[1:])/h[0]) for h in fan_poly.inequalities()]
assert 24*Polyhedron(vertices=polar_vertices, base_ring=QQ).volume() == 83
# v0 is the first standard basis vector. Its dual facet is u0=-1.
facet = Polyhedron(vertices=[u[1:] for u in polar_vertices if u[0]==-1], base_ring=QQ)
assert 6*facet.volume() == 12 and 7*12-83 == 1
