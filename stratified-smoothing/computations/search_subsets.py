"""Discovery search, in SageMath: sage -python search_subsets.py.

Read the saved Paper V census. Inspect four entries with few or no old-node
relations. For 248, 439 and 484, try all subsets of at least two hexagons.
For the exceptional twenty-hexagon entry 575, try all subsets of size 2--5.
The vertex-only verifier check_models.py independently certifies the examples.
"""

from itertools import combinations, combinations_with_replacement, product
import json
from pathlib import Path

from sage.all import QQ, VectorSpace, matrix, vector

HERE = Path(__file__).resolve().parent
INPUT = HERE / 'census_search_inputs.json'


def select(m, faces):
    return m.matrix_from_rows([3*i+j for i in faces for j in range(3)])


def profiles(k, count):
    for marks in product('LP',repeat=count):
        zero = [3*i+j for i,m in enumerate(marks) for j in ([1,2] if m=='L' else [0])]
        yield k.matrix_from_columns(zero).left_kernel().basis_matrix()


def quadratic_rows(k, count):
    return matrix(QQ,[
        [k[a,3*i]*k[b,3*i+j] + (a!=b)*k[b,3*i]*k[a,3*i+j]
         for a,b in combinations_with_replacement(range(k.nrows()),2)]
        for i in range(count) for j in (1,2)
    ])


def search(record):
    number = record['census_index']
    n,seven,h = record['inventory']
    assert not seven
    data = record if 'full_branch_matrix' in record else record['uniform_profile_tests']
    # Quotient the root-row space by the node rows, retaining the full
    # cone image of the unrestricted relation space.
    rows = matrix(QQ,data['full_branch_matrix'])
    nodes = rows[:n,:]
    roots = rows[n:,:] * nodes.right_kernel().basis_matrix().T
    hits = {}
    inspected = 0
    for size in range(2,(5 if number==575 else h)+1):
        for faces in combinations(range(h),size):
            inspected += 1
            k = select(roots,faces).left_kernel().basis_matrix()
            if k.nrows()!=2:
                continue
            qrank = quadratic_rows(k,size).rank()
            if qrank < 2:
                continue
            if any(s.nrows() for s in profiles(k,size)):
                continue
            hits.setdefault(qrank,[0,faces])[0] += 1
    print(f'{number}: {inspected} subsets; pure node relations {n-nodes.rank()}; '
          f'two cone directions, zero selected cone image (count, first subset) {hits}')
    return roots


def failed_master(roots):
    master = list(range(7))
    target = [0,1,2,5,6]
    k = select(roots,master).left_kernel().basis_matrix()
    killed = [3*i+j for i in master if i not in target for j in range(3)]
    a,b = k.matrix_from_columns(killed).left_kernel().basis()
    plane = matrix(QQ,[a,b])*k
    functional, = quadratic_rows(plane,len(master)).right_kernel().basis()

    def sym(x,y):
        return vector(QQ,[(x[i]*y[j]+y[i]*x[j])/2
                          for i,j in combinations_with_replacement(range(len(x)),2)])
    mixed = functional[0]*sym(a,a) + 2*functional[1]*sym(a,b) + functional[2]*sym(b,b)
    tangent_vectors,quadratic_vectors = [],[]
    for s in profiles(k,len(master)):
        tangent_vectors.extend(s.rows())
        quadratic_vectors += [sym(s[i],s[j]) for i,j in combinations_with_replacement(range(s.nrows()),2)]
    tensors = VectorSpace(QQ, len(mixed)).subspace(quadratic_vectors)
    assert matrix(QQ,tangent_vectors).rank()==8
    assert len(mixed)==36 and tensors.dimension()==22 and mixed not in tensors
    print('seven-face master, faces 0..6: tangent span 8, quadratic span 22/36;')
    print('  target (0,1,2,5,6) missing socle tensor lies outside that span')


def main():
    data = json.loads(INPUT.read_text())
    for record in data['records'] + data['skipped']:
        if record['census_index'] in (248,439,484,575):
            roots = search(record)
            if record['census_index']==575:
                failed_master(roots)


if __name__ == '__main__':
    main()
