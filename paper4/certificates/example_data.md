# The fixed Paper 4 example

**Status:** WP0 exact-data record, 20 August 2026.

This file fixes the notation and finite geometric input for Paper 4. Every
enumerative statement below is asserted by
[`fixed_example.py`](fixed_example.py), which uses exact integer or rational
arithmetic. The source theorem is [Paper 3, Section
6](../../paper3/main.tex).

## 1. Lattice and polarity conventions

Let

$$
N=\mathbb Z^4,\qquad M=\operatorname{Hom}(N,\mathbb Z).
$$

For a reflexive polytope $\Delta\subset N_{\mathbb R}$, the papers use

$$
\Delta^\circ
=\{u\in M_{\mathbb R}:\langle u,v\rangle\geq -1
\text{ for every }v\in\Delta\}.
$$

The exact polytope routines use outward normals satisfying
$\langle u,v\rangle\leq1$. Thus the vertices in the paper's polar
convention are the negatives of the primitive facet normals returned by the
code. The checker verifies polarity in both directions.

Write

$$
\Delta=\Delta_{\mathbb F_1},\qquad
\Delta^\vee=\Delta_{\mathbb F_1}^{\circ}.
$$

The symbol $\Delta^\vee$ is used in this working file to keep the polar
polytope visually distinct from the Calabi–Yau $X^\circ$.

## 2. The 22 vertices of $\Delta$

The following order agrees with Paper 3.

| $i$ | $v_i$ | $i$ | $v_i$ |
|---:|---|---:|---|
| 1 | $(1,0,0,0)$ | 12 | $(-1,1,-1,1)$ |
| 2 | $(0,1,0,0)$ | 13 | $(0,-1,0,-1)$ |
| 3 | $(1,-1,0,0)$ | 14 | $(0,-1,-1,0)$ |
| 4 | $(0,0,1,0)$ | 15 | $(-1,1,-1,0)$ |
| 5 | $(0,0,0,1)$ | 16 | $(-1,0,0,-1)$ |
| 6 | $(0,0,1,-1)$ | 17 | $(-1,0,-1,1)$ |
| 7 | $(0,0,-1,1)$ | 18 | $(-1,-1,0,-1)$ |
| 8 | $(0,0,0,-1)$ | 19 | $(-2,1,-1,1)$ |
| 9 | $(0,0,-1,0)$ | 20 | $(-2,1,-1,0)$ |
| 10 | $(-1,1,0,0)$ | 21 | $(-1,-1,-1,0)$ |
| 11 | $(0,-1,0,0)$ | 22 | $(-2,0,-1,0)$ |

Exact consequences:

- $\Delta$ is reflexive, with 22 vertices and 26 facets.
- It has 72 two-faces: 50 smooth triangles, 20 unit squares, one
  $dP_6$ hexagon, and one $\mathbb F_1$ quadrilateral.
- Every dual edge has lattice length one.
- A generic anticanonical $X\subset\mathbb P_\Delta$ has 20 nodes, one
  $dP_6$-cone point, and one nonsmoothable $\mathbb F_1$-cone point.

The last item explains the notation but is not the deformation problem of
Paper 4.

## 3. The 26 vertices of $\Delta^\vee$

The order below is the authoritative Paper 4 ordering
$q_1,\ldots,q_{26}$. Later matrices and face labels must use it.

| $j$ | $q_j$ | $j$ | $q_j$ |
|---:|---|---:|---|
| 1 | $(-1,-1,-1,-1)$ | 14 | $(0,-1,-1,0)$ |
| 2 | $(-1,-1,-1,0)$ | 15 | $(0,-1,0,1)$ |
| 3 | $(-1,-1,0,-1)$ | 16 | $(0,-1,0,0)$ |
| 4 | $(-1,-1,0,1)$ | 17 | $(0,1,-1,-1)$ |
| 5 | $(-1,-1,1,0)$ | 18 | $(0,1,-1,0)$ |
| 6 | $(-1,-1,1,1)$ | 19 | $(0,1,0,-1)$ |
| 7 | $(-1,0,-1,-1)$ | 20 | $(0,1,0,0)$ |
| 8 | $(-1,0,-1,0)$ | 21 | $(1,1,-1,-1)$ |
| 9 | $(-1,0,0,-1)$ | 22 | $(1,0,-1,0)$ |
| 10 | $(-1,0,0,1)$ | 23 | $(0,0,0,-1)$ |
| 11 | $(-1,0,1,0)$ | 24 | $(0,0,0,1)$ |
| 12 | $(-1,0,1,1)$ | 25 | $(0,0,1,0)$ |
| 13 | $(0,-1,-1,-1)$ | 26 | $(0,0,1,1)$ |

Exact consequences:

- $\Delta^\vee$ is reflexive, with 26 vertices and 22 facets.
- It has 68 two-faces: 38 smooth triangles, 26 unit squares, two
  $dP_7$ pentagons, and two $dP_6$ hexagons.
- Every dual edge has lattice length one.

Let

$$
X^\circ\subset\mathbb P_{\Delta^\vee}
$$

be a generic anticanonical hypersurface. Since every singular face has
dual-edge length one, it contributes one singular point. Hence

$$
\operatorname{Sing}(X^\circ)
=\{N_1,\ldots,N_{26},D_{7,1},D_{7,2},D_{6,1},D_{6,2}\}.
$$

Here the $N_i$ are nodes and $D_{k,a}$ denotes both a face and, when no
confusion is possible, the corresponding cone singularity of degree $k$.
The generic defining equation is fixed in
[`global_section.md`](global_section.md). That note gives its 20 normalized
coefficients and the exact point of $X^\circ$ on each of these 30
one-dimensional torus orbits.

## 4. The four divisorial faces

| Label | Vertices | Interior lattice point | Local germ | Local smoothing components |
|---|---|---|---|---:|
| $D_{6,1}$ | $q_1,q_2,q_3,q_4,q_5,q_6$ | $(-1,-1,0,0)$ | anticanonical cone over $dP_6$ | 2 |
| $D_{6,2}$ | $q_7,q_8,q_9,q_{10},q_{11},q_{12}$ | $(-1,0,0,0)$ | anticanonical cone over $dP_6$ | 2 |
| $D_{7,1}$ | $q_1,q_7,q_{13},q_{17},q_{21}$ | $(0,0,-1,-1)$ | anticanonical cone over $dP_7$ | 1 |
| $D_{7,2}$ | $q_2,q_8,q_{14},q_{18},q_{22}$ | $(0,0,-1,0)$ | anticanonical cone over $dP_7$ | 1 |

The checker enumerates every lattice point on each face and verifies that the
displayed point is the unique nonvertex lattice point. Star subdivision at
these four points produces exceptional divisors that will be denoted

$$
E_{6,1},E_{6,2},E_{7,1},E_{7,2}.
$$

This notation does not yet choose a global MPCP subdivision.

The induced edge lattices of these four faces are matched to Altmann's
standard $dP_6$ and $dP_7$ polygons by explicit unimodular matrices in
[`local_deformations.md`](local_deformations.md). The same matches are exact
assertions in `fixed_example.py`. Their ambient Hilbert-character pullbacks
are checked by [`chart_characters.sage`](chart_characters.sage), while the
global Cox grading and all pairwise singular-chart intersections are computed
by [`cox_data.sage`](cox_data.sage).

## 5. The 26 nodal faces

The node labels are the lexicographic ordering of their one-based polar
vertex sets. For a square $N_i=(q_a,q_b,q_c,q_d)$, the checker finds the
unique diagonal relation

$$
\rho_i=\epsilon_a e_a+\epsilon_b e_b+\epsilon_c e_c+\epsilon_d e_d
\in\mathbb Z^{26},\qquad \epsilon_\bullet\in\{\pm1\}.
$$

The sign of an entire row is conventional. The table fixes one convention.

| Node | Square vertices | Positive indices in $\rho_i$ | Negative indices in $\rho_i$ |
|---:|---|---|---|
| $N_1$ | $1,2,7,8$ | $1,8$ | $2,7$ |
| $N_2$ | $1,2,13,14$ | $1,14$ | $2,13$ |
| $N_3$ | $1,3,7,9$ | $1,9$ | $3,7$ |
| $N_4$ | $2,4,8,10$ | $2,10$ | $4,8$ |
| $N_5$ | $2,4,14,15$ | $2,15$ | $4,14$ |
| $N_6$ | $3,5,9,11$ | $3,11$ | $5,9$ |
| $N_7$ | $3,5,13,16$ | $3,16$ | $5,13$ |
| $N_8$ | $3,5,23,25$ | $3,25$ | $5,23$ |
| $N_9$ | $3,9,19,23$ | $3,19$ | $9,23$ |
| $N_{10}$ | $4,6,10,12$ | $4,12$ | $6,10$ |
| $N_{11}$ | $4,10,15,24$ | $4,24$ | $10,15$ |
| $N_{12}$ | $5,6,11,12$ | $5,12$ | $6,11$ |
| $N_{13}$ | $5,6,15,16$ | $5,15$ | $6,16$ |
| $N_{14}$ | $5,6,25,26$ | $5,26$ | $6,25$ |
| $N_{15}$ | $7,8,17,18$ | $7,18$ | $8,17$ |
| $N_{16}$ | $7,9,17,19$ | $7,19$ | $9,17$ |
| $N_{17}$ | $10,12,18,20$ | $10,20$ | $12,18$ |
| $N_{18}$ | $10,12,24,26$ | $10,26$ | $12,24$ |
| $N_{19}$ | $11,12,19,20$ | $11,20$ | $12,19$ |
| $N_{20}$ | $11,12,25,26$ | $11,26$ | $12,25$ |
| $N_{21}$ | $13,14,15,16$ | $13,15$ | $14,16$ |
| $N_{22}$ | $13,16,23,25$ | $13,25$ | $16,23$ |
| $N_{23}$ | $15,16,25,26$ | $15,25$ | $16,26$ |
| $N_{24}$ | $17,18,19,20$ | $17,20$ | $18,19$ |
| $N_{25}$ | $18,20,24,26$ | $18,26$ | $20,24$ |
| $N_{26}$ | $19,20,25,26$ | $19,26$ | $20,25$ |

Let $R$ be the $26\times26$ matrix with rows $\rho_i$. Exact rational
row reduction gives

$$
\operatorname{rank}_{\mathbb Q}R=18.
$$

Deleting row 9 or row 11 lowers the rank to 17; deleting both lowers it to
16; deleting any other single row leaves rank 18. Therefore $\rho_9$ and
$\rho_{11}$ are the two coloops. Equivalently, every relation

$$
\sum_{i=1}^{26}\lambda_i\rho_i=0
$$

has $\lambda_9=\lambda_{11}=0$. In particular, the nodal subsystem has no
relation with every coefficient nonzero.

This proves only that a purely nodal Friedman–Batyrev–Kreuzer rescue is
impossible. It does not prove that $X^\circ$ is nonsmoothable.

## 6. Resolution notation

[`resolution_data.md`](resolution_data.md) now fixes a smooth projective
crepant 30-ray fan $\widehat\Sigma$. Write

- $Y_{\widehat\Sigma}\to\mathbb P_{\Delta^\vee}$ for the toric resolution;
- $\widehat X^\circ\to X^\circ$ for the induced crepant resolution of the
  hypersurface;
- $\Gamma_i$ for the exceptional curve over $N_i$; and
- $E_{6,1},E_{6,2},E_{7,1},E_{7,2}$ for the four exceptional del Pezzo
  surfaces.

All 26 square diagonals are fixed in that note. The primitive circuit row is
oriented so that the two rays of the chosen diagonal have coefficient $-1$;
this is the divisor-intersection row of $\Gamma_i$. The full-star
restrictions on the four divisorial polygons and all numerical restriction
maps are recorded in [`restriction_data.md`](restriction_data.md).

For any MPCP model used in Paper 3,

$$
(h^{1,1},h^{2,1})(\widehat X^\circ)=(26,20),
\qquad
(h^{1,1},h^{2,1})(\widehat X)=(20,26).
$$

## 7. Unconditional and conditional statements

The following are unconditional finite consequences of the printed vertices:

- reflexivity and polarity of the displayed pair;
- both face inventories;
- all singularity types and multiplicities;
- the four interior lattice points;
- the Hodge-number calculation;
- the nodal rank and coloop statements; and
- local smoothability and the number of local smoothing components.

Only the statement that this pair is unique in the complete
Kreuzer–Skarke classification depends on the database hypothesis. Paper 4's
deformation problem does not require that uniqueness claim.

## 8. Reproduction

From `cy_smoothing/`:

```bash
python3 paper4/certificates/fixed_example.py
python3 src/paper3_node_relations.py
python3 src/face_data.py
```

The first command is the Paper 4 gate. The other two are independent
cross-checks inherited from Paper 3.

## 9. Frozen notation rule

The labels $q_j$, $N_i$, $D_{k,a}$, $E_{k,a}$, and $\rho_i$ are
now frozen. Any later change must be accompanied by:

1. an update to `fixed_example.py`;
2. a successful exact rerun;
3. an update to every Paper 4 research note using the labels; and
4. an explanation in the corresponding data note of why the change was necessary.
