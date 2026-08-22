# A fixed smooth projective crepant fan

**Status:** exact WP2/WP0-resolution input, 20 August 2026. This note fixes
one smooth projective crepant subdivision of the face fan of
$\Delta^\vee$. The rays, all maximal cones, the projectivity certificate,
and all 26 square diagonals are checked by
[`resolution_fan.sage`](resolution_fan.sage). The choice makes the
resolution-theoretic route concrete; it does not yet determine the global
deformation obstruction map.

## 1. Rays

Keep the ordering $q_1,\ldots,q_{26}$ of `example_data.md`, and set

$$
\begin{aligned}
q_{27}&=(-1,-1,0,0)=e_{6,1},&
q_{28}&=(-1,0,0,0)=e_{6,2},\\
q_{29}&=(0,0,-1,-1)=e_{7,1},&
q_{30}&=(0,0,-1,0)=e_{7,2}.
\end{aligned}
$$

The complete ray table is therefore:

| $i$ | $q_i$ | $i$ | $q_i$ |
|---:|---|---:|---|
| 1 | $(-1,-1,-1,-1)$ | 16 | $(0,-1,0,0)$ |
| 2 | $(-1,-1,-1,0)$ | 17 | $(0,1,-1,-1)$ |
| 3 | $(-1,-1,0,-1)$ | 18 | $(0,1,-1,0)$ |
| 4 | $(-1,-1,0,1)$ | 19 | $(0,1,0,-1)$ |
| 5 | $(-1,-1,1,0)$ | 20 | $(0,1,0,0)$ |
| 6 | $(-1,-1,1,1)$ | 21 | $(1,1,-1,-1)$ |
| 7 | $(-1,0,-1,-1)$ | 22 | $(1,0,-1,0)$ |
| 8 | $(-1,0,-1,0)$ | 23 | $(0,0,0,-1)$ |
| 9 | $(-1,0,0,-1)$ | 24 | $(0,0,0,1)$ |
| 10 | $(-1,0,0,1)$ | 25 | $(0,0,1,0)$ |
| 11 | $(-1,0,1,0)$ | 26 | $(0,0,1,1)$ |
| 12 | $(-1,0,1,1)$ | 27 | $(-1,-1,0,0)$ |
| 13 | $(0,-1,-1,-1)$ | 28 | $(-1,0,0,0)$ |
| 14 | $(0,-1,-1,0)$ | 29 | $(0,0,-1,-1)$ |
| 15 | $(0,-1,0,1)$ | 30 | $(0,0,-1,0)$ |

Exact lattice-point enumeration gives

$$
\Delta^\vee\cap M
=\{0,q_1,\ldots,q_{30}\}.
$$

Thus these are all nonzero boundary lattice points; there are no hidden
edge- or facet-interior rays.

## 2. Construction and projectivity certificate

Define a common integral perturbation on the 30 rays by

$$
b_i=
\begin{cases}
2^{26-i},&1\leq i\leq26,\\
-2^{30}+2^{30-i},&27\leq i\leq30.
\end{cases}
$$

On every three-dimensional facet of $\Delta^\vee$, take the lower facets
of the lifted points $(q_i,b_i)$. The checker verifies directly that these
are simplicial, that the induced triangulations agree on every intersection
of two facets, and that they use every boundary lattice point. Coning the
boundary tetrahedra to the origin therefore gives a complete fan
$\widehat\Sigma$ refining the face fan $\Sigma$.

The large negative term in $b_{27},\ldots,b_{30}$ forces the restriction to
each of the four divisorial polygons to be its full star triangulation at
the unique interior lattice point. Thus the exceptional surfaces on the
resolved hypersurface are the intended toric del Pezzo surfaces
$E_{6,1},E_{6,2}\simeq dP_6$ and
$E_{7,1},E_{7,2}\simeq dP_7$, rather than toric boundary blowdowns.

Projectivity is certified without a floating-point regularity test. Assign
the integral height

$$
a_i=2{,}181{,}038{,}071+b_i\qquad(1\leq i\leq30)
$$

to $q_i$. For every maximal cone $\sigma$ below, let $m_\sigma$ be the
unique rational linear form satisfying

$$
\langle m_\sigma,q_i\rangle=a_i
\qquad(q_i\in\sigma).
$$

The exact inequalities

$$
a_j-\langle m_\sigma,q_j\rangle>0
\qquad(q_j\notin\sigma)
$$

hold for all $124\cdot26$ outside-ray tests; the minimum margin is exactly
$1$. The checker also derives $2{,}181{,}038{,}071$ as the least positive
integer baseline for which this fixed perturbation passes every inequality.
Hence the $a_i$ define a strictly convex integral support function. The fan
is projective.

Every maximal cone has determinant $1$. Consequently
$Y_{\widehat\Sigma}$ is smooth. Since the only added rays are boundary
lattice points lying at height one in the corresponding original Gorenstein
cones, the toric morphism

$$
Y_{\widehat\Sigma}\longrightarrow Y_\Sigma
=\mathbb P_{\Delta^\vee}
$$

is crepant. It is a smooth MPCP resolution in this example; its generic
anticanonical hypersurface is the fixed $\widehat X^\circ$ used in
[`restriction_data.md`](restriction_data.md).

## 3. The 124 maximal cones

The following tuples give their one-based ray indices.

<!-- MAXIMAL_CONES_BEGIN -->
```python
MAXIMAL_CONES = (
    (1, 2, 7, 27),
    (1, 2, 7, 29),
    (1, 2, 13, 27),
    (1, 2, 13, 29),
    (1, 3, 7, 27),
    (1, 3, 7, 29),
    (1, 3, 13, 27),
    (1, 3, 13, 29),
    (2, 4, 8, 27),
    (2, 4, 8, 30),
    (2, 4, 14, 27),
    (2, 4, 14, 30),
    (2, 7, 8, 27),
    (2, 7, 8, 29),
    (2, 8, 29, 30),
    (2, 13, 14, 27),
    (2, 13, 14, 29),
    (2, 14, 29, 30),
    (3, 5, 9, 23),
    (3, 5, 9, 27),
    (3, 5, 13, 23),
    (3, 5, 13, 27),
    (3, 7, 9, 27),
    (3, 7, 9, 29),
    (3, 9, 23, 29),
    (3, 13, 23, 29),
    (4, 6, 10, 15),
    (4, 6, 10, 27),
    (4, 6, 15, 27),
    (4, 8, 10, 27),
    (4, 8, 10, 30),
    (4, 10, 15, 30),
    (4, 14, 15, 27),
    (4, 14, 15, 30),
    (5, 6, 11, 25),
    (5, 6, 11, 27),
    (5, 6, 16, 25),
    (5, 6, 16, 27),
    (5, 9, 11, 23),
    (5, 9, 11, 27),
    (5, 11, 23, 25),
    (5, 13, 16, 23),
    (5, 13, 16, 27),
    (5, 16, 23, 25),
    (6, 10, 12, 15),
    (6, 10, 12, 27),
    (6, 11, 12, 25),
    (6, 11, 12, 27),
    (6, 12, 15, 26),
    (6, 12, 25, 26),
    (6, 15, 16, 26),
    (6, 15, 16, 27),
    (6, 16, 25, 26),
    (7, 8, 17, 28),
    (7, 8, 17, 29),
    (7, 8, 27, 28),
    (7, 9, 17, 28),
    (7, 9, 17, 29),
    (7, 9, 27, 28),
    (8, 10, 18, 28),
    (8, 10, 18, 30),
    (8, 10, 27, 28),
    (8, 17, 18, 28),
    (8, 17, 18, 29),
    (8, 18, 29, 30),
    (9, 11, 19, 23),
    (9, 11, 19, 28),
    (9, 11, 27, 28),
    (9, 17, 19, 28),
    (9, 17, 19, 29),
    (9, 19, 23, 29),
    (10, 12, 15, 24),
    (10, 12, 18, 24),
    (10, 12, 18, 28),
    (10, 12, 27, 28),
    (10, 15, 24, 30),
    (10, 18, 24, 30),
    (11, 12, 19, 25),
    (11, 12, 19, 28),
    (11, 12, 27, 28),
    (11, 19, 23, 25),
    (12, 15, 24, 26),
    (12, 18, 20, 24),
    (12, 18, 20, 28),
    (12, 19, 20, 25),
    (12, 19, 20, 28),
    (12, 20, 24, 26),
    (12, 20, 25, 26),
    (13, 14, 16, 22),
    (13, 14, 16, 27),
    (13, 14, 22, 29),
    (13, 16, 22, 23),
    (13, 21, 22, 23),
    (13, 21, 22, 29),
    (13, 21, 23, 29),
    (14, 15, 16, 22),
    (14, 15, 16, 27),
    (14, 15, 22, 30),
    (14, 22, 29, 30),
    (15, 16, 22, 26),
    (15, 22, 24, 26),
    (15, 22, 24, 30),
    (16, 22, 23, 25),
    (16, 22, 25, 26),
    (17, 18, 19, 21),
    (17, 18, 19, 28),
    (17, 18, 21, 29),
    (17, 19, 21, 29),
    (18, 19, 20, 21),
    (18, 19, 20, 28),
    (18, 20, 21, 24),
    (18, 21, 22, 24),
    (18, 21, 22, 30),
    (18, 21, 29, 30),
    (18, 22, 24, 30),
    (19, 20, 21, 25),
    (19, 21, 23, 25),
    (19, 21, 23, 29),
    (20, 21, 24, 26),
    (20, 21, 25, 26),
    (21, 22, 23, 25),
    (21, 22, 24, 26),
    (21, 22, 25, 26),
    (21, 22, 29, 30),
)
```
<!-- MAXIMAL_CONES_END -->

## 4. Fixed diagonals at the 26 nodes

For a node face $N_i$, the following pair is the diagonal contained in the
two triangles of the chosen subdivision.

<!-- NODE_DIAGONALS_BEGIN -->
```python
NODE_DIAGONALS = (
    (2, 7),    # N_1
    (2, 13),   # N_2
    (3, 7),    # N_3
    (4, 8),    # N_4
    (4, 14),   # N_5
    (5, 9),    # N_6
    (5, 13),   # N_7
    (5, 23),   # N_8
    (9, 23),   # N_9
    (6, 10),   # N_10
    (10, 15),  # N_11
    (6, 11),   # N_12
    (6, 16),   # N_13
    (6, 25),   # N_14
    (8, 17),   # N_15
    (9, 17),   # N_16
    (12, 18),  # N_17
    (12, 24),  # N_18
    (12, 19),  # N_19
    (12, 25),  # N_20
    (14, 16),  # N_21
    (16, 23),  # N_22
    (16, 26),  # N_23
    (18, 19),  # N_24
    (20, 24),  # N_25
    (20, 25),  # N_26
)
```
<!-- NODE_DIAGONALS_END -->

These choices fix the signs of the exceptional curve classes once an
orientation convention is imposed. In particular, the two coloop nodes use
the diagonals $(q_9,q_{23})$ and $(q_{10},q_{15})$.

## 5. What is fixed and what remains

Let

$$
\pi:\widehat X^\circ\longrightarrow X^\circ
$$

be the restriction of the toric resolution to the proper transform of the
generic anticanonical hypersurface. The four new rays define the exceptional
surfaces $E_{6,1},E_{6,2},E_{7,1},E_{7,2}$, and the 26 diagonals define the
exceptional curves $\Gamma_1,\ldots,\Gamma_{26}$. The standard MPCP
construction gives

$$
(h^{1,1},h^{2,1})(\widehat X^\circ)=(26,20).
$$

The resolution is now fixed strongly enough to compute divisor relations,
curve classes, restrictions to the four exceptional surfaces, and the
resolution-theoretic boundary map. None of those maps is inferred merely
from the existence of the fan; they are the next Route C calculations.
