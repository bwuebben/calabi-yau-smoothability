# The mixed obstruction matrix

**Current reference status, 5 September 2026.** The finite rows and ranks
below are the exact data established on 20 August 2026.
[`mixed_candidate.sage`](mixed_candidate.sage) derives them from the
polygon coordinates and [`restriction_data.md`](restriction_data.md).
Their intrinsic interpretation is in [Paper 4, Sections 8.2–8.3](../main.tex)
and [Paper 5, Section 2](../../paper5/main.tex).
Paper 4's Sections 3–6 prove that the distinguished threefold admits no
smoothing, including arcs with higher-order contact. The earlier
first-order limitation of this computation is explained below; it is no
longer an unresolved smoothability question for that example.

## 1. From Altmann edge parameters to $K_E^\perp$

Let $v_1,\ldots,v_k$ be the cyclic vertices of a unit-edge polygon and put

$$
d_i=v_{i+1}-v_i.
$$

Altmann's integral tangent lattice is

$$
T_Q=
\left\{t=(t_i)\in\mathbb Z^k:
\sum_i t_i d_i=0\right\}
\big/\mathbb Z(1,\ldots,1).
$$

Define the cyclic difference

$$
\lambda_i=t_i-t_{i-1}.
$$

Then

$$
\sum_i\lambda_i=0,
\qquad
\sum_i\lambda_i v_i
=-\sum_i t_i d_i=0.
$$

If $e$ is the interior point and $u_i=v_i-e$ are the fan rays of the
exceptional surface $E$, these equations say

$$
\sum_i\lambda_i u_i=0,
\qquad
\lambda(K_E)=-\sum_i\lambda_i=0.
$$

Dualizing the toric divisor sequence therefore identifies $\lambda$ with a
functional in

$$
K_E^\perp\subset\operatorname{Pic}(E)^\vee.
$$

The cyclic-difference map kills exactly the simultaneous-translation
direction, and its image is the full saturated lattice $K_E^\perp$. Thus it
gives an integral isomorphism

$$
T_Q\simeq K_E^\perp.
$$

This is the local numerical bridge used below. The checker verifies all
edge equations, both vanishing sums, the Picard functional equations, ranks,
and saturation directly.

## 2. The ten divisorial rows

The $t$-vectors are the derivatives of the formulas for $t_i$ in
[`local_deformations.md`](local_deformations.md). The following table gives
the cyclic difference and the resulting functional in the Picard bases of
`restriction_data.md`.

| local coordinate | $t$ in canonical edge order | $\lambda$ in canonical vertex order | row in $\operatorname{Pic}(E)^\vee$ |
|---|---|---|---|
| $D_{6,1}:s_1$ | $(0,0,-1,0,0,-1)$ | $(1,0,-1,1,0,-1)$ | $(0,0,-1,1)$ |
| $D_{6,1}:s_2$ | $(0,-1,0,0,-1,0)$ | $(0,-1,1,0,-1,1)$ | $(-1,-1,1,0)$ |
| $D_{6,1}:s_3$ | $(0,-1,0,-1,0,-1)$ | $(1,-1,1,-1,1,-1)$ | $(-1,1,1,-1)$ |
| $D_{6,2}:s_1$ | $(0,0,-1,0,0,-1)$ | $(1,0,-1,1,0,-1)$ | $(0,0,-1,1)$ |
| $D_{6,2}:s_2$ | $(0,-1,0,0,-1,0)$ | $(0,-1,1,0,-1,1)$ | $(-1,-1,1,0)$ |
| $D_{6,2}:s_3$ | $(0,-1,0,-1,0,-1)$ | $(1,-1,1,-1,1,-1)$ | $(-1,1,1,-1)$ |
| $D_{7,1}:\alpha$ | $(0,0,-1,1,-1)$ | $(1,0,-1,2,-2)$ | $(1,-1,0)$ |
| $D_{7,1}:\beta$ | $(0,-1,0,0,-1)$ | $(1,-1,1,0,-1)$ | $(1,1,-1)$ |
| $D_{7,2}:\alpha$ | $(0,0,-1,1,-1)$ | $(1,0,-1,2,-2)$ | $(-2,0,1)$ |
| $D_{7,2}:\beta$ | $(0,-1,0,0,-1)$ | $(1,-1,1,0,-1)$ | $(-1,-1,1)$ |

For each row $w\in\operatorname{Pic}(E)^\vee$, compose with the exact
restriction matrix

$$
R_E:\operatorname{Pic}(Y_{\widehat\Sigma})
\longrightarrow\operatorname{Pic}(E).
$$

This gives ten rows $wR_E$ in the dual of the common rank-26 ambient Picard
lattice. Their span has rank $10$ and is saturated.

## 3. The candidate relation matrix

Put the 26 chosen node-curve rows first and the ten divisorial rows second:

$$
B_{\mathrm{mix}}=
\begin{pmatrix}
[C_1]\\[-2pt]
\vdots\\[-2pt]
[C_{26}]\\
w_{6,1,s_1}R_{6,1}\\[-2pt]
\vdots\\[-2pt]
w_{7,2,\beta}R_{7,2}
\end{pmatrix}
\in\operatorname{Mat}_{36\times26}(\mathbb Z).
$$

The matrix is completely determined by the tables in this note and
`restriction_data.md`. Exact Smith reduction gives

$$
\operatorname{rank}B_{\mathrm{mix}}=21,
\qquad
\dim_{\mathbb Q}\ker(B_{\mathrm{mix}}^{\mathsf T})=15,
$$

and all 21 nonzero Smith invariants are $1$. Among the full set of 36 rows,
the single coloop is

$$
D_{7,1}:\beta.
$$

Thus even before the nonlinear local branch equations are imposed, every
integral relation among all 36 candidate classes has zero coefficient in the
first $dP_7$ smoothing direction. This asymmetry belongs to the fixed fan and
labeling. The bare coloop label is not invariant under every allowed $dP_7$
change of coordinates, but its forced vanishing on the reduced smoothing
line is invariant by the symmetry comparison below.

## 4. Testing the four reduced smoothing profiles

Write a candidate local vector as

$$
(u_1,\ldots,u_{26};
s_{1,1},s_{2,1},s_{3,1};
s_{1,2},s_{2,2},s_{3,2};
\alpha_1,\beta_1;\alpha_2,\beta_2).
$$

For a reduced arc, $\alpha_1=\alpha_2=0$. On a $dP_6$ line component $L$
one imposes $s_1=s_2=0$ and requires $s_3\neq0$. On a plane component $P$
one imposes $s_3=0$ and requires
$s_1s_2(s_1-s_2)\neq0$. Every node coordinate and both $\beta$ coordinates
must also be nonzero for a transverse first-order simultaneous smoothing.

For each profile, the checker sets the prohibited coordinates to zero,
computes the exact kernel of the restricted transpose matrix, and tests
whether each required coordinate or difference vanishes identically on that
kernel.

| profile at $(D_{6,1},D_{6,2})$ | allowed coordinates | kernel dimension | required quantities forced to zero |
|---|---:|---:|---|
| $LL$ | 30 | 9 | $u_9,u_{11},\beta_1,\beta_2$ |
| $LP$ | 31 | 10 | $u_9,u_{11},\beta_1,\beta_2,s_{3,1}$ |
| $PL$ | 31 | 10 | $u_9,u_{11},\beta_1,\beta_2,s_{3,2}$ |
| $PP$ | 32 | 12 | $u_9,u_{11},\beta_1,\beta_2$ |

Consequently,

> No vector in $\ker(B_{\mathrm{mix}}^{\mathsf T})$ lies in any of the four
> transverse reduced local smoothing profiles.

The calculation also explains why the unrestricted tangent space is
misleading. Before one imposes the reduced Kuranishi branches, incompatible
$dP_6$ tangent directions can participate in relations that remove three of
the eventual forced-zero conditions. Those tangent vectors do not belong to
one reduced deformation component.

## 5. Intrinsic interpretation and the remaining formal task

Friedman--Laza's local theorem, the link Gysin sequence, and the exact
dimension match prove that, at every $dP_6$ and $dP_7$ point, the canonical
quotient

$$
T^1_p\longrightarrow H_2(L_p,\mathbb C)
$$

is an isomorphism. Their global local-cohomology diagram then identifies the
image of tangent localization with the kernel of the map from all 30 link
spaces to $H_2(\widehat X,\mathbb C)$. Section 5 of
`restriction_data.md` proves that the ambient Picard basis is the full
complex Picard space of $\widehat X$. These results establish the intrinsic
topological target of the matrix; the current global comparison is in Paper 4, Sections 8.2–8.3,
and Paper 5, Section 2.

The two local isomorphisms need not agree in the printed coordinates. Put
$C_p=D_p^{-1}\kappa_p$, where $D_p$ is cyclic difference and $\kappa_p$ is
the canonical link map. The local marking comparison in Paper 5, Section 2, identifies the
branch subspaces and their discriminants; the earlier coordinate argument
used this invariance of every $C_p$. If $C=\bigoplus_p C_p$, then
globally induced local vectors satisfy

$$
B_{\mathrm{mix}}^{\mathsf T}Cx=0.
$$

Since $C$ preserves the union of the four transverse smoothing profiles, the
forced-zero calculation is an unconditional first-order theorem: no image
vector is transverse to all 30 reduced local smoothing discriminants. Exact
LOCAL-NORM would identify the displayed coordinates more rigidly, but is not
needed for this conclusion.

The tangent calculation alone excludes transverse first-order smoothings.
An analytic arc can have local parameters beginning in different orders,
so a tangent-space argument by itself does not settle smoothability.
Paper 4 supplies the additional period and surgery argument in Sections 3–6,
and proves that $X^\circ$ is nonsmoothable. In particular, the formerly
separate order-of-contact problem is resolved for this example.

The distinction remains useful when interpreting the finite matrix:
Gross's hull theorem alone does not exclude an equation of the form
$\beta_1+\text{(quadratic terms)}=0$ with $\beta_1$ nonzero at higher order.
The all-orders conclusion uses the period theorem, not this linear-algebra
calculation alone.
