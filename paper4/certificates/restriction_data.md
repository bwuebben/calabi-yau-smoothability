# Divisor restrictions on the fixed resolution

**Status:** exact resolution-level computation, 20 August 2026. The matrices
in this note are derived and checked by
[`restriction_data.sage`](restriction_data.sage) from the rays in
[`example_data.md`](example_data.md) and the diagonals in
[`resolution_data.md`](resolution_data.md). They are numerical geometry on a
fixed crepant resolution. Section 5 proves that the ambient Picard space is
the full Picard space of the resolved Calabi--Yau over $\mathbb Q$ and
$\mathbb C$. The deformation-theoretic comparison is treated in Section 2.4 and
Lemma 4.2 of the paper (`main.tex`).

## 1. Ambient Picard basis

Let $Y_{\widehat\Sigma}$ be the smooth toric fourfold defined by the fixed
30-ray fan. The ray matrix gives an exact sequence

$$
0\longrightarrow M\longrightarrow\mathbb Z^{30}
\longrightarrow\operatorname{Pic}(Y_{\widehat\Sigma})\longrightarrow0.
$$

Its Smith diagonal is $(1,1,1,1)$, so the Picard group is torsion-free of
rank $26$. The four rays

$$
(q_{15},q_{22},q_{24},q_{26})
$$

form a unimodular maximal cone. The classes of the remaining divisors form
the fixed basis

$$
\mathcal B=
(D_1,D_2,\ldots,D_{14},D_{16},\ldots,D_{21},D_{23},D_{25},
D_{27},D_{28},D_{29},D_{30}).
$$

In this basis the anticanonical class is

$$
\begin{aligned}
[-K_Y]_{\mathcal B}={}&(3,2,3,1,2,1,3,2,3,1,2,1,2,1,\\
&1,2,1,2,1,1,2,1,2,2,2,1).
\end{aligned}
$$

The checker constructs the complete $26\times30$ grading matrix and verifies
that it annihilates the ray matrix. The identity columns indexed by
$\mathcal B$ and the four derived pivot columns determine it uniquely.

## 2. The exceptional surfaces

The restriction of the fixed fan to each divisorial polygon is the full star
at its interior ray. In the quotient by that ray, the counterclockwise fan
rays and boundary self-intersections are:

| surface | cyclic boundary rays | boundary self-intersections | $K_E^2$ |
|---|---|---|---:|
| $E_{6,1}$ | $(q_1,q_3,q_5,q_6,q_4,q_2)$ | $(-1,-1,-1,-1,-1,-1)$ | 6 |
| $E_{6,2}$ | $(q_7,q_9,q_{11},q_{12},q_{10},q_8)$ | $(-1,-1,-1,-1,-1,-1)$ | 6 |
| $E_{7,1}$ | $(q_1,q_{13},q_{21},q_{17},q_7)$ | $(-1,0,0,-1,-1)$ | 7 |
| $E_{7,2}$ | $(q_2,q_{14},q_{22},q_{18},q_8)$ | $(-1,-1,0,0,-1)$ | 7 |

Thus the surfaces are the intended toric $dP_6,dP_6,dP_7,dP_7$. All
adjacencies, smooth-cone determinants, self-intersections, and degrees in
this table are recomputed from the ray coordinates.

### 2.1 The two $dP_6$ surfaces

For $E_{6,1}$ use the boundary-divisor basis

$$
([D_3|_E],[D_4|_E],[D_5|_E],[D_6|_E]);
$$

for $E_{6,2}$ use the translated basis with indices $(9,10,11,12)$. With
the boundary columns in increasing ray order, $(1,2,3,4,5,6)$ or
$(7,8,9,10,11,12)$, the common grading matrix is

$$
Q_6=
\begin{pmatrix}
-1&1&1&0&0&0\\
1&-1&0&1&0&0\\
0&1&0&0&1&0\\
1&0&0&0&0&1
\end{pmatrix}.
$$

The canonical class and intersection matrix in the chosen basis are

$$
[K_{E_6}]=(-1,-1,-2,-2),
\qquad
I_6=
\begin{pmatrix}
-1&0&1&0\\
0&-1&0&1\\
1&0&-1&1\\
0&1&1&-1
\end{pmatrix},
$$

and $[K_{E_6}]I_6[K_{E_6}]^{\mathsf T}=6$.

### 2.2 The first $dP_7$ surface

For $E_{7,1}$ use the basis with boundary rays $(q_{13},q_{17},q_{21})$.
With columns $(q_1,q_7,q_{13},q_{17},q_{21})$,

$$
Q_{7,1}=
\begin{pmatrix}
-1&1&1&0&0\\
1&-1&0&1&0\\
1&0&0&0&1
\end{pmatrix},
\quad
[K_{E_{7,1}}]=(-1,-1,-2),
$$

and

$$
I_{7,1}=
\begin{pmatrix}
0&0&1\\
0&-1&1\\
1&1&0
\end{pmatrix}.
$$

### 2.3 The second $dP_7$ surface

For $E_{7,2}$ use the basis with boundary rays $(q_{14},q_{18},q_{22})$.
With columns $(q_2,q_8,q_{14},q_{18},q_{22})$,

$$
Q_{7,2}=
\begin{pmatrix}
-1&1&1&0&0\\
1&-1&0&1&0\\
0&1&0&0&1
\end{pmatrix},
\quad
[K_{E_{7,2}}]=(-1,-1,-2),
$$

and

$$
I_{7,2}=
\begin{pmatrix}
-1&0&1\\
0&0&1\\
1&1&0
\end{pmatrix}.
$$

Both displayed $dP_7$ canonical vectors have square $7$ in their respective
intersection forms.

## 3. Restriction from the ambient Picard lattice

At the generic point of the singular torus orbit, a boundary divisor of the
face restricts to the corresponding toric boundary curve of $E$. The
exceptional divisor restricts to

$$
D_E|_E=K_E
$$

by adjunction on the crepant Calabi--Yau resolution. All other invariant
divisors miss this generic exceptional fiber. Consequently the complete map
from invariant divisors is obtained from $Q_6,Q_{7,1},Q_{7,2}$ by adding the
canonical column at $q_{27},q_{28},q_{29},q_{30}$, respectively, and zero
columns elsewhere.

Each resulting matrix annihilates all four principal-divisor columns and
therefore factors integrally through $\operatorname{Pic}(Y_{\widehat\Sigma})$.
The nonzero columns in the ambient basis $\mathcal B$ are:

| surface | nonzero basis columns |
|---|---|
| $E_{6,1}$ | $D_1=(-1,1,0,1)$, $D_2=(1,-1,1,0)$, $D_3=e_1$, $D_4=e_2$, $D_5=e_3$, $D_6=e_4$, $D_{27}=K_{E_6}$ |
| $E_{6,2}$ | $D_7=(-1,1,0,1)$, $D_8=(1,-1,1,0)$, $D_9=e_1$, $D_{10}=e_2$, $D_{11}=e_3$, $D_{12}=e_4$, $D_{28}=K_{E_6}$ |
| $E_{7,1}$ | $D_1=(-1,1,1)$, $D_7=(1,-1,0)$, $D_{13}=e_1$, $D_{17}=e_2$, $D_{21}=e_3$, $D_{29}=K_{E_{7,1}}$ |
| $E_{7,2}$ | $D_2=(-1,1,0)$, $D_8=(1,-1,1)$, $D_{14}=e_1$, $D_{18}=e_2$, $D_{30}=K_{E_{7,2}}$ |

Here the $e_j$ in a row are the standard basis vectors of that surface's
Picard group; they are unrelated to the exceptional ray notation used in
`resolution_data.md`. The class of $D_{22}|_{E_{7,2}}$ is already encoded by
the four pivot-degree columns of the ambient grading.

The stacked restriction map

$$
R:\operatorname{Pic}(Y_{\widehat\Sigma})
\longrightarrow
\operatorname{Pic}(E_{6,1})\oplus\operatorname{Pic}(E_{6,2})
\oplus\operatorname{Pic}(E_{7,1})\oplus\operatorname{Pic}(E_{7,2})
$$

is a $14\times26$ integral matrix of rank $14$. Its 14 nonzero Smith
invariants are all $1$, so it is surjective as a map of lattices.

## 4. The node curves and the combined numerical map

For each square, orient its primitive circuit relation so that the two rays
of the chosen diagonal have coefficient $-1$. This is the toric divisor
intersection row of the corresponding exceptional curve $C_i$. The resulting
$26\times26$ matrix on $\mathcal B$ has rank $18$, agreeing with the
resolution-independent rank of the node-relation subsystem.

Stacking the curve rows with the surface restrictions gives

$$
\begin{aligned}
\Phi_{\mathrm{num}}:\operatorname{Pic}(Y_{\widehat\Sigma})&\longrightarrow
\mathbb Z^{26}\oplus\bigoplus_E\operatorname{Pic}(E),\\
D&\longmapsto
\bigl((D\cdot C_i)_{i=1}^{26},(D|_E)_E\bigr).
\end{aligned}
$$

Exact integer elimination gives

$$
\operatorname{rank}\Phi_{\mathrm{num}}=25,
\qquad
\ker\Phi_{\mathrm{num}}=\mathbb Z[-K_Y].
$$

All 25 nonzero Smith invariants are $1$. Thus the image is saturated and the
kernel generator is primitive. Geometrically, the exceptional curves and
surfaces detect every ambient Picard direction except the anticanonical
polarization direction, which is numerically trivial on all crepant
exceptional fibers.

## 5. The ambient Picard space is the Calabi--Yau Picard space

Let $\widehat X\subset Y_{\widehat\Sigma}$ be the smooth anticanonical
proper transform. The restriction map

$$
\operatorname{Pic}(Y_{\widehat\Sigma})
\longrightarrow \operatorname{Pic}(\widehat X)
$$

is injective. Indeed, if an ambient class $D$ restricts trivially to
$\widehat X$, it restricts trivially to each exceptional surface and has
degree zero on each exceptional node curve. The exact kernel calculation in
Section 4 then gives $D=m[-K_Y]$. On the other hand,

$$
\bigl((-K_Y)|_{\widehat X}\bigr)^3=(-K_Y)^4
=4!\operatorname{Vol}(\Delta_{\mathbb F_1})=94.
$$

The final equality is recomputed from the 22 printed vertices. Hence
$(-K_Y)|_{\widehat X}$ is not torsion, and $m=0$.

Batyrev's formula gives $h^{1,1}(\widehat X)=26$, while Section 1 gives
$\operatorname{rank}\operatorname{Pic}(Y_{\widehat\Sigma})=26$. Since a
smooth Calabi--Yau threefold has $H^1(\mathcal O)=H^2(\mathcal O)=0$, the
restriction map is therefore an isomorphism after tensoring with $\mathbb Q$
or $\mathbb C$:

$$
\operatorname{Pic}(Y_{\widehat\Sigma})\otimes\mathbb C
\simeq \operatorname{Pic}(\widehat X)\otimes\mathbb C.
$$

An integral index could remain; it is irrelevant to the complex
first-order obstruction map. In particular, the 26 ambient divisor
functionals detect every class in $H_2(\widehat X,\mathbb C)$.

## 6. What this does and does not prove

This calculation supplies the full complex Picard dual in which a mixed
analogue of Friedman's relation lives. It also shows that the four surface
terms add seven independent ranks beyond the rank-18 node subsystem, so the
surface data are not contained in the node detector. This rank statement
alone says nothing about whether the two node rows cease to be coloops in the
smaller $K_E^\perp$ subspaces relevant to local deformation parameters.

It does **not** by itself identify Altmann's chosen deformation coordinates
with the canonical link classes in the local-to-global sequence. The
topological comparison and the branch-preserving normalization are now
Section 2.4 and Lemma 4.2 of the paper (`main.tex`), whose vanishing-period
route settles the verdict at all orders of contact.
