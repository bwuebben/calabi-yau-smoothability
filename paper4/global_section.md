# The fixed generic anticanonical section

**Status:** exact global-equation record, 20 August 2026.

This note fixes the anticanonical hypersurface used throughout Paper 4.  Its
companion checker is [`global_section.sage`](global_section.sage).  The
calculation supplies the defining-section input for a future
cotangent-complex computation; it does not compute the abstract global
deformation ring or any quadratic Kuranishi map.

## 1. The coefficient chart

The Newton polytope \(\Delta\) has exactly 25 lattice points.  On the dense
torus of \(\mathbb P_{\Delta^\vee}\), write

$$
f=\sum_{m\in\Delta\cap N}c_m x^m.
$$

Work on the open set where the coefficients at

$$
0,\quad e_1,\quad e_2,\quad e_3,\quad e_4
$$

are nonzero.  Overall scaling and the four-dimensional dense-torus action
normalize these five coefficients to one.  The exponent differences
\(e_i-0\) form a lattice basis, so this normalization is unique on that open
set.  The remaining coefficients give the following 20-parameter chart.

| Parameter | Exponent | Parameter | Exponent |
|---|---|---|---|
| \(a_{01}\) | \((-2,0,-1,0)\) | \(a_{11}\) | \((-1,1,-1,1)\) |
| \(a_{02}\) | \((-2,1,-1,0)\) | \(a_{12}\) | \((-1,1,0,0)\) |
| \(a_{03}\) | \((-2,1,-1,1)\) | \(a_{13}\) | \((0,-1,-1,0)\) |
| \(a_{04}\) | \((-1,-1,-1,0)\) | \(a_{14}\) | \((0,-1,0,-1)\) |
| \(a_{05}\) | \((-1,-1,0,-1)\) | \(a_{15}\) | \((0,-1,0,0)\) |
| \(a_{06}\) | \((-1,0,-1,0)\) | \(a_{16}\) | \((0,0,-1,0)\) |
| \(a_{07}\) | \((-1,0,-1,1)\) | \(a_{17}\) | \((0,0,-1,1)\) |
| \(a_{08}\) | \((-1,0,0,-1)\) | \(a_{18}\) | \((0,0,0,-1)\) |
| \(a_{09}\) | \((-1,0,0,0)\) | \(a_{19}\) | \((0,0,1,-1)\) |
| \(a_{10}\) | \((-1,1,-1,0)\) | \(a_{20}\) | \((1,-1,0,0)\) |

Set

$$
A=\mathbb Q[a_{01},\ldots,a_{20}],\qquad K=\operatorname{Frac}(A).
$$

The fixed generic Laurent polynomial is

$$
f_{\mathrm{gen}}
=1+x_1+x_2+x_3+x_4
 +\sum_{j=1}^{20}a_jx^{m_j}
\in K[x_1^{\pm1},\ldots,x_4^{\pm1}],
$$

where \(m_j\) is the exponent paired with \(a_j\) in the table.  Here and
in the checker, subscripts are padded to two digits only for sorting.

The dimension 20 agrees numerically with
\(h^{2,1}(\widehat X^\circ)=20\).  This agreement does **not** by itself
identify these coefficient directions with all locally trivial deformations
of \(X^\circ\), nor does it prove that the parameter \(r\) in the formal
support calculation equals 20.  Such an identification requires a comparison
of deformation functors.

## 2. Cox homogenization

Let \(q_1,\ldots,q_{26}\) be the ordered rays fixed in
[`example_data.md`](example_data.md), and let \(X_1,\ldots,X_{26}\) be the
corresponding Cox variables.  The Laurent monomial \(x^m\) homogenizes to

$$
M_m=\prod_{j=1}^{26}X_j^{\langle m,q_j\rangle+1}.
$$

Reflexivity gives nonnegative exponents.  The exact Cox grading from
[`cox_data.sage`](cox_data.sage) verifies

$$
\deg M_m=\deg(-K_{\mathbb P_{\Delta^\vee}})
$$

for all 25 lattice points.  Consequently

$$
F_{\mathrm{gen}}=\sum_{m\in\Delta\cap N}c_mM_m
$$

is a homogeneous anticanonical equation over \(K\).  This is the equation to
be restricted to affine toric charts in the intrinsic computation.

## 3. Generic face nondegeneracy

The face lattice of \(\Delta\) has the exact \(f\)-vector

$$
(f_0,f_1,f_2,f_3,f_4)=(22,68,72,26,1),
$$

so there are 189 nonempty faces.  For a face \(\Theta\), let \(f_\Theta\)
be the corresponding face polynomial.  Its dense-torus nondegeneracy is the
condition

$$
(f_\Theta,
x_1\partial_{x_1}f_\Theta,\ldots,
x_4\partial_{x_4}f_\Theta)=(1)
$$

in the Laurent polynomial ring.  This is a Zariski-open condition on the
coefficient chart.

The checker proves that every one of these 189 opens is nonempty.  It uses
the coefficient vector

$$
(2,3,5,7,11,13,17,19,23,29,
31,37,41,43,47,53,59,61,67,71)
$$

for every proper face and the all-one coefficient vector for the
full-dimensional face; all ideal-membership checks are exact over
\(\mathbb Q\).  These facewise witnesses already suffice: since the
coefficient chart is irreducible and the family of opens is finite, their
intersection is nonempty and contains its generic point.  Thus
\(f_{\mathrm{gen}}\) is \(\Delta\)-regular.

The opt-in command

```bash
sage paper4/global_section.sage --verify-common-witness
```

proves the stronger specialization statement needed by the projective
cocycle calculation.  It checks exactly that the same distinct-prime vector
is also nondegenerate on the full-dimensional face.  Because that vector is
already used on every proper face above, it is one common
\(\Delta\)-regular rational witness for all 189 faces.  The full-face
Gröbner calculation is intentionally opt-in because it takes substantially
longer than the standard facewise proof.

## 4. The 30 marked orbit points

Every singular two-face \(F\subset\Delta^\vee\) has a dual edge
\(F^*=[m_-,m_+]\subset\Delta\) with exactly two lattice points.  With the
orbit coordinate

$$
z=x^{m_+-m_-},
$$

the face restriction is a unit times \(c_-+c_+z\).  Hence the unique point
of the hypersurface on that orbit is

$$
z_F=-\frac{c_-}{c_+}.
$$

The following table fixes the orientation and the resulting root used by
later local expansions.

| Germ | Oriented dual edge \(m_-\to m_+\) | \(z_F\) |
|---|---|---|
| \(N_1\) | \((0,0,1,0)\to(1,0,0,0)\) | \(-1\) |
| \(N_2\) | \((0,0,1,0)\to(0,1,0,0)\) | \(-1\) |
| \(N_3\) | \((0,0,0,1)\to(1,0,0,0)\) | \(-1\) |
| \(N_4\) | \((0,0,1,-1)\to(1,0,0,0)\) | \(-a_{19}\) |
| \(N_5\) | \((0,0,1,-1)\to(0,1,0,0)\) | \(-a_{19}\) |
| \(N_6\) | \((0,0,-1,1)\to(1,0,0,0)\) | \(-a_{17}\) |
| \(N_7\) | \((-1,1,-1,1)\to(0,1,0,0)\) | \(-a_{11}\) |
| \(N_8\) | \((-1,1,-1,1)\to(0,0,-1,1)\) | \(-a_{11}/a_{17}\) |
| \(N_9\) | \((0,0,-1,1)\to(0,0,0,1)\) | \(-a_{17}\) |
| \(N_{10}\) | \((0,0,0,-1)\to(1,0,0,0)\) | \(-a_{18}\) |
| \(N_{11}\) | \((0,0,0,-1)\to(0,0,1,-1)\) | \(-a_{18}/a_{19}\) |
| \(N_{12}\) | \((0,0,-1,0)\to(1,0,0,0)\) | \(-a_{16}\) |
| \(N_{13}\) | \((-1,1,-1,0)\to(0,1,0,0)\) | \(-a_{10}\) |
| \(N_{14}\) | \((-1,1,-1,0)\to(0,0,-1,0)\) | \(-a_{10}/a_{16}\) |
| \(N_{15}\) | \((0,0,1,0)\to(1,-1,0,0)\) | \(-1/a_{20}\) |
| \(N_{16}\) | \((0,0,0,1)\to(1,-1,0,0)\) | \(-1/a_{20}\) |
| \(N_{17}\) | \((0,-1,0,-1)\to(1,-1,0,0)\) | \(-a_{14}/a_{20}\) |
| \(N_{18}\) | \((0,-1,0,-1)\to(0,0,0,-1)\) | \(-a_{14}/a_{18}\) |
| \(N_{19}\) | \((0,-1,-1,0)\to(1,-1,0,0)\) | \(-a_{13}/a_{20}\) |
| \(N_{20}\) | \((0,-1,-1,0)\to(0,0,-1,0)\) | \(-a_{13}/a_{16}\) |
| \(N_{21}\) | \((-1,1,0,0)\to(0,1,0,0)\) | \(-a_{12}\) |
| \(N_{22}\) | \((-2,1,-1,1)\to(-1,1,-1,1)\) | \(-a_{03}/a_{11}\) |
| \(N_{23}\) | \((-2,1,-1,0)\to(-1,1,-1,0)\) | \(-a_{02}/a_{10}\) |
| \(N_{24}\) | \((0,-1,0,0)\to(1,-1,0,0)\) | \(-a_{15}/a_{20}\) |
| \(N_{25}\) | \((-1,-1,0,-1)\to(0,-1,0,-1)\) | \(-a_{05}/a_{14}\) |
| \(N_{26}\) | \((-1,-1,-1,0)\to(0,-1,-1,0)\) | \(-a_{04}/a_{13}\) |
| \(D_{6,1}\) | \((0,1,0,0)\to(1,0,0,0)\) | \(-1\) |
| \(D_{6,2}\) | \((1,-1,0,0)\to(1,0,0,0)\) | \(-a_{20}\) |
| \(D_{7,1}\) | \((0,0,0,1)\to(0,0,1,0)\) | \(-1\) |
| \(D_{7,2}\) | \((0,0,1,-1)\to(0,0,1,0)\) | \(-a_{19}\) |

The 30 dual edges are distinct.  Each root is nonzero over \(K\), and the
edge direction is primitive.  This both locates every marked germ in the
generic fibre and supplies the invertible orbit derivative used in the graph
argument of [`global_map.md`](global_map.md), Lemma 2.1.

## 5. What this settles and what it does not

The computation now fixes, without an unnamed genericity assumption:

1. the 25 monomials and the 20 normalized coefficients;
2. a homogeneous Cox equation in the printed grading;
3. generic nondegeneracy on every face; and
4. the location of all 30 marked singularities on their one-dimensional
   torus orbits.

It therefore removes the defining equation as an ambiguity in Route A of
[`global_map.md`](global_map.md).  The next missing object is an affine-cover
model of the coarse toric ambient space and of this Cartier equation, with
enough syzygies and overlap maps to compute
\(\operatorname{Ext}^1(L_{X^\circ},\mathcal O_{X^\circ})\), its localization,
and the five Picard-paired quadratic classes.  Nothing in the present
calculation supplies those brackets or decides the formal smoothing question.
