# Link periods and deformations of Calabi–Yau threefolds with toric singularities

Bernd Johannes Wuebben · 26 November 2026 · Preprint

[Paper](main.pdf) · [Paper page](https://bwuebben.github.io/papers/toric-link-periods/) · [Zenodo](https://doi.org/10.5281/zenodo.24575797) · [Source and computations](submission-source.tar.gz)

The periods of a relative volume form on the link of an isolated Gorenstein toric threefold singularity embed its semiuniversal deformation germ in the third cohomology of the link. For the homogeneous toric form this embedding is linear in Altmann’s parameters. Multiplication by a holomorphic unit changes the embedding by an automorphism of its source, even over nonreduced bases. The automorphism comes from a common change of the toric height coordinate at all edges of the defining polygon. Let X be a connected normal projective threefold with trivial dualizing sheaf, vanishing irregularity and only such singularities. In the product of the local period images, impose the linear condition that the Poincaré dual link classes sum to zero in the second homology of a projective log resolution. Using a local-to-global lifting theorem for small Artin extensions, we identify the semiuniversal germ of X with the product of this scheme-theoretic linear section and a smooth germ of dimension $h^1(X,\Theta_X)$. We construct a projective example in which a first-order deformation extends to second order and is obstructed at third order.

## Reproduction

Build from this directory with the following commands. The TeX build uses TeX Live 2026. The calculations use exact rational polynomial arithmetic.

~~~sh
env SOURCE_DATE_EPOCH=1795694400 FORCE_SOURCE_DATE=1 latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
uv run --no-project --with sympy==1.14.0 python computations/verify_toric_periods.py
uv run --no-project --with sympy==1.14.0 python computations/octagon.py
uv run --no-project --with sympy==1.14.0 python computations/compact_cubic.py
~~~

© 2026 Bernd Johannes Wuebben. All rights reserved.
