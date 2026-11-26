# Local equations for deformations of Calabi–Yau threefolds with del Pezzo cone points

Bernd Johannes Wuebben · 26 November 2026 · Preprint

[Read the paper](main.pdf) · [Paper page](https://bwuebben.github.io/papers/stratified-cy-deformation-germs/) · [Cite the preprint](https://doi.org/10.5281/zenodo.24573015) · [Source and calculations](submission-source.tar.gz)

We determine the analytic deformation space of a projective Calabi–Yau threefold with ordinary double points and anticanonical cones over smooth del Pezzo surfaces of degrees five and six, admitting a smooth projective crepant resolution small at the nodes and given by the vertex blow-up at the cones. Up to a smooth factor, its semiuniversal base is the scheme-theoretic intersection of the product of the marked local bases with the space of homological relations on the resolution. Link periods identify the strata on which prescribed singularities persist. We construct two projective examples with the same smooth reduced germ and the same tangent dimension. The nilradical has square zero in one and nonzero square in the other. The proof establishes compatibility of link periods with Artin base change and detects global obstructions to lifting by localization and pairing with Cartier divisor classes.

## Reproducing the calculations

From this directory build the manuscript with `latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex`. The same sources and computations are included in the archive. With Python 3.10 or later and SageMath 10.9, run

```sh
python3 computations/check_models.py | diff -u computations/verified.txt -
sage -python computations/search_subsets.py | diff -u computations/search_results.txt -
```

The first computation starts with the lattice vertices and marked hexagons. The second uses the four required inventory entries and curve matrices extracted from the census accompanying Paper V. The two computations reproduce the compact examples and the bounded search, respectively.

© 2026 Bernd Johannes Wuebben. All rights reserved.
