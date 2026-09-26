# Tropical 2-cycles and relations among vanishing spheres in Batyrev mirror families

Paper 6 · Bernd Johannes Wuebben · 190 pages

[PDF](cy-vanishing-spheres.pdf) · [LaTeX source](main.tex) · [Series overview](../README.md)

Morrison proposed that mirror symmetry exchanges the two sides of a conifold
transition. This paper proves that prediction for Batyrev hypersurfaces whose
singular points are nodes and cones over del Pezzo surfaces of degrees six and
seven, under the hypothesis that the crepant partial resolution with only nodes
is induced by a projective toric morphism. Near the maximally degenerate limit,
the vanishing spheres of the mirror family satisfy exactly the integral
relations of the exceptional curves of a small resolution, and the members with
all these nodes form a smooth locus whose codimension is the relative Picard
number. Consequently the threefold is smoothable exactly when the vanishing
spheres of its mirror satisfy a relation with every coefficient nonzero.

Every relation is realised by an explicit integral chain lifted from a relative
tropical 2-cycle on the base of Gross's toric degeneration; that the node
coefficients of these cycles generate the relation lattice is proved by a
homology computation in the complement of the discriminant. For the tropical
2-cycles of Castaño-Bernard and Matessi the paper also constructs the
corresponding chains on the small resolution. A combinatorial condition used in
the construction (a unimodular triangulation with a convex height) is verified
for all admissible reflexive polytopes by a computer-assisted check of the
Kreuzer–Skarke list.

The source consists of `main.tex` and the section files it inputs.

## Build

From this directory:

```sh
latexmk -pdf main.tex
```

## Computations

The programs and data cited in the paper (the unimodular-height certificates
and their two checking programs, the data of the examples, and the counts of
Section 11) will be added to this directory.
