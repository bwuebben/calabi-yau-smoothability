# Mirror Symmetry for Conifold Relations in Calabi–Yau Threefolds

Paper 6 · Bernd Johannes Wuebben · 190 pages

[PDF](cy-vanishing-spheres.pdf) · [LaTeX source](main.tex) · [Series overview](../README.md)

We establish and strengthen Morrison's numerical prediction for conifold transitions
under mirror symmetry for a class of toric Calabi–Yau threefolds. We consider
anticanonical hypersurfaces with nodes and isolated cones over del Pezzo surfaces of
degrees six and seven, at most one singularity on each torus orbit, and a nodal crepant
partial resolution induced by a projective toric morphism. For suitable degenerations of
the Batyrev mirror family near the maximally degenerate limit, we prove that the
exceptional curves of a small resolution and the mirror vanishing spheres have the same
integral relation lattice, and that the spheres generate a torsion-free subgroup. Every
relation is realised by an explicit integral four-chain lifted from a relative tropical
two-cycle on the base of Gross's toric degeneration. The corresponding nodal locus is
smooth of codimension equal to the relative Picard number, as Morrison predicts.
Consequently, the nodal partial resolution is smoothable precisely when the mirror
spheres admit a relation with every coefficient nonzero. The required triangulations are
established by a computer-assisted verification over the specified class of reflexive
polytopes.

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
