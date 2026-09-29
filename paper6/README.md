# Mirror Symmetry for Conifold Relations in Calabi–Yau Threefolds

Paper 6 · Bernd Johannes Wuebben · 206 pages

[PDF](cy-vanishing-spheres.pdf) · [LaTeX source](main.tex) · [Series overview](../README.md)

For a class of Calabi–Yau hypersurfaces in toric fourfolds, we establish Morrison's
numerical prediction for conifold transitions under mirror symmetry, as a count in the
coefficient space of the Batyrev mirror family, and strengthen it to an integral
statement. For suitable degenerations of the Batyrev mirror family near the maximally
degenerate limit, we prove that the exceptional curves of a small resolution and the
mirror vanishing three-spheres have the same integral relation lattice, and that the
spheres generate a torsion-free subgroup. Every relation is realised by an explicit
integral four-chain on the mirror hypersurface, constructed from a relative tropical
two-cycle. Near the mirror members constructed, the locus where all the spheres collapse
is smooth of codimension equal to the relative Picard number. We characterise
smoothability by relations among the mirror spheres and obtain a symplectic conifold
transition on the mirror whenever a smoothing exists. The results apply to nodal crepant
partial resolutions, induced by projective toric morphisms, of a specified class of
anticanonical hypersurfaces with nodes and cones over del Pezzo surfaces of degrees six
and seven. The required triangulations are established by computer-assisted
verification.

The source consists of `main.tex` and the section files it inputs.

## Build

From this directory:

```sh
latexmk -pdf main.tex
```

## Computations

The [computation bundle](COMPUTATIONS.md) contains:

- [Unimodular-height certificates](unimodular/) for all 51,827 admissible
  polytopes, with integer and exact rational checking programs.
- [Exact data and checks for the five worked examples](examples/).
- [Forced-node and projectivity data](forced_nodes/) for the finite input list.
- [Tropical-cycle census summaries](census/README.md), with their verification
  scope and retained-evidence limitations stated explicitly.
- [Classification input provenance](provenance/DATASET_CHECK.md) and a
  [file-integrity manifest](MANIFEST.json).

See [COMPUTATIONS.md](COMPUTATIONS.md) for dependencies and reproduction commands.
The large raw classification files are external; their pinned hashes are included.
