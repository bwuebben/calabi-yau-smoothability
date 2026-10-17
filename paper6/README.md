# Mirror Symmetry for Conifold Relations in Calabi–Yau Threefolds

Paper 6 · Bernd Johannes Wuebben · 214 pages

[PDF](cy-vanishing-spheres.pdf) · [LaTeX source](main.tex) · [Series overview](../README.md)

For a class of Calabi–Yau hypersurfaces in toric fourfolds, and for suitable members of the Batyrev
mirror family near the maximally degenerate limit, we prove that the exceptional curves of a small
resolution and the mirror vanishing three-spheres have the same integral relation lattice, and that
the spheres generate a torsion-free subgroup. Every relation is realised by an integral four-chain
on the mirror hypersurface, constructed from a tropical two-cycle, relative to the discriminant, on
the base of Gross's toric degeneration. Near these members the locus where all the spheres collapse
is smooth of codimension equal to the relative Picard number; this establishes, locally in the
coefficient space of the mirror family, Morrison's prediction that the mirror nodes impose exactly
that number of conditions. We characterise smoothability by relations among the mirror spheres and
obtain a symplectic conifold transition on the mirror whenever a smoothing exists. The results apply
to nodal crepant partial resolutions, induced by projective toric morphisms, of anticanonical
hypersurfaces whose singular points are nodes and cones over del Pezzo surfaces of degrees six and
seven, one on each singular two-face of the polytope. The unimodular triangulations needed for
Gross's construction are verified to exist by computer for the 51,827 such polytopes identified in
the Kreuzer–Skarke list.

Theorem E supplies a supplementary geometric realization on the small resolution:
a continuous map from a completed circle domain, with signed degrees on its boundary
spheres. Section 10 explains the construction and X₉; Appendix F contains the full proof.
Theorem F verifies unimodular heights for the fixed list of 51,827 polytopes. The separate
classification-completeness proposition retains the claim that this list contains every
admissible class, with the coverage of the independent checks stated explicitly.

The source consists of `main.tex` and the section files it inputs.

## Build

From this directory:

```sh
latexmk -pdf main.tex
```

## Computations

The [computation bundle](COMPUTATIONS.md) contains:

- [Unimodular-height certificates](unimodular/) for the fixed list of 51,827 admissible
  polytopes, with integer and exact rational checking programs; completeness of the
  list is a separate assertion supported by the classification and admissibility filter.
- [Exact data and checks for the five worked examples](examples/).
- [Forced-node and projectivity data](forced_nodes/) for the finite input list, including
  projectivity witnesses for every profile of the examples on the branch that admits one.
- [Height-condition results](census/height_condition.json) for the 364 hexagon-free
  polytopes of the list, with a [checker](census/check_height_condition.py).
- [Tropical-cycle census summaries](census/README.md), with their verification
  scope and retained-evidence limitations stated explicitly.
- [Classification input provenance](provenance/DATASET_CHECK.md) and a
  [file-integrity manifest](MANIFEST.json).

See [COMPUTATIONS.md](COMPUTATIONS.md) for dependencies and reproduction commands.
The large raw classification files are external; their pinned hashes are included.
