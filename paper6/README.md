# Mirror Symmetry for Conifold Relations in Calabi–Yau Threefolds

Paper 6 · Bernd Johannes Wuebben · 234 pages

[PDF](cy-vanishing-spheres.pdf) · [LaTeX source](main.tex) · [Series overview](../README.md)

For a class of Batyrev mirror families, we establish Morrison's predicted codimension for
nodal degenerations and strengthen it to an equality of integral relation lattices. The
result holds on every connected component of the nodal coefficient locus meeting our
construction near the maximally degenerate limit. After marking and orienting the nodes, the
relations among exceptional curves of the original small resolution coincide with those
among mirror vanishing three-spheres, which generate a torsion-free subgroup. Each such
component is smooth of codimension equal to the relative Picard number. We realise every
relation near the limit by an integral four-chain lifted from a tropical two-cycle relative
to the discriminant on Gross's degeneration base. The comparison characterises smoothability
of the nodal partial resolution and, when it is smoothable, yields a symplectic conifold
transition on the mirror. We treat specified projective toric partial resolutions of
hypersurfaces with nodes and degree-six or degree-seven del Pezzo cone points. In an
explicit four-node family, every member of the unique codimension-three nodal component
admits a projective small resolution. For general members, these resolutions are birational
to Batyrev mirrors of smoothings of the original nodal threefold, giving the reversed
transition predicted by Morrison in this family.

Theorem E supplies a supplementary geometric realization on the small resolution:
a continuous map from a completed circle domain, with signed degrees on its boundary
spheres. Section 10 explains the construction and X₉; Appendix F contains the full proof.
Theorem F verifies unimodular heights for the fixed list of 51,827 polytopes. The separate
classification-completeness proposition retains the claim that this list contains every
admissible class, with the coverage of the independent checks stated explicitly.

The componentwise extension is Proposition 9.13 (§9.7); the exact-node
locus is smooth everywhere by Lemma 9.17. Theorem 9.21 (§9.9) classifies the
entire codimension-three four-node locus of X₉, including members with zero
coefficients. It is one reached determinantal component, and every member
has a projective small resolution. Other four-node components have
codimension four and no relation. Theorem 9.25 (§9.10) identifies the general
resolved mirror with the Batyrev partner family of a smoothing of the exact
original nodal partial resolution, up to birational equivalence. The
original smoothing is proper and analytic with projective fibres; quantum
compatibility and the geometric reverse transition for the entire
admissible class remain open.

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
- [Exact data and checks for the five worked examples](examples/), including
  the full X₉ four-node component classification and reverse-transition inputs.
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
