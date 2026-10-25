# Mirror Symmetry for Conifold Relations in Calabi–Yau Threefolds

Paper 6 · Bernd Johannes Wuebben · 234 pages

[PDF](cy-vanishing-spheres.pdf) · [LaTeX source](main.tex) · [Series overview](../README.md)

For a class of Calabi–Yau hypersurfaces in toric fourfolds, we identify the integral relations
among exceptional curves of a small resolution with those among vanishing three-spheres in the
Batyrev mirror family, and exclude nonzero torsion combinations of the spheres. Near the
maximally degenerate limit, every relation is realised by an integral four-chain lifted from a
tropical two-cycle, relative to the discriminant, on Gross's degeneration base. The locus of
mirror hypersurfaces with the prescribed number of nodes and no other singularities is smooth,
with local equations given by independent vanishing-cycle periods. After marking and orienting
the nodes, the integral comparison extends along every connected component of this locus
meeting the near-limit construction, and the bounding chains continue by transport. On each
such component the codimension equals the relative Picard number, establishing Morrison's
predicted count there. In an explicit four-node family, we classify the entire codimension-
three locus as one such component; all other components have codimension four and no relation.
Every member of the former component has a projective small resolution, and for general members
these resolutions are birational to Batyrev mirrors of smoothings of the original nodal
threefold. The comparison characterises smoothability and yields a symplectic conifold
transition on the mirror whenever a smoothing exists. The results apply to nodal crepant
partial resolutions, induced by projective toric morphisms, of anticanonical hypersurfaces
whose singular points are nodes and cones over del Pezzo surfaces of degrees six and seven, one
on each singular two-face. The required unimodular triangulations are verified by computer for
the 51,827 such polytopes identified in the Kreuzer–Skarke list.

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
