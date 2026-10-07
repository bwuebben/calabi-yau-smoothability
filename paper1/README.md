# Non-smoothable Calabi–Yau threefolds from reflexive polytopes

Paper 1 · Bernd Johannes Wuebben · 18 pages · [arXiv:2609.28499](https://arxiv.org/abs/2609.28499)

[PDF](cy-non-smoothable.pdf) · [LaTeX source](main.tex) · [Series overview](../README.md)

Explicit projective Calabi–Yau threefolds with one nonsmoothable singular
point: an anticanonical cone over F₁, or a pentagon cone with a
positive-dimensional reduced deformation base and no smoothing component.
The local classification contains 217 isolated Gorenstein toric singularity
classes with primitive edge vectors in [−2,2]² in some lattice basis.

The scan of the verified database copy finds that
39,175,536 polytopes carry the sufficient local obstruction, approximately
8.27% of the Kreuzer–Skarke classification. This gives a lower bound for
global nonsmoothability. Type (R) means reduced-rigid, allowing nonzero
first-order deformations; the smooth unimodular triangle is excluded.
Gross's earlier constructions and the explicit toric examples are
identified separately.

## Build

From this directory:

```sh
latexmk -pdf main.tex
```

The resulting `main.pdf` is the local build output; `cy-non-smoothable.pdf` is the
published PDF.

## Computations

Run from the repository root:

```sh
python3 src/toric_census.py
python3 src/batyrev_global.py
python3 src/hodge_numbers.py
```

The [series reproduction guide](../README.md#reproducing) lists the complete
commands and distinguishes standard Python from SageMath. Database-wide
reproduction and the input manifest are described in [DATASET_CHECK.md](../DATASET_CHECK.md).
