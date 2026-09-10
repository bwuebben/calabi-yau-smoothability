# Doubly isolated Batyrev mirror pairs and non-smoothable Calabi–Yau threefolds

Paper 3 · Bernd Johannes Wuebben · 11 pages

[PDF](cy-mirror-pairs.pdf) · [LaTeX source](main.tex) · [Series overview](../README.md)

The scan of the verified database copy classifies 590
reflexive four-polytopes whose two associated Batyrev hypersurfaces have
at most isolated singularities. Their two-faces form twelve lattice classes.
The F₁ cone occurs in a unique mirror pair in this classification; its
explicit 22/26-vertex polytopes and resolution Hodge numbers (20,26)/(26,20)
are verified independently of database completeness. The self-dual 24-cell
is the unique smooth/smooth pair within the classification.

The local nonsmoothable types are two cyclic quotients and the F₁ cone.
Reduced rigidity is distinguished from absence of first-order deformations.
The mirror whose germs are all locally smoothable is globally nonsmoothable
by the mixed obstruction in [Paper 4](../paper4/README.md).

## Build

From this directory:

```sh
latexmk -pdf main.tex
```

The resulting `main.pdf` is the local build output; `cy-mirror-pairs.pdf` is the
published PDF.

## Computations

Run from the repository root:

```sh
python3 src/both_sides_census.py
python3 src/face_data.py
python3 src/paper3_node_relations.py
```

The [series reproduction guide](../README.md#reproducing) lists the complete
commands and distinguishes standard Python from SageMath. Database-wide
reproduction and the input manifest are described in [DATASET_CHECK.md](../DATASET_CHECK.md).
