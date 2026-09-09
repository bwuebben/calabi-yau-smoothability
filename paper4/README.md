# A vanishing-cycle obstruction to smoothing Calabi–Yau threefolds

Paper 4 · Bernd Johannes Wuebben · 32 pages

[PDF](cy-vanishing-cycle.pdf) · [LaTeX source](main.tex) · [Series overview](../README.md)

A homological necessary condition for smoothings with nodes and exact
anticanonical cones over smooth del Pezzo surfaces of degrees six and seven.
Vanishing periods force nonzero coefficients in the appropriate relation
space for analytic arcs of arbitrary order. The 26-node/four-cone mirror
example has no smoothing, although all its germs are smoothable.

The crepant resolution is small at nodes and has the specified exceptional
surface over each cone point. The equality between the matrix kernel and
the topological link-homology kernel uses these isolated-germ hypotheses;
the more general ambient-divisor spanning assertion is stated separately.
Single-germ obstructions require injectivity into the resolution homology.
The mathematical data accompanying the proofs are in [certificates/](certificates/).

## Build

From this directory:

```sh
latexmk -pdf main.tex
```

The resulting `main.pdf` is the local build output; `cy-vanishing-cycle.pdf` is the
published PDF.

## Computations

Run from the repository root:

```sh
python3 paper4/certificates/global_kernel.py
python3 paper4/certificates/relation_class.py
python3 paper4/certificates/kernel_equality.py
python3 paper4/certificates/classification.py
```

The [series reproduction guide](../README.md#reproducing) lists the complete
commands and distinguishes standard Python from SageMath. Database-wide
reproduction and the input manifest are described in [DATASET_CHECK.md](../DATASET_CHECK.md).
