# Deformations of toric pairs and the smoothing of Batyrev Calabi–Yau threefolds

Paper 2 · Bernd Johannes Wuebben · 19 pages

[PDF](cy-toric-pairs.pdf) · [LaTeX source](main.tex) · [Series overview](../README.md)

Trinomial deformations of toric pairs smooth the points over a chosen
smoothable two-face when the dual edge has an interior lattice point.
This smooths the whole hypersurface when it is the only singular face.
Simultaneous smoothing is unconditional for nodes and otherwise assumes
irreducibility of the semiuniversal deformation space. A one-node example
shows the dual-edge threshold is sharp.

The lone del Pezzo cone cases use the explicitly stated link-homology
necessary condition from [Paper 5](../paper5/README.md). The comparison of
3,774 facewise cases with 30,241 Batyrev–Kreuzer smoothings retains the
manuscript's database hypothesis.

## Build

From this directory:

```sh
latexmk -pdf main.tex
```

The resulting `main.pdf` is the local build output; `cy-toric-pairs.pdf` is the
published PDF.

## Computations

Run from the repository root:

```sh
python3 src/paper2_check.py
python3 src/cascade_check.py
```

The [series reproduction guide](../README.md#reproducing) lists the complete
commands and distinguishes standard Python from SageMath. Database-wide
reproduction and the input manifest are described in [DATASET_CHECK.md](../DATASET_CHECK.md).
