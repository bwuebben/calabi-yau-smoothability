# Smoothing Calabi–Yau threefolds with nodes and del Pezzo cone points

Paper 5 · Bernd Johannes Wuebben · 61 pages

[PDF](cy-mixed-smoothing.pdf) · [LaTeX source](main.tex) · [Series overview](../README.md)

A necessary and sufficient smoothing criterion for connected normal
projective complex threefolds with trivial dualizing sheaf, H¹(O_X)=0,
and nodes or exact anticanonical cones over smooth del Pezzo surfaces of
degrees five, six and seven, or over P¹×P¹. A homological relation must avoid
every local discriminant. The selected global deformation base is smooth
when the relation space projects nontrivially to each selected degree-six
local parameter space; its tangent vectors then integrate convergently.

The paper determines the full reduced analytic deformation ring of a
nonsmoothable example with four smooth components. It proves that X₁₉
smooths and has a smooth 30-dimensional deformation base, and constructs an
ambient smoothing of X₉. Corollary 8.15 specifies positivity on the common
face and counts threshold points: at most 2m+1 sign patterns for m distinct
thresholds, including the m+1 open chambers.

The [sections/](sections/) directory contains required LaTeX inputs,
including the maintained degree-scope table. The theorem does not assert
the same mixed criterion for degrees one through four.

The bounded census uses the at-most-nine-vertex part of the fully verified
database copy. [Input verification and provenance](../DATASET_CHECK.md)
establish completeness and absence of repetitions up to lattice equivalence;
the census counts remain the original computations.

## Build

From this directory:

```sh
latexmk -pdf main.tex
```

The resulting `main.pdf` is the local build output; `cy-mixed-smoothing.pdf` is the
published PDF.

## Computations

Run from the repository root:

```sh
python3 paper5/certificates/locking.py
python3 paper5/certificates/slice_rigidity.py
python3 paper5/certificates/theoremB_class.py
python3 paper5/figures/code/compute01_degree5_local.py
python3 paper5/figures/code/compute01_degree8_local.py
python3 paper5/figures/code/compute01_cubic_scope.py
```

The [series reproduction guide](../README.md#reproducing) lists the complete
commands and distinguishes standard Python from SageMath. Database-wide
reproduction and the input manifest are described in [DATASET_CHECK.md](../DATASET_CHECK.md).
