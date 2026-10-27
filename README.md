# folded-sheet-fracture

Code and data for the paper

> E. Barbieri, *Cracks running into the ridge of a folded sheet* (submitted).

A straight crack runs at right angles towards the ridge of a folded thin sheet.
The theory writes the crack as distributions of dislocations in a shallow shell
and solves singular integral equations. COMSOL shell finite elements check the
theory for modes I, II and III.

All lengths are measured in the ridge length scale

    lc = h / (alpha * sqrt(3 (1 - nu^2)))

where h is the thickness, alpha the slope of each facet and nu Poisson's ratio.
(Some file names still say "fold": they mean the ridge.)

## Contents

| Folder | What is in it |
|---|---|
| `scripts/theory/` | Python solvers. `StageB_solver.py` (mode I), `ModeII_solver.py`, `ModeIII_solver.py`: singular integral equations solved by Gauss–Chebyshev collocation. `make_theory_curves.py`: all theory curves. `compare_fe.py`: finite elements against theory. `williams_fold.py`, `williams_plate_support.py`: Williams exponents at the ridge. `K_modeI.py`: stress intensity factors. `deform/`: crack-face profiles. |
| `scripts/comsol/` | COMSOL 6.4 Java files, one per finite-element test (mode I at a/lc = 10, 100, 1000; mesh and domain refinement; slope test; modes II and III). `Rooftop_Export.java` exports the fields for the 3D figure. |
| `scripts/comsol/figures/` | COMSOL Java files for the deformation figure. |
| `scripts/figs3d/` | `render_fig.py`: draws the 3D deformation figure from the exported fields. |
| `results/theory/` | Theory curves (CSV). `mode2/` and `mode3/` for modes II and III. |
| `results/comsol/` | Finite-element results: run logs (`.txt`), tables (`.csv`), comparisons with theory (`fe_vs_theory*.csv`) and exported fields (`fig_*`). |

The COMSOL model files (`.mph`, up to 250 MB) are not included. Each Java file
rebuilds its model.

## How to rerun

Theory (Python 3 with NumPy, SciPy, SymPy and Matplotlib; see `requirements.txt`):

    cd scripts/theory
    python3 make_theory_curves.py
    python3 compare_fe.py

Finite elements (COMSOL 6.4 with the Structural Mechanics Module):

    cd scripts/comsol
    comsol compile Rooftop_FE_ModeII.java
    comsol batch -inputfile Rooftop_FE_ModeII.class

Each run writes its table to `results/comsol/`. If COMSOL's security settings
block relative paths, set the absolute path to `results/comsol` in `OUT_DIRS`
near the end of the Java file, and allow file access in
Preferences → Security.

## Normalisation

The COMSOL model measures the crack length, the tip distance and the released
segment along x. The theory measures them along the facet, which is
sqrt(1 + alpha^2) times longer. `compare_fe.py` divides the finite-element
G/G0 by 1 + alpha^2 and evaluates the theory at the facet lengths.

## Licence

MIT (see `LICENSE`). If you use the code or data, please cite the paper and
the archived version of this repository (DOI on the Zenodo badge).
