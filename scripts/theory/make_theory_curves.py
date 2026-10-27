"""Theory curves for the rooftop-crack paper.

Uses StageB_solver.py (coupled membrane + bending singular integral
equations, Gauss-Chebyshev). All lengths are in units of the fold length
lc = h / (alpha * sqrt(3 (1 - nu^2))).

Writes CSV files to ../../results/theory/:
  G_vs_delta_a{10,100,1000}.csv  : G/G_flat, Gm/G_flat, Gb/G_flat against delta/lc
  G_at_fold_vs_a.csv             : G/G_flat at the fold (delta -> 0) against a/lc

Run from this folder:  python3 make_theory_curves.py
"""
import os
import numpy as np
import StageB_solver as S

NU = 0.3
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "results", "theory")


def curve(a_over_lc, dmin=1e-3, dmax=None, npts=40):
    """G/G_flat at the near tip for log-spaced tip-ridge distances delta/lc."""
    dmax = dmax or 3 * a_over_lc
    ds = np.logspace(np.log10(dmin), np.log10(dmax), npts)
    rows = [(d, *S.solve_converged(a_over_lc, d, nu=NU)) for d in ds]
    return np.array(rows)


def main():
    os.makedirs(OUT, exist_ok=True)
    hdr = "delta_over_lc,G_over_Gflat,Gm_over_Gflat,Gb_over_Gflat"
    for a in (10, 100, 1000):
        data = curve(a)
        np.savetxt(os.path.join(OUT, f"G_vs_delta_a{a}.csv"), data,
                   delimiter=",", header=hdr, comments="", fmt="%.6g")
        print(f"a/lc = {a}: peak G/G_flat = {data[:, 1].max():.4f} "
              f"at delta/lc = {data[data[:, 1].argmax(), 0]:.3g}; "
              f"G/G_flat at delta = 1e-3 lc: {data[0, 1]:.4f}")

    # Shielding law: G at the fold against crack size
    a_list = [1, 3, 10, 30, 100, 300, 1000, 3000]
    rows = [(a, S.solve_converged(a, 1e-3, nu=NU)[0]) for a in a_list]
    np.savetxt(os.path.join(OUT, "G_at_fold_vs_a.csv"), np.array(rows),
               delimiter=",", header="a_over_lc,G_fold_over_Gflat",
               comments="", fmt="%.6g")
    for (a1, g1), (a2, g2) in zip(rows[:-1], rows[1:]):
        print(f"a/lc {a1:>5} -> {a2:>5}: local exponent "
              f"{np.log(g2 / g1) / np.log(a2 / a1):+.3f}")


if __name__ == "__main__":
    main()
