"""G with the tip on the ridge against crack half-length, for modes I, II and III,
on a fine grid of a/lc (for the continuous curves in the paper).

Writes ../../results/theory/G_at_ridge_vs_a_fine.csv with columns
a_over_lc, G1 (mode I, G/G0), G2 (mode II, G/G0), G3 (mode III, G/G0).
Run from this folder:  python3 make_ridge_curves.py
"""
import os
import numpy as np
import StageB_solver as S, ModeII_solver as M2, ModeIII_solver as M3

DELTA = 1e-3   # tip distance to the ridge, in units of lc (the "tip on the ridge" limit)
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "results", "theory")

def mode3(a, d, tol=1e-4):
    n = 400; g = M3.ratio(a, d, n=n)[0]
    while n < 3200:   # n = 3200 is the memory limit of an ordinary computer
        n *= 2; g2 = M3.ratio(a, d, n=n)[0]
        if abs(g2 - g) < tol * abs(g2): return g2
        g = g2
    return g

a_list = np.unique(np.r_[np.logspace(0, np.log10(3000), 25), [1, 3, 10, 30, 100, 300, 1000, 3000]])
fn = os.path.join(OUT, "G_at_ridge_vs_a_fine.csv")
done = np.loadtxt(fn, delimiter=",", skiprows=1, ndmin=2)[:, 0] if os.path.exists(fn) else []
if not len(done):
    open(fn, "w").write("a_over_lc,G1,G2,G3\n")
for a in a_list:
    if np.any(np.isclose(done, a, rtol=1e-4)): continue   # resume an interrupted run
    row = (a, S.solve_converged(a, DELTA)[0], M2.conv(a, DELTA)[0], mode3(a, DELTA))
    print("%8.2f  %.5f  %.4f  %.5f" % row, flush=True)
    open(fn, "a").write(",".join("%.6g" % v for v in row) + "\n")
