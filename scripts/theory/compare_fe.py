"""Compare COMSOL shell results with the integral-equation theory.

Reads the batch logs (*.txt, written with `tee`) and the older CSV tables
in ../../results/comsol/. For each case it:
  1. divides G_FE/G_flat by the flat-sheet run of the same set
     (the 'FLAT CALIBRATION' line that starts the set);
  2. divides by 1 + alpha^2: the model measures a and the released segment
     Delta along x, the theory along the facet (each is sqrt(1+alpha^2)
     times longer);
  3. evaluates the theory at the midpoint of the released segment,
     delta_mid = (delta - Delta/2) * sqrt(1 + alpha^2), for the unfolded
     crack size a * sqrt(1 + alpha^2).
Prints FE, theory and the difference, and writes fe_vs_theory.csv in
../../results/comsol/.

Run from this folder:  python3 compare_fe.py
"""
import glob
import os
import re
import numpy as np
import StageB_solver as S

HERE = os.path.dirname(os.path.abspath(__file__))
FE = os.path.join(HERE, "..", "..", "results", "comsol")

RE_CASE = re.compile(r"---- alpha=([\d.]+)\s+delta/lc=([\d.]+)\s+(FLAT CALIBRATION)?.*dx/lc=([\d.]+)")
RE_G = re.compile(r"G_energy/G_flat = ([-\d.]+)\s+Gm/G_flat = ([-\d.]+)\s+Gb/G_flat = ([-\d.]+)")
RE_A = re.compile(r"a/lc = ([\d.]+)")
RE_SET = re.compile(r"(REFINE_TEST|DECOUPLE_TEST|A100_TEST|SHEAR_TEST).*")


def read_log(path):
    """Yield dicts: one per non-flat case, with its set's flat calibration."""
    a_over_lc, label, flat, case = 10.0, "", None, None
    for line in open(path, errors="replace"):
        m = RE_SET.search(line)
        if m:
            label = line.strip()
            ma = RE_A.search(line)
            a_over_lc = float(ma.group(1)) if ma else 10.0
            flat = None
            continue
        m = RE_CASE.search(line)
        if m:
            case = dict(alpha=float(m.group(1)), delta=float(m.group(2)),
                        is_flat=bool(m.group(3)), dx=float(m.group(4)))
            continue
        m = RE_G.search(line)
        if m and case:
            g, gb = float(m.group(1)), float(m.group(3))
            if case["is_flat"]:
                flat = g
            elif flat:
                yield dict(file=os.path.basename(path), set=label, a=a_over_lc,
                           G=g, Gb=gb, flat=flat, **case)
            case = None


def read_csv(path):
    """Older COMSOL table export (19 Oct): one flat row per alpha."""
    rows = [l.strip().split(",") for l in open(path) if l.strip() and not l.startswith("%")]
    flat = {}
    for r in rows:
        if r[6] == "1":
            flat[float(r[0])] = float(r[9])
    for r in rows:
        if r[6] == "0":
            al = float(r[0])
            yield dict(file=os.path.basename(path), set="", a=float(r[2]), alpha=al,
                       delta=float(r[3]), dx=float(r[4]), is_flat=False,
                       G=float(r[9]), Gb=float(r[11]), flat=flat[al])


def main():
    cases = []
    for p in sorted(glob.glob(os.path.join(FE, "*.txt"))):
        cases += list(read_log(p))
    for p in sorted(glob.glob(os.path.join(FE, "*.csv"))):
        if "fe_vs_theory" not in p:
            cases += list(read_csv(p))

    out = ["file,set,alpha,a_over_lc,delta_mid_over_lc,G_FE,G_theory,diff_percent,Gb_FE"]
    for c in cases:
        s = np.sqrt(1 + c["alpha"] ** 2)
        g_fe = c["G"] / c["flat"] / s ** 2
        d_mid = (c["delta"] - c["dx"] / 2) * s
        g_th = S.solve_converged(c["a"] * s, d_mid)[0]
        diff = 100 * (g_fe / g_th - 1)
        print(f"{c['file'][:44]:44s} a={c['alpha']:<5} a/lc={c['a']:<6g} "
              f"d={d_mid:<8.4g} FE={g_fe:.4f} th={g_th:.4f} {diff:+.2f}%")
        out.append(f"{c['file']},\"{c['set']}\",{c['alpha']},{c['a']},{d_mid:.5g},"
                   f"{g_fe:.5f},{g_th:.5f},{diff:.3f},{c['Gb'] / c['flat'] / s**2:.5f}")
    with open(os.path.join(FE, "fe_vs_theory.csv"), "w") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
