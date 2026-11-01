"""
Sanity check: does Python's h2_theory formula (h2_pop + sigma_s^2/C_f, C_f = mean(f^2/(1+f^2)))
match the formula and numeric values produced by the R baseline, under the original AR1
rho=0.3 structure (not the new realistic-LD scenarios)?

Formula comparison:
    R (Sam_R_code.R / Sam_R_code_fast_randf.R): C_f = mean(f_vec^2/(1+f_vec^2))
                                                 res_theory = sigma_s^2 / C_f        (bias only)
    Python (do_AR1_sim_Sam.py):                 Cf = np.mean(f**2/(1+f**2))
                                                 h2_theory = h2_pop + sigma_s**2/Cf  (total)
These are the same quantity up to the h2_pop offset: Python's h2_theory should equal
h2_pop + R's res_theory.

Numeric comparison uses data/r_results/Sam_R_code_fast_randf_summary.csv (R, f redrawn every
replicate) against data/python_insample/*.csv (Python, f also redrawn every replicate via
do_AR1_sim_Sam.py) - the two use the same f-redraw-per-replicate convention and the same
Fst grid, so they are directly comparable despite being independent random draws (same n_sim
= 100, same m = 1000, same sigma_s = 0.3).

Usage: python compare_theory_AR1_vs_R.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PROJECT_ROOT.parent
DATA_DIR = REPO_ROOT / "data"

H2_POP = 0.3

# Fst label as used in python_insample filenames: str(Fst).replace(".", "")
FST_LIST = [0, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]


def fst_label(fst):
    # matches do_AR1_sim_Sam.py: Fst = np.float64(sys.argv[1]); str(Fst).replace(".", "")
    return str(np.float64(fst)).replace(".", "")


def main():
    r_summary = pd.read_csv(DATA_DIR / "r_results" / "Sam_R_code_fast_randf_summary.csv")
    # the CSV has a duplicated header row baked in as data (R write quirk) - drop it if present
    r_summary = r_summary[pd.to_numeric(r_summary["F_ST"], errors="coerce").notna()].copy()
    r_summary["F_ST"] = r_summary["F_ST"].astype(float)
    r_summary["res_theory"] = r_summary["res_theory"].astype(float)
    r_summary["res_theory_se"] = r_summary["res_theory_se"].astype(float)

    rows = []
    for fst in FST_LIST:
        py_fp = DATA_DIR / "python_insample" / "AR1_n500_m1000_Fst{lbl}_sigmas03_rho03_h2pop03.csv".format(lbl=fst_label(fst))
        if not py_fp.exists():
            print(f"skip Fst={fst}: no python_insample file at {py_fp}")
            continue
        r_row = r_summary[np.isclose(r_summary["F_ST"], fst)]
        if r_row.empty:
            print(f"skip Fst={fst}: no matching row in R randf summary")
            continue
        r_bias = r_row["res_theory"].iloc[0]
        r_bias_se = r_row["res_theory_se"].iloc[0]

        py_df = pd.read_csv(py_fp)
        py_theory_mean = py_df["h2_theory"].mean()
        py_theory_se = py_df["h2_theory"].std(ddof=1) / np.sqrt(len(py_df))
        py_bias = py_theory_mean - H2_POP

        r_total = H2_POP + r_bias
        abs_diff = py_theory_mean - r_total
        pct_diff = 100 * abs_diff / r_total if np.isfinite(r_total) and r_total != 0 else np.nan

        rows.append({
            "Fst": fst,
            "R_bias_term (sigma_s^2/Cf)": r_bias,
            "R_bias_se": r_bias_se,
            "Python_bias_term (h2_theory - h2_pop)": py_bias,
            "Python_bias_se": py_theory_se,
            "R_total (h2_pop + R_bias)": r_total,
            "Python_total (h2_theory mean)": py_theory_mean,
            "abs_diff": abs_diff,
            "pct_diff": pct_diff,
        })

    out = pd.DataFrame(rows)
    pd.set_option("display.width", 160)
    pd.set_option("display.float_format", lambda v: f"{v:0.5f}")
    print(out.to_string(index=False))

    out_path = DATA_DIR / "r_results" / "compare_theory_AR1_python_vs_R.csv"
    out.to_csv(out_path, index=False)
    print(f"\nwrote {out_path}")

    finite = out[np.isfinite(out["pct_diff"])]
    if not finite.empty:
        print(f"\nmax |pct_diff| (excluding Fst=0, where the formula is undefined): {finite['pct_diff'].abs().max():.3f}%")


if __name__ == "__main__":
    main()
