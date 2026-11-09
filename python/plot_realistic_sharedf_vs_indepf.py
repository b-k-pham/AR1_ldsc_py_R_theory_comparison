"""
Plot shared-f vs. independent-f heritability estimates under REAL (1000G chr1) LD
structure, in two settings: no reference panel, and with a reference panel.

This is Python-only: unlike the AR1 pipeline, there is no R-side equivalent for the
realistic-LD scenarios, so there is nothing to cross-validate against R here.

Reads the four data/python_realistic_* directories written by:
    do_realistic_sim_Sam_noref_sharedf.py   (no ref panel, f shared across replicates)
    do_realistic_sim_Sam_noref_indepf.py    (no ref panel, f redrawn every replicate)
    do_realistic_sim_Sam_ref_fixedf.py      (ref panel, f shared with main replicates)
    do_realistic_sim_Sam_ref_indepf.py      (ref panel, f independent of main replicates)

Produces two figures (n/m suffix, default n=500/m=1000):
    figures/h2_by_fst_realistic_noref_sharedf_vs_indepf_n{n}_m{m}.png
    figures/h2_by_fst_realistic_refpanel_fixedf_vs_indepf_n{n}_m{m}.png

Usage: python plot_realistic_sharedf_vs_indepf.py [n]   (m is always n*2; default n=500)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PROJECT_ROOT.parent
DATA_DIR = REPO_ROOT / "data"
FIGURES_DIR = REPO_ROOT / "figures"

# Per feedback: only F_ST 0-0.1 is relevant, so the "full range" panel now stops at 0.1
# (0.1-0.5 dropped entirely) and the zoomed panel goes further in, to 0-0.01 in 0.001 steps.
FST_FINE = [0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01]
FST_FULL = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
H2_POP = 0.3

SCENARIOS = {
    "no ref panel": {
        "shared f": "python_realistic_noref_sharedf",
        "independent f": "python_realistic_noref_indepf",
    },
    "ref panel": {
        "shared f": "python_realistic_refpanel_fixedf",
        "independent f": "python_realistic_refpanel_indepf",
    },
}

ESTIMATOR_COLS = {
    "h2_ldsc_reg": "LDSC (free intercept)",
    "h2_gwash": "GWASH",
}


def fst_label(fst):
    # matches driver scripts: Fst = np.float64(sys.argv[1]); str(Fst).replace(".", "") -
    # plain str(fst) gives "0" for Fst=0 instead of the actual "00" filename label.
    return str(np.float64(fst)).replace(".", "")


def load_summary(dirname, n, m, fst_list):
    """Read every available Fst CSV in a data/<dirname>/ folder and return one row per
    Fst with mean and standard-error-of-the-mean for each estimator in ESTIMATOR_COLS,
    plus the mean theoretical h2."""
    rows = []
    for fst in fst_list:
        fp = DATA_DIR / dirname / "realisticLD_n{n}_m{m}_Fst{lbl}_sigmas03_h2pop03.csv".format(n=n, m=m, lbl=fst_label(fst))
        if not fp.exists():
            continue
        df = pd.read_csv(fp)
        row = {"Fst": fst, "n_reps": len(df)}
        for col in ESTIMATOR_COLS:
            row[col + "_mean"] = df[col].mean()
            row[col + "_se"] = df[col].std(ddof=1) / np.sqrt(len(df))
        row["h2_theory_mean"] = df["h2_theory"].mean()
        rows.append(row)
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows)


def build_long_df(scenario_dict, n, m, fst_list):
    """scenario_dict: {"shared f": dirname, "independent f": dirname} -> long-format
    DataFrame with columns [Fst, f_sharing, estimator, mean_h2, se_h2, h2_theory_mean]."""
    pieces = []
    for f_sharing_label, dirname in scenario_dict.items():
        summary = load_summary(dirname, n, m, fst_list)
        if summary.empty:
            print(f"warning: no data found yet in data/{dirname}/")
            continue
        for col, estimator_label in ESTIMATOR_COLS.items():
            piece = summary[["Fst", "h2_theory_mean", col + "_mean", col + "_se"]].copy()
            piece = piece.rename(columns={col + "_mean": "mean_h2", col + "_se": "se_h2"})
            piece["estimator"] = estimator_label
            piece["f_sharing"] = f_sharing_label
            pieces.append(piece)
    if not pieces:
        return pd.DataFrame()
    return pd.concat(pieces, axis=0, ignore_index=True)


FOCUS_ESTIMATOR = "LDSC (free intercept)"


def make_plot(fine_long_df, full_long_df, title, out_path, n, m, focus_estimator=FOCUS_ESTIMATOR):
    colors = {"shared f": "#1b9e77", "independent f": "#d95f02"}
    fine_sub = fine_long_df[fine_long_df["estimator"] == focus_estimator].sort_values("Fst")
    full_sub = full_long_df[full_long_df["estimator"] == focus_estimator].sort_values("Fst")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    panel_specs = [
        (axes[0], fine_sub, "Zoomed: F_ST 0 to 0.01", 0.001),
        (axes[1], full_sub, "F_ST 0 to 0.1", 0.01),
    ]

    for ax, sub, panel_title, tick_step in panel_specs:
        ax.axhline(H2_POP, linestyle="--", color="grey", linewidth=1, label="h2_pop")
        # y-limit covers both the observed estimates AND the theory line (no clipping) for
        # every Fst value shown in this panel.
        observed_max = (sub["mean_h2"] + 1.96 * sub["se_h2"]).max() if not sub.empty else 1
        # At Fst=0 (f fixed at 0), C_f = mean(f^2/(1+f^2)) = 0, so the theory formula's
        # sigma_s^2/C_f term divides by zero and is +inf - a formula artifact at that one
        # point (no real stratification signal when f=0), not a real divergence. Excluded
        # from the y-limit so it doesn't blow up the axis; the point itself just doesn't
        # render for that one Fst value.
        finite_theory = sub["h2_theory_mean"][np.isfinite(sub["h2_theory_mean"])] if not sub.empty else pd.Series([], dtype=float)
        theory_max = finite_theory.max() if not finite_theory.empty else observed_max
        panel_max = max(observed_max, theory_max)
        for f_sharing, color in colors.items():
            s = sub[sub["f_sharing"] == f_sharing]
            if s.empty:
                continue
            ax.errorbar(
                s["Fst"], s["mean_h2"], yerr=1.96 * s["se_h2"],
                marker="o", color=color, label=f_sharing, capsize=3,
            )
            # theory line: h2_pop + sigma_s^2/C_f, using this scenario's own per-replicate
            # f draws (same formula as Sam_R_code.R's C_f = mean(f^2/(1+f^2))). The y-limit
            # below is sized to include this line in full (it diverges as F_ST -> 0, so it
            # dominates the panel's lowest Fst value, but is never cut off).
            ax.plot(
                s["Fst"], s["h2_theory_mean"],
                linestyle=":", color=color, alpha=0.7, linewidth=1.5,
                label=f"{f_sharing} (theory)",
            )
        ax.set_ylim(0, panel_max * 1.1)
        ax.set_xlabel("F_ST")
        ax.set_ylabel("Estimated h2")
        ax.set_title(panel_title)
        ax.set_xticks(np.arange(0, sub["Fst"].max() + 1e-9, tick_step) if not sub.empty else [0])
        ax.tick_params(axis="x", labelrotation=45)
        ax.legend(title="Stratification vector f", fontsize=8)

    fig.suptitle(f"{title} ({focus_estimator})\nn={n}, m={m}")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"wrote {out_path}")


def main():
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    n = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    m = n * 2
    suffix = f"_n{n}_m{m}"

    noref_fine = build_long_df(SCENARIOS["no ref panel"], n, m, FST_FINE)
    noref_full = build_long_df(SCENARIOS["no ref panel"], n, m, FST_FULL)
    if not noref_fine.empty and not noref_full.empty:
        make_plot(
            noref_fine, noref_full,
            "Real chr1 LD, no reference panel: shared vs. independent f",
            FIGURES_DIR / f"h2_by_fst_realistic_noref_sharedf_vs_indepf{suffix}.png",
            n, m,
        )

    refpanel_fine = build_long_df(SCENARIOS["ref panel"], n, m, FST_FINE)
    refpanel_full = build_long_df(SCENARIOS["ref panel"], n, m, FST_FULL)
    if not refpanel_fine.empty and not refpanel_full.empty:
        make_plot(
            refpanel_fine, refpanel_full,
            "Real chr1 LD, with reference panel: shared vs. independent f",
            FIGURES_DIR / f"h2_by_fst_realistic_refpanel_fixedf_vs_indepf{suffix}.png",
            n, m,
        )


if __name__ == "__main__":
    main()
