"""
Plot shared-f vs. independent-f heritability estimates under REAL (1000G chr22) LD
structure, in two settings: no reference panel, and with a reference panel.

Chr22-specific twin of plot_realistic_sharedf_vs_indepf.py (which plots chr1): reads the
chr22-suffixed data directories/filenames written by the _chr22 driver scripts, and writes
chr22-suffixed figures, so chr1 and chr22 results stay fully separate and traceable.

This is Python-only: unlike the AR1 pipeline, there is no R-side equivalent for the
realistic-LD scenarios, so there is nothing to cross-validate against R here.

Reads the four data/python_realistic_*_chr22 directories written by:
    do_realistic_sim_Sam_noref_sharedf_chr22.py   (no ref panel, f shared across replicates)
    do_realistic_sim_Sam_noref_indepf_chr22.py    (no ref panel, f redrawn every replicate)
    do_realistic_sim_Sam_ref_fixedf_chr22.py      (ref panel, f shared with main replicates)
    do_realistic_sim_Sam_ref_indepf_chr22.py      (ref panel, f independent of main replicates)

Produces two figures (n/m suffix, default n=500/m=1000):
    figures/h2_by_fst_realistic_chr22_noref_sharedf_vs_indepf_n{n}_m{m}.png
    figures/h2_by_fst_realistic_chr22_refpanel_fixedf_vs_indepf_n{n}_m{m}.png

Usage: python plot_realistic_sharedf_vs_indepf_chr22.py [n]   (m is always n*2; default n=500)
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

FST_LIST = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
H2_POP = 0.3

SCENARIOS = {
    "no ref panel": {
        "shared f": "python_realistic_noref_sharedf_chr22",
        "independent f": "python_realistic_noref_indepf_chr22",
    },
    "ref panel": {
        "shared f": "python_realistic_refpanel_fixedf_chr22",
        "independent f": "python_realistic_refpanel_indepf_chr22",
    },
}

ESTIMATOR_COLS = {
    "h2_ldsc_reg": "LDSC (free intercept)",
    "h2_gwash": "GWASH",
}


def fst_label(fst):
    return str(fst).replace(".", "")


def load_summary(dirname, n, m):
    """Read every available Fst CSV in a data/<dirname>/ folder and return one row per
    Fst with mean and standard-error-of-the-mean for each estimator in ESTIMATOR_COLS,
    plus the mean theoretical h2."""
    rows = []
    for fst in FST_LIST:
        fp = DATA_DIR / dirname / "realisticLD_chr22_n{n}_m{m}_Fst{lbl}_sigmas03_h2pop03.csv".format(n=n, m=m, lbl=fst_label(fst))
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


def build_long_df(scenario_dict, n, m):
    """scenario_dict: {"shared f": dirname, "independent f": dirname} -> long-format
    DataFrame with columns [Fst, f_sharing, estimator, mean_h2, se_h2, h2_theory_mean]."""
    pieces = []
    for f_sharing_label, dirname in scenario_dict.items():
        summary = load_summary(dirname, n, m)
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
ZOOM_FST_MAX = 0.05


def make_plot(long_df, title, out_path, n, m, focus_estimator=FOCUS_ESTIMATOR, zoom_fst_max=ZOOM_FST_MAX):
    colors = {"shared f": "#1b9e77", "independent f": "#d95f02"}
    sub_all = long_df[long_df["estimator"] == focus_estimator].sort_values("Fst")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    panel_specs = [
        (axes[0], sub_all[sub_all["Fst"] <= zoom_fst_max], f"Zoomed: F_ST 0 to {zoom_fst_max}"),
        (axes[1], sub_all, "Full range: F_ST 0 to 0.5"),
    ]

    for ax, sub, panel_title in panel_specs:
        ax.axhline(H2_POP, linestyle="--", color="grey", linewidth=1, label="h2_pop")
        # y-limit covers both the observed estimates AND the theory line (no clipping) for
        # every Fst value shown in this panel.
        observed_max = (sub["mean_h2"] + 1.96 * sub["se_h2"]).max() if not sub.empty else 1
        theory_max = sub["h2_theory_mean"].max() if not sub.empty else observed_max
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
        if panel_title.startswith("Zoomed"):
            ax.set_xticks(np.arange(0, zoom_fst_max + 1e-9, 0.01))
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

    noref_long = build_long_df(SCENARIOS["no ref panel"], n, m)
    if not noref_long.empty:
        make_plot(
            noref_long,
            "Real chr22 LD, no reference panel: shared vs. independent f",
            FIGURES_DIR / f"h2_by_fst_realistic_chr22_noref_sharedf_vs_indepf{suffix}.png",
            n, m,
        )

    refpanel_long = build_long_df(SCENARIOS["ref panel"], n, m)
    if not refpanel_long.empty:
        make_plot(
            refpanel_long,
            "Real chr22 LD, with reference panel: shared vs. independent f",
            FIGURES_DIR / f"h2_by_fst_realistic_chr22_refpanel_fixedf_vs_indepf{suffix}.png",
            n, m,
        )


if __name__ == "__main__":
    main()
