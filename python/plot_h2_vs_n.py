"""
Plots mean h2 vs. sample size n (m = 2n), at a fixed F_ST=0.01 and sigma_s=0.3, on real chr1
LD structure. One figure per scenario (no ref panel / ref panel x shared f / independent f),
each overlaying all six LDSC estimator variants. chr22 is not included - it has not been run
at n=10000, and this plot is specifically about how estimates behave as sample size scales up.

Usage: python plot_h2_vs_n.py
Writes figures/h2_vs_n_chr1_Fst001_<scenario>.png (one per scenario)
"""
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from generate_realistic_sharedf_vs_indepf_report import add_end_labels, csv_path
from generate_ldsc_estimator_variants_report import (
    ESTIMATOR_COLORS,
    ESTIMATOR_COLS,
    ESTIMATOR_LINESTYLES,
    ESTIMATOR_MARKERS,
    ESTIMATOR_SHORT_LABELS,
    SCENARIOS,
)

PROJECT_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PROJECT_ROOT.parent
FIGURES_DIR = REPO_ROOT / "figures"

CHROM = "chr1"
FST = 0.01
NS = [500, 1000, 2000, 5000, 10000]
H2_POP = 0.3


def load_point(dirname, n, col):
    fp = csv_path(CHROM, dirname, n, n * 2, FST)
    if not fp.exists():
        return None
    df = pd.read_csv(fp)
    mean = df[col].mean()
    se = df[col].std(ddof=1) / np.sqrt(len(df))
    return mean, se


def load_series(dirname):
    """dict: estimator col -> (ns, means, ses), only over n where data exists."""
    series = {}
    for col in ESTIMATOR_COLS:
        ns_have, means, ses = [], [], []
        for n in NS:
            pt = load_point(dirname, n, col)
            if pt is None:
                continue
            ns_have.append(n)
            means.append(pt[0])
            ses.append(pt[1])
        if ns_have:
            series[col] = (ns_have, means, ses)
    return series


def build_figure(series, title):
    fig, ax = plt.subplots(figsize=(8.0, 6.2))
    ax.axhline(H2_POP, linestyle="-.", color="#333333", linewidth=1.6, alpha=0.9, zorder=2.5, label="h2_pop = 0.3")
    end_label_entries = []
    for col, (ns_have, means, ses) in series.items():
        ax.errorbar(
            ns_have, means, yerr=1.96 * np.array(ses), capsize=3,
            color=ESTIMATOR_COLORS[col], marker=ESTIMATOR_MARKERS[col],
            markersize=5, markeredgewidth=0.6, markeredgecolor="black",
            linestyle=ESTIMATOR_LINESTYLES[col], label=ESTIMATOR_COLS[col],
        )
        end_label_entries.append((ESTIMATOR_SHORT_LABELS[col], means[-1], ESTIMATOR_COLORS[col]))
    ax.set_xscale("log")
    ax.set_xticks(NS)
    ax.set_xticklabels([str(n) for n in NS])
    ax.set_xlabel("n (m = 2n)")
    ax.set_ylabel("Estimated h2")
    ax.set_title(title)
    # Direct end-of-line labels (h2_pop excluded - it's a flat reference line, not an estimator
    # variant - mirroring the companion reports' F_ST-sweep figures).
    add_end_labels(ax, NS[-1], end_label_entries)
    # Shared legend outside the plot area, not overlapping data.
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=2,
        fontsize=7, frameon=False, handlelength=2.5, columnspacing=1.0, labelspacing=0.5,
    )
    fig.tight_layout(rect=[0, 0.2, 0.86, 1])
    return fig


def out_path_for(scenario_label):
    slug = re.sub(r"[^a-z0-9]+", "_", scenario_label.lower()).strip("_")
    return FIGURES_DIR / f"h2_vs_n_chr1_Fst001_{slug}.png"


def main():
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for dirname, scenario_label in SCENARIOS.items():
        series = load_series(dirname)
        fig = build_figure(series, f"Real {CHROM} LD, {scenario_label}: h2 vs. n at F_ST={FST}, sigma_s=0.3")
        out_path = out_path_for(scenario_label)
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
