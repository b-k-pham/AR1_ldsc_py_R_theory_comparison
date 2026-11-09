"""
Builds report/realistic_LD_sharedf_vs_indepf_summary.pdf from the data/python_realistic_*
CSVs: shared-f vs. independent-f heritability estimates under real (1000G chr1/chr22) LD
structure, with and without a reference panel, across every available (chromosome, scale).

This PDF previously existed in report/ with no generator script anywhere in the repo (it must
have been built in an earlier, uncommitted session). This script reconstructs it from the
underlying data so that adding a new scale - e.g. n=10000/m=20000, added here - is just a
matter of re-running this script once that scale's CSVs exist in data/, instead of
hand-editing a PDF. A scale/scenario with no data for a given (chrom, n) is skipped for
figures and shown as "-" in tables, so this stays usable mid-sweep.

Note: as of n=10000, chr22 has not been run at that scale (chr1-only), and the n=10000 data
only covers the fine 0-0.1 F_ST grid (not the older 0.2/0.3/0.4/0.5 points used in Table 1),
so n=10000 only appears in Table 2 and the figures sections, not Table 1.

Usage: python generate_realistic_sharedf_vs_indepf_report.py
"""
from io import BytesIO
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    Image, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle,
)

# ReportLab's base-14 fonts fall back to the (unembedded) Symbol font for Greek letters like
# sigma, which some PDF renderers can't resolve and silently render blank. Embedding DejaVu Sans
# (bundled with matplotlib, so no extra dependency) avoids that - used wherever this report needs
# an actual sigma glyph rather than the word "sigma".
pdfmetrics.registerFont(TTFont("DejaVuSans", str(Path(matplotlib.get_data_path()) / "fonts" / "ttf" / "DejaVuSans.ttf")))
SIGMA_S2 = '<font name="DejaVuSans">σ</font><sub>s</sub><super>2</super>'

PROJECT_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PROJECT_ROOT.parent
DATA_DIR = REPO_ROOT / "data"
REPORT_DIR = REPO_ROOT / "report"
OUT_PDF = REPORT_DIR / "realistic_LD_sharedf_vs_indepf_summary.pdf"

H2_POP = 0.3
FST_FINE = [0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01]
FST_FULL = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
FST_WIDE = [0.3, 0.5]       # Table 1 - from the older, wider grid (not run at n=10000)
FST_LOW = [0, 0.005, 0.01]  # Table 2

ESTIMATOR_COLS = {"h2_ldsc_reg": "LDSC (free intercept)", "h2_gwash": "GWASH"}
FOCUS_ESTIMATOR = "LDSC (free intercept)"

# All six LDSC estimator variants, for the shared-f/independent-f figures below. Single source
# of truth for this encoding - generate_ldsc_estimator_variants_report.py imports these too, so
# both reports stay visually consistent. Color encodes intercept type (red=free, pink=fixed),
# using the exact hex values from the user's previously-published figures (LDSC Free/Fixed
# Intercept legend swatches: #FF0000 / #FB93C4) for visual continuity; linestyle encodes SNP
# weighting (solid=weighted, dashed=unweighted); marker shape encodes estimation mode
# (circle=two-step, triangle=single-step, square=fixed-intercept, which has no re-estimation
# step). Each of the three channels maps to exactly one concept, so no two variants share all
# three.
LDSC_VARIANT_COLS = {
    "h2_ldsc_reg": "LDSC free intercept, weighted, two-step",
    "h2_ldsc_fixed": "LDSC fixed intercept, weighted",
    "h2_ldsc_fixed_nw": "LDSC fixed intercept, unweighted",
    "h2_ldsc_nw": "LDSC free intercept, unweighted, two-step",
    "h2_ldsc_nts": "LDSC free intercept, weighted, single-step",
    "h2_ldsc_nwts": "LDSC free intercept, unweighted, single-step",
}
LDSC_VARIANT_COLORS = {
    "h2_ldsc_reg": "#FF0000",       # free intercept (red)
    "h2_ldsc_nw": "#FF0000",        # free intercept (red)
    "h2_ldsc_nts": "#FF0000",       # free intercept (red)
    "h2_ldsc_nwts": "#FF0000",      # free intercept (red)
    "h2_ldsc_fixed": "#FB93C4",     # fixed intercept (pink)
    "h2_ldsc_fixed_nw": "#FB93C4",  # fixed intercept (pink)
}
LDSC_VARIANT_LINESTYLES = {
    "h2_ldsc_reg": "-",        # weighted
    "h2_ldsc_fixed": "-",      # weighted
    "h2_ldsc_nts": "-",        # weighted
    "h2_ldsc_nw": "--",        # unweighted
    "h2_ldsc_fixed_nw": "--",  # unweighted
    "h2_ldsc_nwts": "--",      # unweighted
}
LDSC_VARIANT_MARKERS = {
    "h2_ldsc_reg": "o",        # two-step
    "h2_ldsc_nw": "o",         # two-step
    "h2_ldsc_nts": "^",        # single-step
    "h2_ldsc_nwts": "^",       # single-step
    "h2_ldsc_fixed": "s",      # fixed intercept
    "h2_ldsc_fixed_nw": "s",   # fixed intercept
}
THEORY_LABEL = r"theory (h2_pop + $\sigma_{s}^{2}$/C_f)"

# Compact end-of-line labels (see add_end_labels below) - short enough to sit directly next to
# each line without a separate legend lookup, still self-decodable from the three dimensions
# (intercept type, weighting, step) rather than needing a key.
LDSC_VARIANT_SHORT_LABELS = {
    "h2_ldsc_reg": "Free·Wt·2S",
    "h2_ldsc_fixed": "Fixed·Wt",
    "h2_ldsc_fixed_nw": "Fixed·Unwt",
    "h2_ldsc_nw": "Free·Unwt·2S",
    "h2_ldsc_nts": "Free·Wt·1S",
    "h2_ldsc_nwts": "Free·Unwt·1S",
}


def add_end_labels(ax, last_x, entries, min_gap=0.055):
    """Direct end-of-line labels at the right edge of `ax`, vertically de-collided so close
    lines don't overlap. entries: list of (label, y_value, color)."""
    entries = sorted(entries, key=lambda e: -e[1])
    y0, y1 = ax.get_ylim()
    prev_frac = None
    placed = []
    for label, y, color in entries:
        frac = min(max((y - y0) / (y1 - y0), 0.0), 1.0)
        if prev_frac is not None and prev_frac - frac < min_gap:
            frac = prev_frac - min_gap
        prev_frac = frac
        placed.append((label, y, frac, color))
    for label, y, frac, color in placed:
        ax.annotate(
            label, xy=(last_x, y), xycoords="data",
            xytext=(1.02, frac), textcoords="axes fraction",
            fontsize=6, color=color, va="center", ha="left",
            arrowprops=dict(arrowstyle="-", color=color, lw=0.5, alpha=0.5, shrinkA=0, shrinkB=2),
            annotation_clip=False,
        )

CHROMS = {
    "chr1": {"dir_suffix": "", "fname": "realisticLD_n{n}_m{m}_Fst{lbl}_sigmas03_h2pop03.csv"},
    "chr22": {"dir_suffix": "_chr22", "fname": "realisticLD_chr22_n{n}_m{m}_Fst{lbl}_sigmas03_h2pop03.csv"},
}

# n=10000/m=20000 added here, chr1 only - chr22 has not been run at this scale.
SCALES_BY_CHROM = {
    "chr1": [(500, 1000), (1000, 2000), (2000, 4000), (5000, 10000), (10000, 20000)],
    "chr22": [(500, 1000), (1000, 2000), (2000, 4000), (5000, 10000)],
}
# Scales with the older wide grid (0.2-0.5) available, for Table 1.
WIDE_GRID_SCALES = {
    "chr1": [(500, 1000), (1000, 2000), (2000, 4000), (5000, 10000)],
    "chr22": [(500, 1000), (1000, 2000), (2000, 4000), (5000, 10000)],
}

SCENARIO_DIRS = {
    "no ref panel": {"shared f": "python_realistic_noref_sharedf", "independent f": "python_realistic_noref_indepf"},
    "ref panel": {"shared f": "python_realistic_refpanel_fixedf", "independent f": "python_realistic_refpanel_indepf"},
}


def fst_label(fst):
    return str(np.float64(fst)).replace(".", "")


def csv_path(chrom, dirname, n, m, fst):
    spec = CHROMS[chrom]
    return DATA_DIR / (dirname + spec["dir_suffix"]) / spec["fname"].format(n=n, m=m, lbl=fst_label(fst))


def load_summary(chrom, dirname, n, m, fst_list):
    rows = []
    for fst in fst_list:
        fp = csv_path(chrom, dirname, n, m, fst)
        if not fp.exists():
            continue
        df = pd.read_csv(fp)
        row = {"Fst": fst, "n_reps": len(df)}
        for col in ESTIMATOR_COLS:
            row[col + "_mean"] = df[col].mean()
            row[col + "_se"] = df[col].std(ddof=1) / np.sqrt(len(df))
        row["h2_theory_mean"] = df["h2_theory"].mean()
        rows.append(row)
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def load_ldsc_variant_summary(chrom, dirname, n, m, fst_list):
    rows = []
    for fst in fst_list:
        fp = csv_path(chrom, dirname, n, m, fst)
        if not fp.exists():
            continue
        df = pd.read_csv(fp)
        row = {"Fst": fst}
        for col in LDSC_VARIANT_COLS:
            row[col + "_mean"] = df[col].mean()
            row[col + "_se"] = df[col].std(ddof=1) / np.sqrt(len(df))
        row["h2_theory_mean"] = df["h2_theory"].mean()
        rows.append(row)
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def ldsc_free_mean(chrom, dirname, n, m, fst):
    fp = csv_path(chrom, dirname, n, m, fst)
    if not fp.exists():
        return None
    return pd.read_csv(fp)["h2_ldsc_reg"].mean()


def count_all_csvs():
    total = 0
    for chrom, scales in SCALES_BY_CHROM.items():
        for n, m in scales:
            for scen in SCENARIO_DIRS.values():
                for dirname in scen.values():
                    d = DATA_DIR / (dirname + CHROMS[chrom]["dir_suffix"])
                    if not d.exists():
                        continue
                    pattern = f"realisticLD_{'chr22_' if chrom == 'chr22' else ''}n{n}_m{m}_Fst*_sigmas03_h2pop03.csv"
                    total += len(list(d.glob(pattern)))
    return total


# ---------- figures ----------

def make_plot_image(chrom, dirname, n, m, title_prefix):
    fine = load_ldsc_variant_summary(chrom, dirname, n, m, FST_FINE)
    full = load_ldsc_variant_summary(chrom, dirname, n, m, FST_FULL)
    if fine.empty or full.empty:
        return None

    # Fixed (not data-driven); theory lines run off the top of the chart where they exceed it.
    # Panels use different maxima - the left (zoomed) panel's data stays well under 1 at this
    # F_ST range, so it gets its own tighter scale rather than sharing the right panel's 4.5.
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    panel_specs = [
        (axes[0], fine, "Zoomed: F_ST 0 to 0.01", 0.001, 1),
        (axes[1], full, "F_ST 0 to 0.1", 0.01, 4.5),
    ]
    for ax, sub, panel_title, tick_step, panel_y_max in panel_specs:
        ax.set_axisbelow(True)
        ax.grid(axis="y", color="#e2e2e2", linewidth=0.6, zorder=0)
        ax.axhline(H2_POP, linestyle="-.", color="#333333", linewidth=1.6, alpha=0.9, zorder=2.5, label="h2_pop")
        for col, label in LDSC_VARIANT_COLS.items():
            ax.errorbar(
                sub["Fst"], sub[col + "_mean"], yerr=1.96 * sub[col + "_se"],
                marker=LDSC_VARIANT_MARKERS[col], markersize=5, markeredgewidth=0.6, markeredgecolor="black",
                color=LDSC_VARIANT_COLORS[col], linestyle=LDSC_VARIANT_LINESTYLES[col], label=label,
                capsize=0, linewidth=1.1, elinewidth=0.6, alpha=0.9,
            )
        ax.plot(sub["Fst"], sub["h2_theory_mean"], linestyle=":", color="black", alpha=0.5, linewidth=1.1, label=THEORY_LABEL)
        ax.set_ylim(0, panel_y_max)
        ax.set_xlabel("F_ST")
        ax.set_ylabel("Estimated h2")
        ax.set_title(panel_title)
        ax.set_xticks(np.arange(0, sub["Fst"].max() + 1e-9, tick_step))
        ax.tick_params(axis="x", labelrotation=45)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

    # Direct end-of-line labels on the full-range (right) panel, so a reader scanning that panel
    # doesn't need to cross-reference the legend below to know which line is which. The zoomed
    # (left) panel is too crowded for this (lines nearly coincide near F_ST=0), so it still
    # relies on the shared legend. h2_pop and theory are left unlabeled here - they're reference
    # lines, not estimator variants, and labeling them added noise without disambiguation value.
    end_label_entries = [
        (LDSC_VARIANT_SHORT_LABELS[col], full[col + "_mean"].iloc[-1], LDSC_VARIANT_COLORS[col])
        for col in LDSC_VARIANT_SHORT_LABELS
    ]
    add_end_labels(axes[1], full["Fst"].iloc[-1], end_label_entries)

    # One shared legend below both panels (not duplicated per-panel, not overlapping data).
    # handlelength is widened so dash/dash-dot/dotted patterns are actually legible at this size.
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=4,
        fontsize=7, frameon=False, handlelength=3.0, columnspacing=1.3, labelspacing=0.5,
    )
    fig.suptitle(f"{title_prefix}\nn={n}, m={m}")
    fig.tight_layout(rect=[0, 0.16, 0.87, 1])
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=150)
    plt.close(fig)
    buf.seek(0)
    return buf


# ---------- tables ----------

def fmt(v, nd=3, signed=False):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "-"
    return f"{v:+.{nd}f}" if signed else f"{v:.{nd}f}"


def build_wide_gap_table():
    header = ["Chrom", "n", "m", "Setting",
              "Shared f\n(F_ST=0.3)", "Indep f\n(F_ST=0.3)",
              "Shared f\n(F_ST=0.5)", "Indep f\n(F_ST=0.5)", "Gap\n(F_ST=0.5)"]
    rows = [header]
    for chrom, scales in WIDE_GRID_SCALES.items():
        for n, m in scales:
            for setting, scen in SCENARIO_DIRS.items():
                s03_shared = ldsc_free_mean(chrom, scen["shared f"], n, m, 0.3)
                s03_indep = ldsc_free_mean(chrom, scen["independent f"], n, m, 0.3)
                s05_shared = ldsc_free_mean(chrom, scen["shared f"], n, m, 0.5)
                s05_indep = ldsc_free_mean(chrom, scen["independent f"], n, m, 0.5)
                gap = None if (s05_shared is None or s05_indep is None) else s05_indep - s05_shared
                rows.append([
                    chrom, str(n), str(m), setting,
                    fmt(s03_shared), fmt(s03_indep), fmt(s05_shared), fmt(s05_indep), fmt(gap, signed=True),
                ])
    return rows


def build_low_fst_table():
    header = ["Chrom", "n", "Setting",
              "Shared f\n(F_ST=0)", "Indep f\n(F_ST=0)",
              "Shared f\n(F_ST=0.005)", "Indep f\n(F_ST=0.005)",
              "Shared f\n(F_ST=0.01)", "Indep f\n(F_ST=0.01)"]
    rows = [header]
    for chrom, scales in SCALES_BY_CHROM.items():
        for n, m in scales:
            for setting, scen in SCENARIO_DIRS.items():
                vals = []
                for fst in FST_LOW:
                    vals.append(ldsc_free_mean(chrom, scen["shared f"], n, m, fst))
                    vals.append(ldsc_free_mean(chrom, scen["independent f"], n, m, fst))
                rows.append([chrom, str(n), setting] + [fmt(v) for v in vals])
    return rows


def styled_table(data, col_widths=None):
    t = Table(data, colWidths=col_widths, repeatRows=1)
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1F3864")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 7.5),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F2F2F2")]),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#CCCCCC")),
        ("ALIGN", (1, 0), (-1, -1), "CENTER"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    return t


# ---------- document ----------

def main():
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    styles = getSampleStyleSheet()
    h1 = ParagraphStyle("H1", parent=styles["Heading1"], fontSize=18, spaceAfter=4)
    h2 = ParagraphStyle("H2", parent=styles["Heading2"], fontSize=15, spaceBefore=6, spaceAfter=6)
    h3 = ParagraphStyle("H3", parent=styles["Heading3"], fontSize=12, spaceBefore=4, spaceAfter=4)
    sub = ParagraphStyle("Sub", parent=styles["Normal"], fontSize=9, textColor=colors.HexColor("#444444"), spaceAfter=12)
    body = ParagraphStyle("Body", parent=styles["Normal"], fontSize=9.5, leading=13, spaceAfter=8)
    caption = ParagraphStyle("Caption", parent=styles["Normal"], fontSize=8.5, textColor=colors.HexColor("#666666"), spaceBefore=4)
    table_title = ParagraphStyle("TableTitle", parent=styles["Normal"], fontSize=11, leading=14, alignment=1, spaceAfter=4)

    story = []

    total_csvs = count_all_csvs()

    # --- Title / Executive Summary ---
    story.append(Paragraph("Shared vs. Independent Stratification Vector (f) Under Real LD Structure", h1))
    story.append(Paragraph("A cross-scale, cross-chromosome comparison using 1000 Genomes chr1 and chr22 LD structure", sub))
    story.append(Paragraph("Executive Summary", h2))
    story.append(Paragraph(
        "This analysis extends the project's AR1-based heritability simulation framework to use "
        "<b>real linkage-disequilibrium (LD) structure</b> drawn from 1000 Genomes chromosome 1 and "
        "chromosome 22 genotypes (in place of synthetic AR1 correlation), and tests how a reference "
        "panel's stratification vector f behaves when it is <b>shared</b> with the main sample versus "
        "drawn <b>independently</b>. The comparison was run across five sample/SNP scales on chr1 "
        "(n=500/m=1000, n=1000/m=2000, n=2000/m=4000, n=5000/m=10000, n=10000/m=20000) and four scales "
        f"on chr22 (n=500 through n=5000 - chr22 has not yet been run at n=10000), for a total of "
        f"{total_csvs} simulated datasets (100 replicates each) across all scenarios, scales and chromosomes.",
        body))
    story.append(Paragraph(
        "<b>Headline finding:</b> without a reference panel, whether f is shared or redrawn per replicate "
        "makes almost no difference to the LDSC free-intercept h2 estimate across every scale tested. "
        "<b>With a reference panel</b>, however, an f-mismatched reference panel produces estimates that "
        "increasingly diverge from the shared-f case as F_ST grows - and this divergence gets "
        "<b>substantially worse</b> as sample size and SNP count scale up, not better (see Table 1 for "
        "n=500-5000; F_ST=0.3/0.5 were not re-run at n=10000, but the same widening pattern is already "
        "visible within the 0-0.1 range shown in the n=10000 figures below).",
        body))
    story.append(Paragraph(
        "This is a different failure mode than was seen in the project's earlier AR1-based "
        "reference-panel work: there, an f-mismatched reference panel caused the free-intercept "
        "estimator to <i>collapse toward zero</i>. Under real LD structure, the same mismatch instead "
        "causes the estimator to <i>inflate</i>, and the inflation grows with scale.",
        body))
    story.append(PageBreak())

    # --- Methods ---
    story.append(Paragraph("Methods", h2))
    story.append(Paragraph(
        "Genotypes were simulated with LD structure taken from real 1000 Genomes EUR chr1 (106,136 SNPs) "
        "and chr22 (18,613 SNPs) PLINK files, using the data_generation class's realistic=True path: "
        "synthetic individuals are drawn as random linear combinations of the real genotype matrix's "
        "columns, preserving its LD structure while allowing arbitrary sample size. A genetic-drift/"
        "stratification vector f (parameterized by F_ST) splits the sample into two subpopulations, "
        "exactly as in the project's AR1 pipeline.",
        body))
    story.append(Paragraph("Four scenarios were run at every scale and chromosome (chr22 through n=5000; chr1 additionally at n=10000):", body))
    for line in [
        "<b>No ref panel, shared f</b> - LD scores computed in-sample; f drawn once per F_ST and reused across all 100 replicates.",
        "<b>No ref panel, independent f</b> - LD scores computed in-sample; f redrawn fresh every replicate.",
        "<b>Ref panel, shared f</b> - LD scores/GWASH statistics taken from a separately-simulated reference panel; the reference panel's f is shared with the main replicates.",
        "<b>Ref panel, independent f</b> - same reference panel LD scores, but the main replicates draw their own f, independent of (mismatched with) the reference panel's.",
    ]:
        story.append(Paragraph("&bull; " + line, ParagraphStyle("Bullet", parent=body, leftIndent=14)))
    story.append(Paragraph(
        "The reported estimator is LDSC with a free intercept (h2_ldsc_reg), averaged over 100 replicates "
        "per F_ST value. The figures below focus on F_ST in [0, 0.1] (F_ST in {0, 0.001, ..., 0.01, 0.02, "
        "..., 0.1}); Table 1 further below still reports F_ST=0.3/0.5 for the scales where that older, "
        "wider grid was run (n=500 through n=5000), for reference against earlier findings - n=10000 only "
        "has the 0-0.1 grid and so is omitted from Table 1. A theoretical reference line "
        f"h2_pop + {SIGMA_S2}/C_f (C_f = mean(f^2/(1+f^2))) is overlaid on each plot.",
        body))
    story.append(Paragraph("Fixed parameters throughout: h2_pop=0.3, sigma_s=0.3, seed=123. m is always 2n.", body))
    story.append(PageBreak())

    # --- Table 1 ---
    story.append(Paragraph("Summary: h2 estimate gap (independent f &minus; shared f) at F_ST=0.5", h2))
    story.append(Paragraph(
        "Mean LDSC free-intercept h2 estimate over 100 replicates at F_ST=0.3 and F_ST=0.5, for shared-f "
        "and independent-f, in both the no-ref-panel and ref-panel settings. n=10000 is not shown here - "
        "that scale was only run on the newer 0-0.1 F_ST grid.",
        body))
    story.append(styled_table(build_wide_gap_table()))
    story.append(PageBreak())

    # --- Table 2 ---
    story.append(Paragraph("Summary: h2 estimate at low F_ST (F_ST = 0, 0.005, 0.01)", h2))
    story.append(Paragraph(
        "Mean LDSC free-intercept h2 estimate over 100 replicates at F_ST=0, F_ST=0.005, and F_ST=0.01, "
        "for shared-f and independent-f, in both the no-ref-panel and ref-panel settings, across every "
        "available scale including the new n=10000/m=20000 chr1 rows at the end.",
        body))
    story.append(styled_table(build_low_fst_table()))
    story.append(PageBreak())

    # --- Figures ---
    SETTING_TITLES = {
        "ref panel": "Reference-Panel Setting (the main finding)",
        "no ref panel": "No-Reference-Panel Setting (for contrast)",
    }

    def figure_section(chrom, chrom_title):
        story.append(Paragraph(f"Figures: {chrom_title}", h2))
        story.append(Paragraph(
            "Each figure shows mean h2 (&plusmn;95% CI) vs. F_ST for all six LDSC estimator "
            "variants, within one (reference-panel setting, f-sharing) combination, with the "
            "theoretical line overlaid (dotted). Color encodes intercept type (red=free "
            "intercept, pink=fixed intercept); linestyle encodes SNP weighting (solid=weighted, "
            "dashed=unweighted); marker shape encodes estimation mode (circle=two-step, "
            "triangle=single-step, square=fixed-intercept, which has no re-estimation step). "
            "Left panel zooms into F_ST 0-0.01; right panel shows F_ST 0-0.1.",
            body))
        for setting_key, setting_title in SETTING_TITLES.items():
            story.append(Paragraph(setting_title, h3))
            for f_sharing_label, dirname in SCENARIO_DIRS[setting_key].items():
                story.append(Paragraph(f"{setting_key}, {f_sharing_label}", h3))
                for n, m in SCALES_BY_CHROM[chrom]:
                    img_buf = make_plot_image(
                        chrom, dirname, n, m,
                        f"Real {chrom} LD, {setting_key}, {f_sharing_label}",
                    )
                    if img_buf is None:
                        continue
                    story.append(Image(img_buf, width=6.8 * inch, height=3.06 * inch))
                    new_tag = " (new)" if n == 10000 else ""
                    story.append(Paragraph(f"{chrom}, n={n}, m={m} &mdash; {setting_key}, {f_sharing_label}{new_tag}", caption))
                    story.append(PageBreak())

    figure_section("chr1", "Chromosome 1")
    figure_section("chr22", "Chromosome 22")

    # --- Discussion ---
    story.append(Paragraph("Discussion", h2))
    for line in [
        "<b>Why no-ref-panel is robust to f-sharing:</b> when LD scores and GWASH moments are recomputed "
        "in-sample every replicate, the estimator always sees a self-consistent sample - whether that "
        "sample's own f happens to be fixed across replicates or redrawn each time is irrelevant to any "
        "single replicate's estimate. This is confirmed across all scales and both chromosomes, including "
        "the new n=10000 chr1 scale.",
        "<b>Why ref-panel independent-f gets worse with scale:</b> a mismatched reference panel supplies "
        "LD scores computed from a <i>different</i> realization of population structure than the main "
        "sample's. As n and m grow, the main sample's own signal becomes more precisely estimated (tighter "
        "confidence intervals, visible in the figures), but it is being regressed against LD scores that "
        "reflect the wrong stratification pattern. The systematic bias this introduces does not shrink "
        "with more data - sample size reduces variance, not the bias from using a structurally mismatched "
        "reference panel - so the gap between shared-f and independent-f widens as the estimates "
        "themselves become more precise. The n=10000/m=20000 chr1 figure shows this same widening "
        "starting within the 0-0.1 range already (shared f saturating while independent f keeps climbing "
        "toward F_ST=0.1), even without re-running the older 0.2-0.5 grid at this scale.",
        "<b>Inflation vs. collapse:</b> this contrasts with the earlier AR1-based reference-panel finding, "
        "where an f-mismatched reference panel caused the free-intercept estimator to collapse toward "
        "zero rather than inflate. The sign of the failure mode appears to depend on the LD structure "
        "used (synthetic AR1 vs. real chromosome LD), which is itself a notable result: it is not safe to "
        "assume the direction of bias from a mismatched reference panel without knowing the underlying LD "
        "structure.",
        "<b>chr1 vs. chr22:</b> the near-identical gap trajectories (within ~0.1-0.15 of each other at "
        "every shared scale) across two chromosomes of very different size (106,136 vs. 18,613 SNPs) and "
        "gene density suggest this is primarily a sample-size effect of the simulation design, not an "
        "artifact of a particular chromosome's LD architecture.",
        "<b>Practical implication:</b> these results reinforce that a reference panel is only useful for "
        "this kind of heritability estimation if its stratification structure genuinely matches the "
        "target sample's - and that mismatches do not wash out with larger studies. If anything, larger "
        "reference-panel-based studies have more to lose from an f-mismatch than smaller ones.",
    ]:
        story.append(Paragraph(line, body))
    story.append(PageBreak())

    # --- Appendix ---
    story.append(Paragraph("Appendix: Data and Reproducibility", h2))
    story.append(Paragraph(
        f"All {total_csvs} simulation result CSVs (4 scenarios &times; chr1 scales n=500-10000 + chr22 "
        "scales n=500-5000, 100 replicates each) are saved under data/python_realistic_* (chr1) and "
        "data/python_realistic_*_chr22 (chr22) in the project repository, and are not regenerated by the "
        "plotting scripts - replotting or re-analysis can be done directly from these files.",
        body))
    story.append(Paragraph("Driver scripts (python/):", body))
    for line in [
        "do_realistic_sim_Sam_noref_sharedf[_chr22].py, do_realistic_sim_Sam_noref_indepf[_chr22].py",
        "do_realistic_sim_Sam_ref_fixedf[_chr22].py, do_realistic_sim_Sam_ref_indepf[_chr22].py",
        "make_ref_panel_realistic[_chr22].py",
        "plot_realistic_sharedf_vs_indepf[_chr22].py &lt;n&gt; (regenerates figures from saved CSVs only)",
        "generate_realistic_sharedf_vs_indepf_report.py (regenerates this report PDF from saved CSVs only)",
    ]:
        story.append(Paragraph("&bull; " + line, ParagraphStyle("Bullet2", parent=body, leftIndent=14)))
    story.append(Paragraph(
        "Each driver script accepts F_ST and an optional sample size n (m = 2n, default n=500) as "
        "command-line arguments. n=10000/m=20000 chr1 noref_indepf results were generated on a separate "
        "Unix machine (to avoid repeated Windows crashes during a multi-day run) and verified numerically "
        "equivalent to this machine's output for an identical combo before being merged in "
        "(see python/compare_platform_outputs.py and python/logs/platform_benchmark/).",
        body))

    doc = SimpleDocTemplate(
        str(OUT_PDF), pagesize=letter,
        leftMargin=0.6 * inch, rightMargin=0.6 * inch, topMargin=0.6 * inch, bottomMargin=0.6 * inch,
    )
    doc.build(story)
    print(f"wrote {OUT_PDF}")


if __name__ == "__main__":
    main()
