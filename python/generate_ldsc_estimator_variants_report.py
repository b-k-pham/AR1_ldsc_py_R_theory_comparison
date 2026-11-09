"""
Builds report/ldsc_estimator_variants_comparison.pdf from the data/python_realistic_* CSVs:
compares all six LDSC heritability estimator variants saved per replicate (free/fixed
intercept x weighted/unweighted x two-step/single-step, where two-step only applies when the
intercept is free) under real chr1/chr22 LD structure, across every available (chromosome,
scale, scenario).

This PDF previously existed in report/ with no generator script anywhere in the repo. This
script reconstructs it from the underlying data, reusing the companion report's data-loading
constants (python/generate_realistic_sharedf_vs_indepf_report.py) so the two reports stay
consistent about which scales/chromosomes have data. A scale/scenario with no data for a given
(chrom, n) is skipped for figures and shown as "-" in the table, so this stays usable mid-sweep.

Note: chr1 goes through n=10000/m=20000; chr22 has not been run at that scale. Figures are
grouped by chromosome first (all of chr1, then all of chr22), each with its four scenario
subsections, matching the companion report's layout.

Usage: python generate_ldsc_estimator_variants_report.py
"""
from io import BytesIO

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import Image, PageBreak, Paragraph, SimpleDocTemplate, Table, TableStyle

from generate_realistic_sharedf_vs_indepf_report import (
    CHROMS, DATA_DIR, FST_FINE, FST_FULL, H2_POP, REPORT_DIR, SCALES_BY_CHROM, SIGMA_S2, THEORY_LABEL,
    add_end_labels, csv_path,
    LDSC_VARIANT_COLS as ESTIMATOR_COLS,
    LDSC_VARIANT_COLORS as ESTIMATOR_COLORS,
    LDSC_VARIANT_LINESTYLES as ESTIMATOR_LINESTYLES,
    LDSC_VARIANT_MARKERS as ESTIMATOR_MARKERS,
    LDSC_VARIANT_SHORT_LABELS as ESTIMATOR_SHORT_LABELS,
)

OUT_PDF = REPORT_DIR / "ldsc_estimator_variants_comparison.pdf"

# dir name -> (section label, scenario label used in figure captions/titles)
SCENARIOS = {
    "python_realistic_refpanel_fixedf": "Ref panel, shared f",
    "python_realistic_refpanel_indepf": "Ref panel, independent f",
    "python_realistic_noref_sharedf": "No ref panel, shared f",
    "python_realistic_noref_indepf": "No ref panel, independent f",
}

# ESTIMATOR_COLS/COLORS/LINESTYLES/MARKERS are imported above from the companion report module,
# which is the single source of truth for this encoding (both reports show all six LDSC variants
# and must stay visually consistent).


def load_summary(chrom, dirname, n, m, fst_list):
    rows = []
    for fst in fst_list:
        fp = csv_path(chrom, dirname, n, m, fst)
        if not fp.exists():
            continue
        df = pd.read_csv(fp)
        row = {"Fst": fst}
        for col in ESTIMATOR_COLS:
            row[col + "_mean"] = df[col].mean()
            row[col + "_se"] = df[col].std(ddof=1) / np.sqrt(len(df))
        row["h2_theory_mean"] = df["h2_theory"].mean()
        rows.append(row)
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def estimator_mean(chrom, dirname, n, m, fst, col):
    fp = csv_path(chrom, dirname, n, m, fst)
    if not fp.exists():
        return None
    return pd.read_csv(fp)[col].mean()


def count_all_csvs():
    total = 0
    for chrom, scales in SCALES_BY_CHROM.items():
        for n, m in scales:
            for dirname in SCENARIOS:
                d = DATA_DIR / (dirname + CHROMS[chrom]["dir_suffix"])
                if not d.exists():
                    continue
                pattern = f"realisticLD_{'chr22_' if chrom == 'chr22' else ''}n{n}_m{m}_Fst*_sigmas03_h2pop03.csv"
                total += len(list(d.glob(pattern)))
    return total


# ---------- figures ----------

def make_plot_image(chrom, dirname, n, m, title):
    fine = load_summary(chrom, dirname, n, m, FST_FINE)
    full = load_summary(chrom, dirname, n, m, FST_FULL)
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
        for col, label in ESTIMATOR_COLS.items():
            ax.errorbar(
                sub["Fst"], sub[col + "_mean"], yerr=1.96 * sub[col + "_se"],
                marker=ESTIMATOR_MARKERS[col], markersize=5, markeredgewidth=0.6, markeredgecolor="black",
                color=ESTIMATOR_COLORS[col], linestyle=ESTIMATOR_LINESTYLES[col], label=label,
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
        (ESTIMATOR_SHORT_LABELS[col], full[col + "_mean"].iloc[-1], ESTIMATOR_COLORS[col])
        for col in ESTIMATOR_SHORT_LABELS
    ]
    add_end_labels(axes[1], full["Fst"].iloc[-1], end_label_entries)

    # One shared legend below both panels (not duplicated per-panel, not overlapping data).
    # handlelength is widened so dash/dash-dot/dotted patterns are actually legible at this size.
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=4,
        fontsize=7, frameon=False, handlelength=3.0, columnspacing=1.3, labelspacing=0.5,
    )
    fig.suptitle(title)
    fig.tight_layout(rect=[0, 0.16, 0.87, 1])
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=150)
    plt.close(fig)
    buf.seek(0)
    return buf


# ---------- table ----------

def fmt(v, nd=3):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "-"
    return f"{v:.{nd}f}"


def build_fst01_table():
    header = ["Chrom", "n", "Free/Wt/2step", "Fixed/Wt", "Fixed/Unwt", "Free/Unwt/2step", "Free/Wt/1step", "Free/Unwt/1step"]
    rows = [header]
    dirname = "python_realistic_refpanel_indepf"
    for chrom, scales in SCALES_BY_CHROM.items():
        for n, m in scales:
            vals = [estimator_mean(chrom, dirname, n, m, 0.1, col) for col in ESTIMATOR_COLS]
            rows.append([chrom, str(n)] + [fmt(v) for v in vals])
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
    bullet = ParagraphStyle("Bullet", parent=body, leftIndent=14)
    caption = ParagraphStyle("Caption", parent=styles["Normal"], fontSize=8.5, textColor=colors.HexColor("#666666"), spaceBefore=4)

    story = []
    total_csvs = count_all_csvs()

    # --- Title / Executive Summary ---
    story.append(Paragraph("LDSC Estimator Variant Comparison", h1))
    story.append(Paragraph(
        "All six LDSC variants (free/fixed intercept, weighted/unweighted, and two-step/single-step "
        "where the intercept is free) compared together, under real chr1/chr22 LD structure", sub))
    story.append(Paragraph("Executive Summary", h2))
    story.append(Paragraph(
        "Every simulation in this project's realistic-LD pipeline saves six LDSC heritability estimator "
        "variants per replicate, varying two binary choices: intercept (free vs. fixed at 1) and SNP "
        "weighting (the standard LDSC weights vs. uniform weights). A third choice, the two-step "
        "procedure (which first estimates the intercept excluding high-chi-square SNPs, then "
        "re-estimates the slope with that intercept held fixed) vs. a single joint regression pass with "
        "both parameters free, only applies when the intercept is free - LDSC's own implementation does "
        "not permit combining two-step estimation with a fixed intercept, since there is no intercept "
        "left to estimate in a first step. So the six variants are really: two free-intercept/two-step "
        "variants (weighted and unweighted), two free-intercept/single-step variants (weighted and "
        "unweighted), and two fixed-intercept variants (weighted and unweighted, each necessarily a "
        "single joint regression pass). This report compares all six together, across every scale "
        f"(n=500 through n=10000 on chr1, n=500 through n=5000 on chr22 - chr22 has not yet been run at "
        f"n=10000), chromosome (chr1, chr22), and the project's four core scenarios (reference panel "
        f"present or not, stratification vector f shared with the reference panel or drawn "
        f"independently), for a total of {total_csvs} simulated datasets (100 replicates each).",
        body))
    story.append(Paragraph(
        "<b>Headline finding:</b> in three of the four scenarios (no reference panel with either "
        "f-sharing choice, and a reference panel with <i>shared</i> f), all six variants stay within "
        "roughly a factor of 2 of each other and of h2_pop across the F_ST range shown here (0-0.1; "
        "per reviewer feedback, wider F_ST values are no longer the focus). But in the pathological "
        "scenario identified in the companion report - a reference panel whose stratification vector f "
        "is <i>mismatched</i> with the main sample - the six variants <b>split into two completely "
        "different failure modes</b>, and the split is not between two-step and single-step. The two "
        "variants that <b>collapse</b> as F_ST and scale grow (down to ~0.02-0.03 by F_ST=0.1 at n=5000, "
        "roughly 7-10% of h2_pop) are specifically the two with a <i>free intercept estimated in a "
        "single joint pass</i> (no two-step, no fixed intercept). Every other variant - the two "
        "free-intercept/two-step variants and both fixed-intercept variants - instead <b>inflates</b> "
        "with F_ST and scale, reaching up to 1.7 (roughly 5.8&times; h2_pop) by F_ST=0.1 at the largest "
        "scale, with SNP weighting controlling how far it inflates (weighted variants reach up to ~1.7; "
        "their unweighted counterparts stay below ~0.65).",
        body))
    story.append(Paragraph(
        "This means the real driver of collapse-vs-inflation is whether the intercept is "
        "<i>stabilized</i> in some way - either fixed outright at 1, or estimated robustly first via "
        "the two-step procedure - before the slope is fit. Only when the intercept is both free "
        "<i>and</i> fit jointly with the slope in one pass does the regression have an unconstrained "
        "degree of freedom that can absorb the f-mismatch signal, driving the slope (and hence h2) "
        "toward zero. Fixing the intercept and using two-step are two different ways of preventing that "
        "collapse, not two points on the same axis - and SNP weighting is a largely separate, secondary "
        "axis that scales the degree of inflation among the non-collapsing variants.",
        body))
    story.append(PageBreak())

    # --- Methods ---
    story.append(Paragraph("Methods", h2))
    story.append(Paragraph(
        "See the companion report (realistic_LD_sharedf_vs_indepf_summary.pdf) for full simulation "
        "methodology. The six LDSC variants compared here, all computed from the same per-replicate "
        "data (same X, y, LD scores):", body))
    for col, label in ESTIMATOR_COLS.items():
        story.append(Paragraph(f"&bull; <b>{col}</b>: {label.replace('LDSC ', '')}", bullet))
    story.append(Paragraph(
        "All panels use a linear y-axis. The dark dash-dot line is h2_pop=0.3, the true simulated genetic "
        f"heritability. The black dotted line is the theoretical total h2 (h2_pop + {SIGMA_S2}/C_f, C_f = "
        "mean(f^2/(1+f^2))), using that scenario's own per-replicate f draws - the same formula "
        "validated against the project's R implementation in the companion report. Note that in the "
        "pathological reference-panel/independent-f scenario (within the F_ST 0-0.1 range shown here), "
        "the value range spans roughly 0.02 to 1.7; on a linear axis the two free-intercept/single-step "
        "variants (which shrink toward zero as F_ST and scale grow) can appear visually compressed "
        "against the baseline next to the other four variants (which grow well above h2_pop).",
        body))
    story.append(PageBreak())

    # --- Summary table ---
    story.append(Paragraph("Summary: All Six Estimators at F_ST=0.1, Reference Panel + Independent f", h2))
    story.append(Paragraph(
        "Mean h2 estimate (100 replicates) for each of the six LDSC variants at F_ST=0.1, the upper "
        "edge of the F_ST range this report now focuses on (per reviewer feedback), in the scenario "
        "where the divergence is largest: a reference panel whose f does not match the main sample's. "
        "Includes the new chr1 n=10000/m=20000 row - chr22 has not been run at that scale.",
        body))
    story.append(styled_table(build_fst01_table()))
    story.append(PageBreak())

    # --- Figures, grouped by chromosome ---
    def figure_section(chrom, chrom_title):
        story.append(Paragraph(f"Figures: {chrom_title}", h2))
        story.append(Paragraph(
            "Each figure overlays all six LDSC variants (mean h2 &plusmn;95% CI) vs. F_ST, with the "
            "theoretical line overlaid (dotted). Color encodes intercept type (red=free "
            "intercept, pink=fixed intercept); linestyle encodes SNP weighting (solid=weighted, "
            "dashed=unweighted); marker shape encodes estimation mode (circle=two-step, "
            "triangle=single-step, square=fixed-intercept, which has no re-estimation step). "
            "Left panel zooms into F_ST 0-0.01; right panel shows F_ST 0-0.1.", body))
        for dirname, scenario_label in SCENARIOS.items():
            story.append(Paragraph(scenario_label, h3))
            for n, m in SCALES_BY_CHROM[chrom]:
                img_buf = make_plot_image(chrom, dirname, n, m, f"{chrom}, {scenario_label}")
                if img_buf is None:
                    continue
                story.append(Image(img_buf, width=6.8 * inch, height=3.06 * inch))
                new_tag = " (new)" if n == 10000 else ""
                story.append(Paragraph(f"{chrom}, n={n}, m={m} &mdash; {scenario_label}{new_tag}", caption))
                story.append(PageBreak())

    figure_section("chr1", "Chromosome 1")
    figure_section("chr22", "Chromosome 22")

    # --- Discussion ---
    story.append(Paragraph("Discussion", h2))
    for line in [
        "<b>Correction to an earlier draft of this report:</b> an initial version of this analysis "
        "mislabeled the two fixed-intercept variants as \"two-step\" and attributed the collapse-vs-"
        "inflation split to the two-step/single-step choice itself. LDSC's two-step procedure is only "
        "defined when the intercept is free - the project's own LDSC implementation "
        "(ldsc_barebones.do_ldsc_regression) raises an error if two-step and a fixed intercept are both "
        "requested, and the fixed-intercept driver calls never pass a twostep argument, so those two "
        "variants actually run as a single joint regression pass with the intercept held at 1. The real "
        "split, confirmed below, is between \"intercept stabilized\" (fixed, or free with two-step) and "
        "\"intercept free and jointly fit in one pass\" - not between two-step and single-step as such.",
        "<b>Why the free-intercept/single-step variants collapse:</b> with the intercept both free and "
        "estimated jointly with the slope in one regression pass, the regression has an unconstrained "
        "degree of freedom that can absorb the reference panel's f-mismatch signal through the "
        "intercept term instead of the slope, driving the fitted slope (and hence h2) toward zero. The "
        "other four variants each remove that degree of freedom one way or another: two-step estimates "
        "the intercept first (on a chi-square-filtered SNP subset) and holds it fixed before fitting the "
        "slope; a literally fixed intercept removes the degree of freedom outright. Either safeguard is "
        "enough to prevent collapse - only the combination of no safeguard at all produces it.",
        "<b>Weighting is a separate, secondary axis:</b> among the four non-collapsing variants, SNP "
        "weighting controls the degree of inflation almost independently of the intercept choice: by "
        "F_ST=0.1 at the largest scale, the two weighted variants (free+two-step, fixed) reach "
        "substantially higher estimates than their unweighted counterparts, regardless of whether the "
        "intercept is free+two-step or fixed. Weighting and intercept-stabilization are two different "
        "knobs, not two readings of the same axis.",
        "<b>Consistent with the well-behaved scenarios:</b> in the three scenarios without an "
        "f-mismatched reference panel, all six variants track each other far more closely (see the "
        "figures for 'No ref panel' and 'Ref panel, shared f') - this split is not a general property of "
        "these six estimators, it is specifically triggered by the f-mismatch failure mode.",
        "<b>Practical implication:</b> if a reference panel's stratification structure is suspected to be "
        "mismatched with the target sample, checking whether a free-intercept, single-step LDSC fit "
        "collapses toward zero relative to the two-step or fixed-intercept fits is a cheap diagnostic - "
        "in a well-specified analysis, they should not disagree this dramatically. Note that a fixed "
        "intercept does not protect against inflation the way it protects against collapse - in this "
        "scenario it produced the single largest estimate of all six variants at every scale tested, "
        "from n=500 (0.96, vs. 0.28 for the next-highest variant) through the new n=10000 row (2.06, vs. "
        "1.77) - though the gap to the weighted two-step variant narrows sharply with scale, from "
        "roughly 3.4&times; at n=500 to roughly 1.2&times; at n=10000, suggesting the two could plausibly "
        "cross over at larger scales still.",
    ]:
        story.append(Paragraph(line, body))

    doc = SimpleDocTemplate(
        str(OUT_PDF), pagesize=letter,
        leftMargin=0.6 * inch, rightMargin=0.6 * inch, topMargin=0.6 * inch, bottomMargin=0.6 * inch,
    )
    doc.build(story)
    print(f"wrote {OUT_PDF}")


if __name__ == "__main__":
    main()
