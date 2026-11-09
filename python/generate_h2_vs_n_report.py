"""
Builds report/h2_vs_n_chr1_Fst001_summary.pdf: mean h2 as a function of sample size n (m = 2n),
at a fixed F_ST=0.01 and sigma_s=0.3, on real chr1 LD structure. One section per scenario (no
ref panel / ref panel x shared f / independent f), each comparing all six LDSC estimator
variants. chr22 is not included - it has not been run at n=10000.

This is a scale sweep at one fixed, low F_ST, in contrast to the companion reports
(realistic_LD_sharedf_vs_indepf_summary.pdf, ldsc_estimator_variants_comparison.pdf), which fix
n per figure and sweep F_ST. Reuses python/plot_h2_vs_n.py for the figures and data loading.

Usage: python generate_h2_vs_n_report.py
"""
from io import BytesIO

import matplotlib.pyplot as plt
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import Image, PageBreak, Paragraph, SimpleDocTemplate, Table, TableStyle

from generate_ldsc_estimator_variants_report import ESTIMATOR_COLS, SCENARIOS
from generate_realistic_sharedf_vs_indepf_report import REPORT_DIR
from plot_h2_vs_n import CHROM, FST, H2_POP, NS, build_figure, load_series

OUT_PDF = REPORT_DIR / "h2_vs_n_chr1_Fst001_summary.pdf"


def fmt(v, nd=3):
    return "-" if v is None else f"{v:.{nd}f}"


def build_table(series):
    header = ["n"] + ["Free/Wt/2step", "Fixed/Wt", "Fixed/Unwt", "Free/Unwt/2step", "Free/Wt/1step", "Free/Unwt/1step"]
    rows = [header]
    cols = list(ESTIMATOR_COLS)
    for n in NS:
        row = [str(n)]
        for col in cols:
            if col not in series:
                row.append("-")
                continue
            ns_have, means, _ = series[col]
            row.append(fmt(means[ns_have.index(n)]) if n in ns_have else "-")
        rows.append(row)
    return rows


def styled_table(data, col_widths=None):
    t = Table(data, colWidths=col_widths, repeatRows=1)
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1F3864")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F2F2F2")]),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#CCCCCC")),
        ("ALIGN", (1, 0), (-1, -1), "CENTER"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    return t


def main():
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    styles = getSampleStyleSheet()
    h1 = ParagraphStyle("H1", parent=styles["Heading1"], fontSize=18, spaceAfter=4)
    h2 = ParagraphStyle("H2", parent=styles["Heading2"], fontSize=15, spaceBefore=6, spaceAfter=6)
    h3 = ParagraphStyle("H3", parent=styles["Heading3"], fontSize=12, spaceBefore=4, spaceAfter=4)
    sub = ParagraphStyle("Sub", parent=styles["Normal"], fontSize=9, textColor=colors.HexColor("#444444"), spaceAfter=12)
    body = ParagraphStyle("Body", parent=styles["Normal"], fontSize=9.5, leading=13, spaceAfter=8)
    caption = ParagraphStyle("Caption", parent=styles["Normal"], fontSize=8.5, textColor=colors.HexColor("#666666"), spaceBefore=4)

    all_series = {dirname: load_series(dirname) for dirname in SCENARIOS}
    story = []

    # --- Title / Executive Summary ---
    story.append(Paragraph("h2 Estimate vs. Sample Size at Fixed, Low F_ST, by LDSC Variant", h1))
    story.append(Paragraph(
        f"Real chr1 LD structure, F_ST={FST}, sigma_s=0.3, n=500 through n=10000 (m=2n), all six "
        f"LDSC estimator variants, one section per scenario", sub))
    story.append(Paragraph("Executive Summary", h2))
    story.append(Paragraph(
        f"The companion reports fix sample size per figure and sweep F_ST from 0 to 0.1 "
        f"(realistic_LD_sharedf_vs_indepf_summary.pdf), or fix F_ST=0.1 and compare all six LDSC "
        f"estimator variants (ldsc_estimator_variants_comparison.pdf). This report fixes "
        f"F_ST={FST} (a small but nonzero stratification level) and sigma_s=0.3, and sweeps "
        f"sample size n from 500 to 10000 (m=2n) on real chr1 LD structure, comparing all six "
        f"LDSC estimator variants within each of the four scenarios (reference panel present or "
        f"not, stratification vector f shared with the reference panel or drawn independently).",
        body))
    story.append(Paragraph(
        "<b>Headline finding:</b> the six variants split into three groups with distinct "
        "n-scaling behavior, and this split holds across all four scenarios. The two "
        "<b>fixed-intercept</b> variants climb steadily and roughly monotonically with n in every "
        "scenario, from ~0.33-0.39 at n=500 to ~0.53-0.56 at n=10000 - essentially unaffected by "
        "reference-panel choice or f-sharing at this F_ST. The two <b>free-intercept/two-step</b> "
        "variants (including the project's default, h2_ldsc_reg) stay close to flat near h2_pop "
        "through n=2000 in every scenario, then rise more modestly to ~0.37-0.40 by n=10000. The "
        "two <b>free-intercept/single-step</b> variants are the most scenario-dependent: in three "
        "of the four scenarios they stay near baseline through n=5000 and then spike sharply at "
        "n=10000, reaching ~0.51-0.70 (higher than even the fixed-intercept variants in the "
        "weighted case) - but in the one scenario with <i>both</i> a reference panel and an "
        "independently-drawn f, they instead stay flat near h2_pop (~0.28-0.30) all the way to "
        "n=10000, with no late blow-up at all.",
        body))
    story.append(Paragraph(
        "This is a different angle on the same mechanism the companion report identifies: there, "
        "the ref-panel/independent-f scenario is the one where free-intercept/single-step "
        "estimation <i>collapses toward zero</i> as F_ST grows at fixed n. Here, at a fixed, low "
        "F_ST and growing n instead, that same scenario's free-intercept/single-step estimates "
        "stay anchored near h2_pop while the <i>same estimators</i> inflate sharply with scale in "
        "every other scenario. The intercept-stabilizing mechanism that causes collapse at high "
        "F_ST appears to also suppress this n-driven blow-up at low F_ST - the ref-panel/"
        "independent-f scenario is consistently the outlier for these two estimators, just in "
        "opposite directions depending on which parameter is being varied.",
        body))
    story.append(PageBreak())

    # --- Methods ---
    story.append(Paragraph("Methods", h2))
    story.append(Paragraph(
        "See the companion reports for full simulation methodology and the six LDSC variant "
        "definitions. Each point here is the same F_ST=0.01 dataset already saved for the main "
        "report's low-F_ST table, re-sliced to hold F_ST fixed and vary n instead of holding n "
        "fixed and varying F_ST. No new simulations were run for this report. chr22 is not "
        "included - it has not been run at n=10000.",
        body))
    story.append(PageBreak())

    # --- Per-scenario sections ---
    for dirname, scenario_label in SCENARIOS.items():
        series = all_series[dirname]
        story.append(Paragraph(f"Scenario: {scenario_label}", h2))
        story.append(Paragraph(
            f"Mean h2 estimate over 100 replicates at F_ST={FST}, by sample size n and LDSC "
            f"variant. h2_pop={H2_POP}.", body))
        story.append(styled_table(build_table(series)))
        story.append(Paragraph("", body))

        fig = build_figure(series, f"Real {CHROM} LD, {scenario_label}: h2 vs. n at F_ST={FST}, sigma_s=0.3")
        buf = BytesIO()
        fig.savefig(buf, format="png", dpi=150)
        plt.close(fig)
        buf.seek(0)
        story.append(Image(buf, width=6.8 * inch, height=5.27 * inch))
        story.append(Paragraph(f"chr1, F_ST={FST}, sigma_s=0.3 &mdash; {scenario_label}", caption))
        story.append(PageBreak())

    # --- Discussion ---
    story.append(Paragraph("Discussion", h2))
    for line in [
        "<b>Fixed-intercept variants inflate with scale regardless of scenario:</b> both "
        "fixed-intercept variants roughly double from n=500 to n=10000 in every one of the four "
        "scenarios, including the two with no reference panel at all. Since a fixed intercept "
        "removes the regression's main degree of freedom for absorbing an f-mismatch signal, this "
        "scale-driven inflation cannot be an f-mismatch effect - it is a property of fitting LDSC "
        "with a fixed intercept to a larger sample at this F_ST, independent of reference-panel "
        "structure.",
        "<b>Free-intercept/two-step variants are the most stable:</b> these track closest to "
        "h2_pop at small-to-moderate n across all four scenarios, consistent with the companion "
        "report's finding that the two-step procedure is one of the two ways (along with a fixed "
        "intercept) to stabilize the intercept and avoid the worst scale/F_ST-driven distortions.",
        "<b>Free-intercept/single-step variants reveal the sharpest scenario split:</b> these are "
        "the two variants with no safeguard on the intercept at all. In three scenarios they spike "
        "to the highest h2 estimates of all six variants by n=10000 (roughly 0.51-0.70, i.e. 1.7-"
        "2.3&times; h2_pop); in the ref-panel/independent-f scenario specifically, they instead "
        "stay essentially flat near h2_pop through n=10000. This is the same pair of estimators "
        "the companion report shows collapsing toward zero as F_ST grows at fixed n in that same "
        "scenario - the two findings together suggest the ref-panel/independent-f mismatch "
        "specifically interferes with whatever mechanism drives these two estimators' values "
        "upward with scale elsewhere, rather than simply adding a separate bias on top of it.",
        "<b>Practical implication:</b> a free intercept estimated jointly with the slope in a "
        "single pass is the least scale-robust choice among the six variants tested, by a wide "
        "margin, in every scenario except the one where a mismatched reference panel happens to "
        "suppress that instability. Relying on that apparent stability in the mismatched-"
        "reference-panel case would be a mistake - it is a side effect of a known failure mode "
        "(the collapse documented in the companion report), not evidence that the estimator is "
        "well-behaved there.",
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
