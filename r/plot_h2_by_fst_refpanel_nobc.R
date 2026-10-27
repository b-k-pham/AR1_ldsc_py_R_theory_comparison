library(ggplot2)
library(ggh4x)

wd = getwd()

## --- Build the combined ref-panel comparison table ---
r_raw = read.csv(file.path(wd, "data", "r_results", "Sam_R_code_fast_refpanel_fixedf_nobc_raw_results.csv"))
r_means = aggregate(cbind(res_reg_free_ref, res_reg_fixed_ref) ~ F_ST, data = r_raw, FUN = mean)

fst_order = c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
fst_label_map = c("0"="00", "0.01"="001", "0.02"="002", "0.03"="003", "0.04"="004", "0.05"="005",
                   "0.1"="01", "0.2"="02", "0.3"="03", "0.4"="04", "0.5"="05")

full = data.frame(F_ST = fst_order)
full = merge(full, r_means, by = "F_ST", all.x = TRUE)
full$h2_ldsc_nw_py = NA_real_
full$h2_ldsc_fixed_nw_py = NA_real_
full$h2_theory_py = NA_real_

for (fst in names(fst_label_map)) {
  fst_num = as.numeric(fst)
  lbl = fst_label_map[[fst]]
  py = read.csv(file.path(wd, "data", "python_refpanel_fixedf", sprintf("AR1_n500_m1000_Fst%s_sigmas03_rho03_h2pop03.csv", lbl)))
  full$h2_ldsc_nw_py[full$F_ST == fst_num]       = mean(py$h2_ldsc_nw)
  full$h2_ldsc_fixed_nw_py[full$F_ST == fst_num] = mean(py$h2_ldsc_fixed_nw)
  full$h2_theory_py[full$F_ST == fst_num]        = mean(py$h2_theory)
}

write.csv(full, file.path(wd, "data", "r_results", "compare_refpanel_fixedf_nobc_full.csv"), row.names = FALSE)

## --- Plot (same two-facet style as plot_h2_by_fst.R) ---
full = full[match(fst_order, full$F_ST), ]

facet_zoom  = "F_ST: 0 - 0.05 (zoomed)"
facet_coarse = "F_ST: 0 - 0.5 (by 0.1)"
facet_levels = c(facet_zoom, facet_coarse)

fst_zoom   = c(0, 0.01, 0.02, 0.03, 0.04, 0.05)
fst_coarse = fst_order
fst_coarse_breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)

add_facet = function(d, facet_name) {
  if (nrow(d) == 0) return(data.frame(F_ST=numeric(0), value=numeric(0), series=character(0), facet=character(0)))
  d$facet = facet_name
  d
}

make_long = function(df, value_col, series_name) {
  base = data.frame(F_ST = df$F_ST, value = df[[value_col]], series = series_name)
  base = base[is.finite(base$value), ]
  rbind(
    add_facet(base[base$F_ST %in% fst_zoom, ],   facet_zoom),
    add_facet(base[base$F_ST %in% fst_coarse, ], facet_coarse)
  )
}

series_colors = c(
  "res_reg_free (R, free, ref panel)"             = "#2a78d6",  # slot1 blue
  "res_reg_fixed (R, fixed=1, ref panel)"         = "#eb6834",  # slot2 orange
  "h2_ldsc_nw (python, free, ref panel)"          = "#1baf7a",  # slot3 aqua
  "h2_ldsc_fixed_nw (python, fixed=1, ref panel)" = "#eda100"   # slot4 yellow
)

long = rbind(
  make_long(full, "res_reg_free_ref",     "res_reg_free (R, free, ref panel)"),
  make_long(full, "res_reg_fixed_ref",    "res_reg_fixed (R, fixed=1, ref panel)"),
  make_long(full, "h2_ldsc_nw_py",        "h2_ldsc_nw (python, free, ref panel)"),
  make_long(full, "h2_ldsc_fixed_nw_py",  "h2_ldsc_fixed_nw (python, fixed=1, ref panel)")
)
long$series = factor(long$series, levels = names(series_colors))
long$facet = factor(long$facet, levels = facet_levels)

theory_zoom = make_long(full, "h2_theory_py", "h2_theory")
theory_zoom = theory_zoom[theory_zoom$facet == facet_zoom, ]

theory_coarse_src = data.frame(F_ST = full$F_ST, value = full$h2_theory_py, series = "h2_theory")
theory_coarse_src$value[is.infinite(theory_coarse_src$value)] = 100  # sentinel only for true Inf (F_ST=0 divide-by-zero)
theory_coarse_src = theory_coarse_src[!is.na(theory_coarse_src$value), ]  # drop rows python's ref panel never had (not Inf, just missing)
theory_coarse = transform(theory_coarse_src, facet = facet_coarse)

theory_long = rbind(theory_zoom, theory_coarse)
theory_long$facet = factor(theory_long$facet, levels = facet_levels)

p = ggplot(long, aes(x = F_ST, y = value, color = series, group = series)) +
  geom_line(data = theory_long, aes(x = F_ST, y = value, linetype = series, group = series),
            inherit.aes = FALSE, color = "#898781", linewidth = 0.6) +
  geom_point(data = theory_long, aes(x = F_ST, y = value, group = series),
             inherit.aes = FALSE, color = "#898781", shape = 4, size = 2) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  facet_wrap(~facet, scales = "free", nrow = 1) +
  facetted_pos_scales(
    x = list(
      scale_x_continuous(),
      scale_x_continuous(breaks = fst_coarse_breaks)
    ),
    y = list(
      scale_y_continuous(),
      scale_y_continuous(limits = c(0, 7), oob = scales::squish)
    )
  ) +
  scale_color_manual(values = series_colors, name = NULL) +
  scale_linetype_manual(values = c("h2_theory" = "dashed"), name = NULL) +
  guides(color = guide_legend(order = 1, nrow = 2, byrow = TRUE),
         linetype = guide_legend(order = 2, nrow = 1)) +
  labs(
    x = "F_ST",
    y = "mean h2 estimate",
    title = "Mean h2 by F_ST: free- vs fixed-intercept estimators, R vs Python (reference panel, f shared/fixed, R uncorrected)",
    subtitle = "Gray x markers = h2_theory (python). Mean over 100 sims per F_ST.\nR's reference-panel LD score here has NO bias correction; python's is still bias-corrected.",
    caption = "n = 500, m = 1000, AR1 rho = 0.3, sigma_s = 0.3. f drawn once per F_ST and shared with the reference panel on both sides."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#e1e0d9"),
    panel.spacing = unit(1.5, "lines"),
    strip.text = element_text(face = "bold", size = 10.5),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 12.5),
    plot.subtitle = element_text(color = "#52514e", size = 10),
    plot.caption = element_text(color = "#898781", size = 9, hjust = 0)
  )

ggsave(file.path(wd, "figures", "h2_by_fst_refpanel_fixedf_nobc.png"), p, width = 11, height = 6.5, dpi = 150)
cat("saved h2_by_fst_refpanel_fixedf_nobc.png\n")
