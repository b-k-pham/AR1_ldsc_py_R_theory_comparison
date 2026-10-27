library(ggplot2)
library(ggh4x)

wd = getwd()
full = read.csv(file.path(wd, "data", "r_results", "compare_all_results_full.csv"))

fst_order = c(0, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
full = full[match(fst_order, full$F_ST), ]

# Two facets: a zoomed-in view of the small/unstable F_ST range, and a coarse
# full-range view at 0.1 increments. F_ST=0 is the shared left edge of both.
facet_zoom  = "F_ST: 0 - 0.05 (zoomed)"
facet_coarse = "F_ST: 0 - 0.5 (by 0.1)"
facet_levels = c(facet_zoom, facet_coarse)

fst_zoom   = c(0, 0.005, 0.01, 0.05)
# The coarse panel plots the FULL range (including 0.005/0.01/0.05) so the real peak around
# F_ST=0.01 is visible for every series, not just theory - a straight 0 -> 0.1 segment was
# hiding that the estimators actually peak much higher in between. Tick labels stay at the
# coarse 0.1-step values only (set via scale_x_continuous breaks below).
fst_coarse = fst_order
fst_coarse_breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)

make_long = function(df, value_col, series_name) {
  base = data.frame(F_ST = df$F_ST, value = df[[value_col]], series = series_name)
  rbind(
    transform(base[base$F_ST %in% fst_zoom, ],   facet = facet_zoom),
    transform(base[base$F_ST %in% fst_coarse, ], facet = facet_coarse)
  )
}

series_colors = c(
  "res_reg_free (R, free)"             = "#2a78d6",  # slot1 blue
  "res_reg_fixed (R, fixed=1)"         = "#eb6834",  # slot2 orange
  "h2_ldsc_nw (python, free)"          = "#1baf7a",  # slot3 aqua
  "h2_ldsc_fixed_nw (python, fixed=1)" = "#eda100"   # slot4 yellow
)

## Python always draws f randomly per sim (unchanged, per instruction) - only R's f-design varies.
make_plot = function(r_free_col, r_fixed_col, f_design_label, out_file) {

  long = rbind(
    make_long(full, r_free_col,            "res_reg_free (R, free intercept)"),
    make_long(full, r_fixed_col,            "res_reg_fixed (R, fixed intercept=1)"),
    make_long(full, "h2_ldsc_nw_py",        "h2_ldsc_nw (python, free intercept)"),
    make_long(full, "h2_ldsc_fixed_nw_py",  "h2_ldsc_fixed_nw (python, fixed=1)")
  )
  long$series = factor(long$series, levels = c(
    "res_reg_free (R, free intercept)",
    "res_reg_fixed (R, fixed intercept=1)",
    "h2_ldsc_nw (python, free intercept)",
    "h2_ldsc_fixed_nw (python, fixed=1)"
  ))
  levels(long$series) = c("res_reg_free (R, free)", "res_reg_fixed (R, fixed=1)",
                           "h2_ldsc_nw (python, free)", "h2_ldsc_fixed_nw (python, fixed=1)")
  long$facet = factor(long$facet, levels = facet_levels)

  ## h2_theory is undefined at F_ST=0 (C_f=0, division by zero -> Inf). For the zoom panel we
  ## just omit that point (line starts at its first real value, F_ST=0.005). For the coarse
  ## panel we want a continuous line instead: give F_ST=0 a large sentinel value so the line
  ## draws, and feed in the REAL computed values at F_ST=0.005/0.01/0.05 too (not just the
  ## 0.1-step points) so the segment from 0 to 0.1 traces the actual theory curve, clipped at
  ## the panel's y=5 limit via oob=squish rather than breaking.
  theory_zoom = make_long(full, "h2_theory_py", "h2_theory")
  theory_zoom = theory_zoom[theory_zoom$facet == facet_zoom, ]
  theory_zoom = theory_zoom[is.finite(theory_zoom$value), ]

  theory_coarse_src = data.frame(F_ST = full$F_ST, value = full$h2_theory_py, series = "h2_theory")
  theory_coarse_src$value[!is.finite(theory_coarse_src$value)] = 100  # sentinel for F_ST=0, clipped by oob=squish
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
        scale_x_continuous(),                               # left panel: default ticks (unchanged)
        scale_x_continuous(breaks = fst_coarse_breaks)       # right panel: ticks at 0, 0.1, ..., 0.5 only
      ),
      y = list(
        scale_y_continuous(),                   # left panel (facet_zoom): auto range
        scale_y_continuous(limits = c(0, 7), oob = scales::squish)  # right panel: fixed 0-7, clip don't drop
      )
    ) +
    scale_color_manual(values = series_colors, name = NULL) +
    scale_linetype_manual(values = c("h2_theory" = "dashed"), name = NULL) +
    guides(color = guide_legend(order = 1, nrow = 2, byrow = TRUE),
           linetype = guide_legend(order = 2, nrow = 1)) +
    labs(
      x = "F_ST",
      y = "mean h2 estimate",
      title = sprintf("Mean h2 by F_ST: free- vs fixed-intercept estimators, R vs Python (R: %s)", f_design_label),
      subtitle = "Gray x markers = h2_theory. Mean over 100 sims per F_ST.\nEach panel has its own axis scale - note the y-axis range differs between panels.",
      caption = "n = 500, m = 1000, AR1 rho = 0.3, sigma_s = 0.3. No reference panel."
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
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "#52514e", size = 10),
      plot.caption = element_text(color = "#898781", size = 9, hjust = 0)
    )

  ggsave(file.path(wd, "figures", out_file), p, width = 11, height = 6.5, dpi = 150)
  cat(sprintf("saved %s\n", out_file))
}

make_plot("res_reg_free_randf",  "res_reg_fixed_randf",  "f random per sim",   "h2_by_fst_randf.png")
make_plot("res_reg_free_fixedf", "res_reg_fixed_fixedf", "f fixed per F_ST",   "h2_by_fst_fixedf.png")
