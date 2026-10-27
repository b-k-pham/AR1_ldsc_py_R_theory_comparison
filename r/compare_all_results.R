## Master comparison: R (fixed-f draw vs random-f-per-sim draw) vs Python,
## across all estimators we've looked at, for all 9 F_ST values.

wd = getwd()
h2_pop = 0.3

fixedf_raw = read.csv(file.path(wd, "data", "r_results", "Sam_R_code_fast_raw_results.csv"))
randf_raw  = read.csv(file.path(wd, "data", "r_results", "Sam_R_code_fast_randf_raw_results.csv"))

fixedf_means = aggregate(cbind(res_reg_free, res_reg_fixed) ~ F_ST, data = fixedf_raw, FUN = mean)
names(fixedf_means)[2:3] = c("res_reg_free_fixedf", "res_reg_fixed_fixedf")

randf_means = aggregate(cbind(res_reg_free, res_reg_fixed) ~ F_ST, data = randf_raw, FUN = mean)
names(randf_means)[2:3] = c("res_reg_free_randf", "res_reg_fixed_randf")

fixedf_summary = read.csv(file.path(wd, "data", "r_results", "Sam_R_code_fast_summary.csv"))
fixedf_summary$h2_theory_fixedf = h2_pop + fixedf_summary$res_theory

randf_summary = read.csv(file.path(wd, "data", "r_results", "Sam_R_code_fast_randf_summary.csv"))
randf_summary$h2_theory_randf = h2_pop + randf_summary$res_theory

fst_label_map = c("0"="00", "0.005"="0005", "0.01"="001", "0.05"="005",
                   "0.1"="01", "0.2"="02", "0.3"="03", "0.4"="04", "0.5"="05")

rows = list()
for (fst in fixedf_means$F_ST) {
  lbl = fst_label_map[as.character(fst)]
  fname = file.path(wd, "data", "python_insample", sprintf("AR1_n500_m1000_Fst%s_sigmas03_rho03_h2pop03.csv", lbl))
  py = read.csv(fname)

  rows[[length(rows)+1]] = data.frame(
    F_ST = fst,

    res_reg_free_fixedf  = fixedf_means$res_reg_free_fixedf[fixedf_means$F_ST == fst],
    res_reg_free_randf   = randf_means$res_reg_free_randf[randf_means$F_ST == fst],

    res_reg_fixed_fixedf = fixedf_means$res_reg_fixed_fixedf[fixedf_means$F_ST == fst],
    res_reg_fixed_randf  = randf_means$res_reg_fixed_randf[randf_means$F_ST == fst],

    h2_theory_R_fixedf = fixedf_summary$h2_theory_fixedf[fixedf_summary$F_ST == fst],
    h2_theory_R_randf  = randf_summary$h2_theory_randf[randf_summary$F_ST == fst],
    h2_theory_py       = mean(py$h2_theory),

    h2_ldsc_fixed_nw_py = mean(py$h2_ldsc_fixed_nw),
    h2_ldsc_nw_py       = mean(py$h2_ldsc_nw),
    h2_ldsc_nwts_py     = mean(py$h2_ldsc_nwts),
    h2_gwash_py         = mean(py$h2_gwash)
  )
}

full = do.call(rbind, rows)
full = full[order(full$F_ST), ]

cat("\n======================= R: f drawn fixed-per-F_ST vs random-per-sim =======================\n")
print(full[, c("F_ST","res_reg_free_fixedf","res_reg_free_randf",
                "res_reg_fixed_fixedf","res_reg_fixed_randf")], row.names = FALSE, digits = 5)

cat("\n======================= h2_theory: R (both f designs) vs Python =======================\n")
print(full[, c("F_ST","h2_theory_R_fixedf","h2_theory_R_randf","h2_theory_py")], row.names = FALSE, digits = 6)

cat("\n======================= All estimators vs theory, side by side =======================\n")
print(full[, c("F_ST","res_reg_free_randf","res_reg_fixed_randf","h2_ldsc_fixed_nw_py",
                "h2_ldsc_nw_py","h2_ldsc_nwts_py","h2_gwash_py","h2_theory_py")], row.names = FALSE, digits = 5)

cat("\n======================= Differences vs Python's h2_theory (the common reference) =======================\n")
diffs = data.frame(
  F_ST = full$F_ST,
  res_reg_free_fixedf  = full$res_reg_free_fixedf  - full$h2_theory_py,
  res_reg_free_randf   = full$res_reg_free_randf   - full$h2_theory_py,
  res_reg_fixed_fixedf = full$res_reg_fixed_fixedf - full$h2_theory_py,
  res_reg_fixed_randf  = full$res_reg_fixed_randf  - full$h2_theory_py,
  h2_ldsc_fixed_nw_py  = full$h2_ldsc_fixed_nw_py  - full$h2_theory_py,
  h2_ldsc_nw_py        = full$h2_ldsc_nw_py        - full$h2_theory_py,
  h2_ldsc_nwts_py      = full$h2_ldsc_nwts_py      - full$h2_theory_py,
  h2_gwash_py          = full$h2_gwash_py          - full$h2_theory_py,
  h2_theory_R_fixedf   = full$h2_theory_R_fixedf   - full$h2_theory_py,
  h2_theory_R_randf    = full$h2_theory_R_randf    - full$h2_theory_py
)
print(diffs, row.names = FALSE, digits = 4)

write.csv(full, file.path(wd, "data", "r_results", "compare_all_results_full.csv"), row.names = FALSE)
write.csv(diffs, file.path(wd, "data", "r_results", "compare_all_results_diffs.csv"), row.names = FALSE)
