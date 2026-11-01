## n=2000, m=4000 scale-up of Sam_R_code_fast_refpanel.R ("f independent": the reference
## panel draws its own f, each of the 100 main reps draws a fresh independent f; R's
## reference LD score bias-corrected to match Python). Trimmed to compute only the
## reference-panel estimators (res_reg_free_ref / res_reg_fixed_ref) - the in-sample
## columns from the original script aren't used by the comparison plots and their
## per-replicate crossprod(X) (m x m) would dominate runtime at m=4000.
##
## gen_X_block uses the same recursive AR(1) column generator as
## Sam_R_code_fast_refpanel_fixedf_n2000.R - see that file's header for why.

m=4000
n=2000

h2=0.3

sigm_s=0.3

set.seed(123)

rho=0.3

n1 = n/2

#n_sim=10^3
n_sim=100
args = commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) n_sim = as.integer(args[1])
do.scale=TRUE

F_ST_=c(0,0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5)
n_F_ST=length(F_ST_)

gen_X_block = function(f_vec) {
  X = matrix(0, nrow = n, ncol = m)
  X[,1] = rnorm(n)
  sdw = sqrt(1 - rho^2)
  for (j in 2:m) {
    X[,j] = rho*X[,j-1] + rnorm(n, sd = sdw)
  }
  denom = sqrt(1+f_vec^2)
  X[1:n1,] = sweep(X[1:n1,,drop=FALSE], 2, f_vec, "-")
  X[(n1+1):n,] = sweep(X[(n1+1):n,,drop=FALSE], 2, f_vec, "+")
  X = sweep(X, 2, denom, "/")
  X
}

res_m_sim_ref=rep(NA,n_F_ST); res_s_sim_ref=rep(NA,n_F_ST)
res_m_sim_fixed_ref=rep(NA,n_F_ST); res_s_sim_fixed_ref=rep(NA,n_F_ST)
res_theory=rep(NA,n_F_ST); res_theory_se=rep(NA,n_F_ST)

raw_results = data.frame(F_ST = numeric(n_F_ST*n_sim),
                          rep = integer(n_F_ST*n_sim),
                          res_reg_free_ref = numeric(n_F_ST*n_sim),
                          res_reg_fixed_ref = numeric(n_F_ST*n_sim),
                          theory_bias_term = numeric(n_F_ST*n_sim))
row_idx = 0

for (v in 1:n_F_ST)
{
  F_ST=F_ST_[v]

  ## --- Build the external reference panel for this F_ST (once, held fixed across all
  ## reps), drawing its OWN independent f - NOT shared with the main reps below ---
  f_ref = rnorm(m)*sqrt(F_ST)
  X_ref = scale(gen_X_block(f_ref))
  R_ref = crossprod(X_ref)/(n-1)
  hatell_ref = colSums(R_ref^2)
  hatell_ref_corrected = hatell_ref - (m - hatell_ref)/(n - 2)   # bias correction, matches python's do_ldscore_bias_correct
  ldscore_ref_scaled = n/m * hatell_ref_corrected

  res_reg_free_ref=rep(NA,n_sim); res_reg_fixed_ref=rep(NA,n_sim)
  theory_bias_term_v=rep(NA,n_sim)

  for (r in 1:n_sim)
  {
    f_vec=(rnorm(m)*sqrt(F_ST))  # fresh, independent of f_ref - matches Sam_R_code_fast_randf.R / do_AR1_sim_Sam_ref_indepf.py

    X = gen_X_block(f_vec)

    bet=c(rnorm(m,mean=0,sd=sqrt((h2)/m)))
    xi=c(rep(-sigm_s,n/2),rep(sigm_s,n/2))
    Y=X %*% bet+xi+rnorm(n,mean=0,sd=sqrt(1-h2-sigm_s^2))

    X = scale(X)
    Y = scale(Y)

    u2=(crossprod(X,Y))^2/(n-1)

    ## reference-panel LD score (same u2, different predictor)
    res_reg_free_ref[r]=mean(u2*(ldscore_ref_scaled-mean(ldscore_ref_scaled))) / mean((ldscore_ref_scaled-mean(ldscore_ref_scaled))^2)
    res_reg_fixed_ref[r]=mean(ldscore_ref_scaled*(u2-1)) / mean(ldscore_ref_scaled^2)

    C_f_r = mean(f_vec^2/(1+f_vec^2))
    theory_bias_term_v[r] = sigm_s^2/C_f_r

    row_idx = row_idx + 1
    raw_results$F_ST[row_idx] = F_ST
    raw_results$rep[row_idx] = r
    raw_results$res_reg_free_ref[row_idx] = res_reg_free_ref[r]
    raw_results$res_reg_fixed_ref[row_idx] = res_reg_fixed_ref[r]
    raw_results$theory_bias_term[row_idx] = theory_bias_term_v[r]
  }

  cat(sprintf("F_ST=%s done\n", F_ST))

  res_m_sim_ref[v]=mean(res_reg_free_ref)-h2;     res_s_sim_ref[v]=sd(res_reg_free_ref)/sqrt(n_sim)
  res_m_sim_fixed_ref[v]=mean(res_reg_fixed_ref)-h2; res_s_sim_fixed_ref[v]=sd(res_reg_fixed_ref)/sqrt(n_sim)
  res_theory[v]=mean(theory_bias_term_v);         res_theory_se[v]=sd(theory_bias_term_v)/sqrt(n_sim)
}

summary_results = data.frame(F_ST = F_ST_,
                              res_m_sim_ref = res_m_sim_ref, res_s_sim_ref = res_s_sim_ref,
                              res_m_sim_fixed_ref = res_m_sim_fixed_ref, res_s_sim_fixed_ref = res_s_sim_fixed_ref,
                              res_theory = res_theory, res_theory_se = res_theory_se)

safe_write_csv = function(x, file, ...) {
  for (attempt in 1:5) {
    ok = tryCatch({ write.csv(x, file, ...); TRUE }, error = function(e) FALSE)
    if (ok) return(invisible(NULL))
    Sys.sleep(1)
  }
  write.csv(x, file, ...)
}

safe_write_csv(raw_results, "data/r_results/Sam_R_code_fast_refpanel_n2000_raw_results.csv", row.names = FALSE)
safe_write_csv(summary_results, "data/r_results/Sam_R_code_fast_refpanel_n2000_summary.csv", row.names = FALSE)

print(summary_results)
