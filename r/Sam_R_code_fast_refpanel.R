## Variant of Sam_R_code_fast_randf.R that ALSO computes R's h2 estimators using an
## EXTERNAL REFERENCE PANEL for the LD score, instead of recomputing it in-sample from
## each replicate's own X - mirroring how Python's do_AR1_sim_Sam.py always works
## (compute_u2_ldscores uses self.ref_ldscores, a single static LD-score vector loaded
## from a separately-simulated reference dataset and reused across all 100 reps).
##
## For a clean, fully-matched comparison against Python, the reference panel's LD score
## here gets the SAME bias correction Python applies in make_ref_panel.py
## (do_ldscore_bias_correct=True): ldscore_corrected = ldscore - (M-ldscore)/(Nbar-2).
##
## f (for the MAIN per-replicate data) is drawn random-per-sim, matching
## Sam_R_code_fast_randf.R - the design we established agrees with Python's h2_theory.
## The reference panel itself is built once per F_ST from its own independent f/X draw,
## exactly as Python's make_ref_panel.py builds one reference dataset per F_ST and reuses it.

m=1000
n=500

h2=0.3

sigm_s=0.3

set.seed(123)

rho=0.3

n1 = n/2

Sigm = rho^abs(outer(1:(m), 1:(m), "-"))

out.Sigm=eigen(Sigm)

Sigm_half=out.Sigm$vectors %*% diag(out.Sigm$values^{1/2}) %*% t(out.Sigm$vectors)


#n_sim=10^3
n_sim=100
args = commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) n_sim = as.integer(args[1])
do.scale=TRUE

F_ST_=c(0,0.005,0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5)
n_F_ST=length(F_ST_)

gen_X_block = function(F_ST, f_vec) {
  Z = array(rnorm(n*m), dim=c(n,m))
  X = Z %*% Sigm_half
  denom = sqrt(1+f_vec^2)
  X[1:n1,] = sweep(X[1:n1,,drop=FALSE], 2, f_vec, "-")
  X[(n1+1):n,] = sweep(X[(n1+1):n,,drop=FALSE], 2, f_vec, "+")
  X = sweep(X, 2, denom, "/")
  X
}

res_m_sim=rep(NA,n_F_ST); res_s_sim=rep(NA,n_F_ST)
res_m_sim_fixed=rep(NA,n_F_ST); res_s_sim_fixed=rep(NA,n_F_ST)
res_m_sim_ref=rep(NA,n_F_ST); res_s_sim_ref=rep(NA,n_F_ST)
res_m_sim_fixed_ref=rep(NA,n_F_ST); res_s_sim_fixed_ref=rep(NA,n_F_ST)
res_theory=rep(NA,n_F_ST); res_theory_se=rep(NA,n_F_ST)

raw_results = data.frame(F_ST = numeric(n_F_ST*n_sim),
                          rep = integer(n_F_ST*n_sim),
                          res_reg_free = numeric(n_F_ST*n_sim),       # in-sample LD score (as before)
                          res_reg_fixed = numeric(n_F_ST*n_sim),
                          res_reg_free_ref = numeric(n_F_ST*n_sim),   # NEW: reference-panel LD score
                          res_reg_fixed_ref = numeric(n_F_ST*n_sim),
                          theory_bias_term = numeric(n_F_ST*n_sim))
row_idx = 0

for (v in 1:n_F_ST)
{
  F_ST=F_ST_[v]

  ## --- Build the external reference panel for this F_ST (once, held fixed across all reps) ---
  f_ref = rnorm(m)*sqrt(F_ST)
  X_ref = gen_X_block(F_ST, f_ref)
  X_ref = scale(X_ref)
  R_ref = crossprod(X_ref)/(n-1)
  hatell_ref = colSums(R_ref^2)
  hatell_ref_corrected = hatell_ref - (m - hatell_ref)/(n - 2)   # bias correction, matches python's do_ldscore_bias_correct
  ldscore_ref_scaled = n/m * hatell_ref_corrected

  res_reg_free=rep(NA,n_sim); res_reg_fixed=rep(NA,n_sim)
  res_reg_free_ref=rep(NA,n_sim); res_reg_fixed_ref=rep(NA,n_sim)
  theory_bias_term_v=rep(NA,n_sim)

  for (r in 1:n_sim)
  {
    f_vec=(rnorm(m)*sqrt(F_ST))  # random per sim, matches Sam_R_code_fast_randf.R

    X = gen_X_block(F_ST, f_vec)

    bet=c(rnorm(m,mean=0,sd=sqrt((h2)/m)))
    xi=c(rep(-sigm_s,n/2),rep(sigm_s,n/2))
    Y=X %*% bet+xi+rnorm(n,mean=0,sd=sqrt(1-h2-sigm_s^2))

    X = scale(X)
    Y = scale(Y)

    u2=(crossprod(X,Y))^2/(n-1)
    R=crossprod(X)/(n-1)
    hatell=colSums(R^2)
    ldscore_scaled = n/m*hatell

    ## in-sample (unchanged from Sam_R_code_fast_randf.R)
    res_reg_free[r]=mean(u2*(ldscore_scaled-mean(ldscore_scaled))) / mean((ldscore_scaled-mean(ldscore_scaled))^2)
    res_reg_fixed[r]=mean(ldscore_scaled*(u2-1)) / mean(ldscore_scaled^2)

    ## reference-panel LD score (same u2, different predictor)
    res_reg_free_ref[r]=mean(u2*(ldscore_ref_scaled-mean(ldscore_ref_scaled))) / mean((ldscore_ref_scaled-mean(ldscore_ref_scaled))^2)
    res_reg_fixed_ref[r]=mean(ldscore_ref_scaled*(u2-1)) / mean(ldscore_ref_scaled^2)

    C_f_r = mean(f_vec^2/(1+f_vec^2))
    theory_bias_term_v[r] = sigm_s^2/C_f_r

    row_idx = row_idx + 1
    raw_results$F_ST[row_idx] = F_ST
    raw_results$rep[row_idx] = r
    raw_results$res_reg_free[row_idx] = res_reg_free[r]
    raw_results$res_reg_fixed[row_idx] = res_reg_fixed[r]
    raw_results$res_reg_free_ref[row_idx] = res_reg_free_ref[r]
    raw_results$res_reg_fixed_ref[row_idx] = res_reg_fixed_ref[r]
    raw_results$theory_bias_term[row_idx] = theory_bias_term_v[r]
  }

  cat(sprintf("F_ST=%s done\n", F_ST))

  res_m_sim[v]=mean(res_reg_free)-h2;             res_s_sim[v]=sd(res_reg_free)/sqrt(n_sim)
  res_m_sim_fixed[v]=mean(res_reg_fixed)-h2;      res_s_sim_fixed[v]=sd(res_reg_fixed)/sqrt(n_sim)
  res_m_sim_ref[v]=mean(res_reg_free_ref)-h2;     res_s_sim_ref[v]=sd(res_reg_free_ref)/sqrt(n_sim)
  res_m_sim_fixed_ref[v]=mean(res_reg_fixed_ref)-h2; res_s_sim_fixed_ref[v]=sd(res_reg_fixed_ref)/sqrt(n_sim)
  res_theory[v]=mean(theory_bias_term_v);         res_theory_se[v]=sd(theory_bias_term_v)/sqrt(n_sim)
}

summary_results = data.frame(F_ST = F_ST_,
                              res_m_sim = res_m_sim, res_s_sim = res_s_sim,
                              res_m_sim_fixed = res_m_sim_fixed, res_s_sim_fixed = res_s_sim_fixed,
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

safe_write_csv(raw_results, "data/r_results/Sam_R_code_fast_refpanel_raw_results.csv", row.names = FALSE)
safe_write_csv(summary_results, "data/r_results/Sam_R_code_fast_refpanel_summary.csv", row.names = FALSE)

print(summary_results)
