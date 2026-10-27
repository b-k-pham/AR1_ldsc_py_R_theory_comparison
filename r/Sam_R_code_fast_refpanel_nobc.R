## Same as Sam_R_code_fast_refpanel.R, EXCEPT the reference panel's LD score is left
## RAW/uncorrected (no (m-hatell)/(n-2) bias-correction term). f is independent between
## the reference panel and each replicate (matches the original _refpanel.R design, i.e.
## what produced Figure 4 / h2_by_fst_refpanel_indepf.png) - only the bias-correction choice
## changes here.

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

F_ST_=c(0,0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5)
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

res_m_sim_ref=rep(NA,n_F_ST); res_s_sim_ref=rep(NA,n_F_ST)
res_m_sim_fixed_ref=rep(NA,n_F_ST); res_s_sim_fixed_ref=rep(NA,n_F_ST)
res_theory=rep(NA,n_F_ST)

raw_results = data.frame(F_ST = numeric(n_F_ST*n_sim),
                          rep = integer(n_F_ST*n_sim),
                          res_reg_free_ref = numeric(n_F_ST*n_sim),
                          res_reg_fixed_ref = numeric(n_F_ST*n_sim))
row_idx = 0

for (v in 1:n_F_ST)
{
  F_ST=F_ST_[v]

  ## --- Build the external reference panel for this F_ST (once, held fixed across all reps,
  ## its own INDEPENDENT f draw, RAW/uncorrected LD score) ---
  f_ref = rnorm(m)*sqrt(F_ST)
  X_ref = scale(gen_X_block(F_ST, f_ref))
  R_ref = crossprod(X_ref)/(n-1)
  hatell_ref = colSums(R_ref^2)
  ldscore_ref_scaled = n/m * hatell_ref   # <-- no bias correction (the only change from _refpanel.R)

  res_reg_free_ref=rep(NA,n_sim); res_reg_fixed_ref=rep(NA,n_sim)

  for (r in 1:n_sim)
  {
    f_vec=(rnorm(m)*sqrt(F_ST))  # random per sim, independent from f_ref

    X = gen_X_block(F_ST, f_vec)

    bet=c(rnorm(m,mean=0,sd=sqrt((h2)/m)))
    xi=c(rep(-sigm_s,n/2),rep(sigm_s,n/2))
    Y=X %*% bet+xi+rnorm(n,mean=0,sd=sqrt(1-h2-sigm_s^2))

    X = scale(X)
    Y = scale(Y)

    u2=(crossprod(X,Y))^2/(n-1)

    res_reg_free_ref[r]=mean(u2*(ldscore_ref_scaled-mean(ldscore_ref_scaled))) / mean((ldscore_ref_scaled-mean(ldscore_ref_scaled))^2)
    res_reg_fixed_ref[r]=mean(ldscore_ref_scaled*(u2-1)) / mean(ldscore_ref_scaled^2)

    row_idx = row_idx + 1
    raw_results$F_ST[row_idx] = F_ST
    raw_results$rep[row_idx] = r
    raw_results$res_reg_free_ref[row_idx] = res_reg_free_ref[r]
    raw_results$res_reg_fixed_ref[row_idx] = res_reg_fixed_ref[r]
  }

  cat(sprintf("F_ST=%s done\n", F_ST))

  res_m_sim_ref[v]=mean(res_reg_free_ref)-h2;        res_s_sim_ref[v]=sd(res_reg_free_ref)/sqrt(n_sim)
  res_m_sim_fixed_ref[v]=mean(res_reg_fixed_ref)-h2; res_s_sim_fixed_ref[v]=sd(res_reg_fixed_ref)/sqrt(n_sim)
  C_f=mean(f_ref^2/(1+f_ref^2))
  res_theory[v]=sigm_s^2/C_f
}

summary_results = data.frame(F_ST = F_ST_,
                              res_m_sim_ref = res_m_sim_ref, res_s_sim_ref = res_s_sim_ref,
                              res_m_sim_fixed_ref = res_m_sim_fixed_ref, res_s_sim_fixed_ref = res_s_sim_fixed_ref,
                              res_theory = res_theory)

safe_write_csv = function(x, file, ...) {
  for (attempt in 1:5) {
    ok = tryCatch({ write.csv(x, file, ...); TRUE }, error = function(e) FALSE)
    if (ok) return(invisible(NULL))
    Sys.sleep(1)
  }
  write.csv(x, file, ...)
}

safe_write_csv(raw_results, "data/r_results/Sam_R_code_fast_refpanel_nobc_raw_results.csv", row.names = FALSE)
safe_write_csv(summary_results, "data/r_results/Sam_R_code_fast_refpanel_nobc_summary.csv", row.names = FALSE)

print(summary_results)
