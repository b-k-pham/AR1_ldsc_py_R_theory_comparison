## n=2000, m=4000 scale-up of Sam_R_code_fast_refpanel_fixedf.R ("f fixed, shared with
## reference panel", R's reference LD score bias-corrected to match Python).
##
## gen_X_block here uses the SAME recursive AR(1) column generator Python's gen_X()
## already uses (data_generation.py, the `old=False` branch: X[,1] ~ N(0,1), then
## X[,j] = rho*X[,j-1] + N(0,1-rho^2)) instead of the original eigen-decomposition +
## dense m x m matrix-square-root approach. Both produce the exact same stationary
## AR(1)(rho) covariance structure (Var=1, Corr(i,j)=rho^|i-j|); the recursive form is
## O(n*m) instead of O(n*m^2 + m^3), which is the difference between ~5 hours and
## ~10 minutes at m=4000 with R's reference BLAS. Verified empirically to reproduce the
## same marginal variance and lag-correlations before swapping it in here.

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

res_m_sim=rep(NA,n_F_ST); res_s_sim=rep(NA,n_F_ST)
res_m_sim_fixed=rep(NA,n_F_ST); res_s_sim_fixed=rep(NA,n_F_ST)
res_theory=rep(NA,n_F_ST)

raw_results = data.frame(F_ST = numeric(n_F_ST*n_sim),
                          rep = integer(n_F_ST*n_sim),
                          res_reg_free_ref = numeric(n_F_ST*n_sim),
                          res_reg_fixed_ref = numeric(n_F_ST*n_sim))
row_idx = 0

for (v in 1:n_F_ST)
{
  F_ST=F_ST_[v]

  ## f drawn ONCE per F_ST, shared by the reference panel AND every one of the 100 reps
  f_vec = rnorm(m)*sqrt(F_ST)

  ## reference panel: SAME f_vec, independent individuals, bias-corrected
  X_ref = scale(gen_X_block(f_vec))
  R_ref = crossprod(X_ref)/(n-1)
  hatell_ref = colSums(R_ref^2)
  hatell_ref_corrected = hatell_ref - (m - hatell_ref)/(n - 2)
  ldscore_ref_scaled = n/m * hatell_ref_corrected

  res_reg_free_ref=rep(NA,n_sim); res_reg_fixed_ref=rep(NA,n_sim)

  for (r in 1:n_sim)
  {
    X = gen_X_block(f_vec)  # SAME f_vec as the reference panel, fresh individuals

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

  res_m_sim[v]=mean(res_reg_free_ref)-h2;          res_s_sim[v]=sd(res_reg_free_ref)/sqrt(n_sim)
  res_m_sim_fixed[v]=mean(res_reg_fixed_ref)-h2;   res_s_sim_fixed[v]=sd(res_reg_fixed_ref)/sqrt(n_sim)
  C_f=mean(f_vec^2/(1+f_vec^2))
  res_theory[v]=sigm_s^2/C_f
}

summary_results = data.frame(F_ST = F_ST_, res_m_sim = res_m_sim, res_s_sim = res_s_sim,
                              res_m_sim_fixed = res_m_sim_fixed, res_s_sim_fixed = res_s_sim_fixed,
                              res_theory = res_theory)

safe_write_csv = function(x, file, ...) {
  for (attempt in 1:5) {
    ok = tryCatch({ write.csv(x, file, ...); TRUE }, error = function(e) FALSE)
    if (ok) return(invisible(NULL))
    Sys.sleep(1)
  }
  write.csv(x, file, ...)
}

safe_write_csv(raw_results, "data/r_results/Sam_R_code_fast_refpanel_fixedf_n2000_raw_results.csv", row.names = FALSE)
safe_write_csv(summary_results, "data/r_results/Sam_R_code_fast_refpanel_fixedf_n2000_summary.csv", row.names = FALSE)

print(summary_results)
