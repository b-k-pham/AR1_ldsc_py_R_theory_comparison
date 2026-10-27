## Variant of Sam_R_code_fast_refpanel.R where f_vec is drawn ONCE per F_ST and SHARED
## between the reference panel AND all 100 main replicates - matching the Python fix
## (do_AR1_sim_Sam_ref.py now passes ref_f=f_ref into data_generation, and gen_X uses
## self.ref_f when provided instead of drawing its own f - data_generation.py:333-334,377-378).
##
## This is the "f fixed, shared with reference panel" design, as opposed to:
##   - Sam_R_code_fast.R / Sam_R_code.R: f fixed per F_ST, but NO reference panel at all
##   - Sam_R_code_fast_randf.R: f random per sim, no reference panel
##   - Sam_R_code_fast_refpanel.R: f random per sim for main reps, reference panel uses its
##     OWN independently-drawn f (the version that showed free-intercept collapsing to ~0)

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

gen_X_block = function(f_vec) {
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

safe_write_csv(raw_results, "data/r_results/Sam_R_code_fast_refpanel_fixedf_raw_results.csv", row.names = FALSE)
safe_write_csv(summary_results, "data/r_results/Sam_R_code_fast_refpanel_fixedf_summary.csv", row.names = FALSE)

print(summary_results)
