m=1000
n=500

h2=0.3

sigm_s=0.3

set.seed(123)

rho=0.3

n1 = n/2  # first half of individuals (population 1), matches original script's n/2 split

Sigm = rho^abs(outer(1:(m), 1:(m), "-"))

out.Sigm=eigen(Sigm)  # (duplicate eigen() call removed - it was computed twice in the original script)

Sigm_half=out.Sigm$vectors %*% diag(out.Sigm$values^{1/2}) %*% t(out.Sigm$vectors)


#n_sim=10^3
n_sim=100
args = commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) n_sim = as.integer(args[1])
do.scale=TRUE

F_ST_=c(0,0.005,0.01,0.05,0.1,0.2,0.3,0.4,0.5)
n_F_ST=length(F_ST_)

res_m_sim=rep(NA,n_F_ST)
res_s_sim=rep(NA,n_F_ST)
res_m_sim_fixed=rep(NA,n_F_ST)
res_s_sim_fixed=rep(NA,n_F_ST)
res_theory=rep(NA,n_F_ST)

# raw per-replicate results across ALL Fst values (not overwritten each outer loop like the original script)
raw_results = data.frame(F_ST = numeric(n_F_ST*n_sim),
                          rep = integer(n_F_ST*n_sim),
                          res_reg_free = numeric(n_F_ST*n_sim),
                          res_reg_fixed = numeric(n_F_ST*n_sim),
                          mean_hatell = numeric(n_F_ST*n_sim))
row_idx = 0

for (v in 1:n_F_ST)
{
  F_ST=F_ST_[v]
  f_vec=(rnorm(m)*sqrt(F_ST))

  res_reg_free=rep(NA,n_sim)
  res_reg_fixed=rep(NA,n_sim)

  for (r in 1:n_sim)
  {
    Z=array(rnorm(n*m),dim=c(n,m))

    # Vectorized replacement for the original per-row loop:
    #   for (i in 1:n) X[i,] = (Sigm_half %*% Z[i,] -/+ f_vec) / sqrt(1+f_vec^2)
    # Sigm_half is symmetric, so Sigm_half %*% Z[i,] (as a column) is the same as
    # Z[i,] %*% Sigm_half (as a row) -> do it for all rows at once with one matrix product.
    X = Z %*% Sigm_half
    denom = sqrt(1+f_vec^2)
    X[1:n1,] = sweep(X[1:n1,,drop=FALSE], 2, f_vec, "-")
    X[(n1+1):n,] = sweep(X[(n1+1):n,,drop=FALSE], 2, f_vec, "+")
    X = sweep(X, 2, denom, "/")

    bet=c(rnorm(m,mean=0,sd=sqrt((h2)/m)))

    xi=c(rep(-sigm_s,n/2),rep(sigm_s,n/2))

    Y=X %*% bet+xi+rnorm(n,mean=0,sd=sqrt(1-h2-sigm_s^2))

    if (do.scale)
    {
      X = scale(X)
      Y = scale(Y)

      u2=(crossprod(X,Y))^2/(n-1)
      R=crossprod(X)/(n-1)
    } else {
      u2=(crossprod(X,Y))^2/n
      R=crossprod(X)/n
    }

    hatell=colSums(R^2)

    ldscore_scaled = n/m*hatell
    res_reg_free[r]=mean(u2*(ldscore_scaled-mean(ldscore_scaled))) / mean((ldscore_scaled-mean(ldscore_scaled))^2)

    # fixed-intercept version: regress (u2-1) on ldscore_scaled through the origin
    # (i.e. intercept constrained to 1), instead of the free-intercept slope above.
    res_reg_fixed[r]=mean(ldscore_scaled*(u2-1)) / mean(ldscore_scaled^2)

    row_idx = row_idx + 1
    raw_results$F_ST[row_idx] = F_ST
    raw_results$rep[row_idx] = r
    raw_results$res_reg_free[row_idx] = res_reg_free[r]
    raw_results$res_reg_fixed[row_idx] = res_reg_fixed[r]
    raw_results$mean_hatell[row_idx] = mean(hatell)
  }

  cat(sprintf("F_ST=%s done\n", F_ST))

  res_m_sim[v]=mean(res_reg_free)-h2
  res_s_sim[v]=sd(res_reg_free)/sqrt(n_sim)
  res_m_sim_fixed[v]=mean(res_reg_fixed)-h2
  res_s_sim_fixed[v]=sd(res_reg_fixed)/sqrt(n_sim)
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
    Sys.sleep(1)  # transient lock (e.g. AV/indexer scanning the just-written file) - retry
  }
  write.csv(x, file, ...)  # final attempt: let it raise if still locked
}

safe_write_csv(raw_results, "data/r_results/Sam_R_code_fast_raw_results.csv", row.names = FALSE)
safe_write_csv(summary_results, "data/r_results/Sam_R_code_fast_summary.csv", row.names = FALSE)

print(summary_results)
