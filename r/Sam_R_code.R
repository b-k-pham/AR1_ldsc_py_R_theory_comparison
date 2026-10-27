m=1000
n=500

h2=0.3

sigm_s=0.3

set.seed(123)

rho=0.3


Sigm = rho^abs(outer(1:(m), 1:(m), "-"))

out.Sigm=eigen(Sigm)



out.Sigm=eigen(Sigm)
#v1=out.Sigm$vectors[,1]
#v1

Sigm_half=out.Sigm$vectors %*% diag(out.Sigm$values^{1/2}) %*% t(out.Sigm$vectors)




#n_sim=10^3
n_sim=100
args = commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) n_sim = as.integer(args[1])
res_reg_free=rep(NA,n_sim)
res_reg_fixed=rep(NA,n_sim)
do.scale=TRUE

F_ST_=c(0,0.005,0.01,0.05,0.1,0.2,0.3,0.4,0.5) #0.05,0.1,0.2,0.3,0.4,0.5
#F_ST_=c(0)
n_F_ST=length(F_ST_)

res_m_sim=rep(NA,n_F_ST)
res_s_sim=rep(NA,n_F_ST)
res_m_sim_fixed=rep(NA,n_F_ST)
res_s_sim_fixed=rep(NA,n_F_ST)

res_theory=rep(NA,n_F_ST)

# raw per-replicate results across ALL Fst values, for saving to disk
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
  
  
  for (r in 1:n_sim)
  {
    print(r)
    
    X=array(dim=c(n,m))
    
    X=array(dim=c(n,m))
    Z=array(rnorm(n*m),dim=c(n,m))
    for (i in 1:n)
      if (i<(n/2 +1)) X[i,]=(Sigm_half %*% Z[i,]-f_vec)/(sqrt(1+f_vec^2)) else X[i,]=(Sigm_half %*% Z[i,]+f_vec)/(sqrt(1+f_vec^2))
    
    bet=c(rnorm(m,mean=0,sd=sqrt((h2)/m)))
    
    xi=c(rep(-sigm_s,n/2),rep(sigm_s,n/2))
    
    Y=X %*% bet+xi+rnorm(n,mean=0,sd=sqrt(1-h2-sigm_s^2))
    
    V_y=var(Y)
    
    if (do.scale)
    {
      
      #X=scale(X,center=FALSE)
      X = scale(X)
      Y=scale(Y)
      #Y=Y/as.numeric(sqrt(V_y))
      
      u2=(t(X) %*% Y)^2/(n-1)
      R=t(X) %*% X/(n-1)
      #u2=(t(X) %*% Y)^2/(n)
      #R=t(X) %*% X/(n)
    } else {
      u2=(t(X) %*% Y)^2/n
      R=t(X) %*% X/n
      
    }
    
    d2=diag(R)
    hatell=apply(R^2,2,sum)
    print(mean(hatell))
    
    res_reg_free[r]=mean(u2*(n/m*hatell-mean(n/m*hatell))) /mean((n/m*hatell-mean(n/m*hatell))^2)

    # fixed-intercept version: regress (u2-1) on x=n/m*hatell through the origin
    # (i.e. intercept constrained to 1), instead of the free-intercept slope above.
    res_reg_fixed[r]=mean((n/m*hatell)*(u2-1)) / mean((n/m*hatell)^2)

    #res_den_free[r]=mean((hatell-mean(hatell))^2)

    row_idx = row_idx + 1
    raw_results$F_ST[row_idx] = F_ST
    raw_results$rep[row_idx] = r
    raw_results$res_reg_free[row_idx] = res_reg_free[r]
    raw_results$res_reg_fixed[row_idx] = res_reg_fixed[r]
    raw_results$mean_hatell[row_idx] = mean(hatell)

  }

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

safe_write_csv(raw_results, "data/r_results/Sam_R_code_raw_results.csv", row.names = FALSE)
safe_write_csv(summary_results, "data/r_results/Sam_R_code_summary.csv", row.names = FALSE)

print(summary_results)