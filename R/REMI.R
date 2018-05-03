REMI <-  function(stringname1, stringname2, nstep = 100, epsilon = 0.05, bandwidth = 200,
   shrinkagefactor = 0.9){
  nvars=as.integer(length(betah))
  this.call=match.call() 

  #REMI_run <- function(stringname1, stringname2, method = "lasso", nstep, c(6,3),
  #        epsilon , bandwidth, shrinkagefactor);
  fit = REMI_run(stringname1, stringname2, method = "lasso", nstep, c(6,3),
          epsilon , bandwidth, shrinkagefactor);
  beta = fit$beta; 
  A1 = fit$A1_r;
  rsname = fit$rsname;

  out1 = cbind(fit$df,fit$ll,fit$lam, fit$bic);
  colnames(out1) = c("Df","Loglik","Lambda", "BIC")
  outlist=list(diagnostics = out1, beta = beta, A1 = A1,rsname=rsname)
  outlist$call=this.call
  
  class(outlist)="REMI"
  outlist
}
