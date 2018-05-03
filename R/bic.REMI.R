bic.ssLasso <- function(stringname1, stringname2, method = "lasso", gamseq = c(10,6), nstep = 100,
  epsilon = 0.4, bandwidth = 200, shrinkagefactor = 0.9){
  this.call=match.call() 
  #fit = .Call('rssLasso_rsslasso_run', PACKAGE = 'ssLasso', stringname1, stringname2, stringname3, nstep, epsilon, n)
   
  #outlist=getcoef(fit)
  #class(outlist)="ssLasso.c"
  fit = sslasso_run(stringname1, stringname2, method, nstep, gamseq, epsilon, bandwidth, shrinkagefactor)
  if (method == "lasso"){
    fit=getcoef(fit);
  }
  else if ( method == "mcp"){
    fit=getcoef2(fit);
  }
  #class(outlist)="ssLasso.c"
  #outlist
  #outlist
  fit$call=this.call
  class(fit)=c("bic.REMI")
  fit
}

