getcoef2=function(fit){
  #stopi = fit$stopi + 1;
  lamseq=fit$lam#[1:stopi];
  ngamma = length(fit$beta);
  bic = fit$bic;
  df = fit$df;
  ll = fit$ll;
  beta =  vector("list",ngamma);#fit$beta[,1:stopi,];
  for( j in 1:ngamma){
     beta[[j]] = fit$beta[[j]];
  }
  nsnp = dim(beta)[1];

  list(lamseq=lamseq,beta=beta,bic = bic, df=df,ll=ll,nsnp = nsnp)
}
