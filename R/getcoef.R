getcoef=function(fit){
  lamseq=fit$lam#[1:stopi];
  bic = fit$bic#[1:stopi];
  df = fit$df#[1:stopi];
  ll = fit$ll#[1:stopi];
  beta = fit$beta#[,1:stopi];
  nsnp = dim(beta)[1];
  
  list(lamseq=lamseq,beta=beta,bic = bic, df=df,ll=ll,nsnp = nsnp)
}
