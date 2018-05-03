plot.ssLasso=function(out, xvar="lambda",label=FALSE,sign.lambda=1,...){
  xvar=match.arg(xvar)
  #plotCoef(x$beta,lambda=x$lambda,df=x$df,dev=x$dev.ratio,label=label,xvar=xvar,...)
  par(mfrow = c(2,1))
  plotCoef(out$beta,lambda=out$lam,df=out$df,label=label,xvar=xvar,...)

  cvobj=out
  xlab="log(Lambda)"
  if(sign.lambda<0)xlab=paste("-",xlab,sep="")
  plot.args=list(x=sign.lambda*log(cvobj$lam),y=cvobj$bic,ylim=range(min(cvobj$bic),max(cvobj$bic)),
  xlab=xlab,ylab="BIC",type="n")
  new.args=list(...)
  if(length(new.args))plot.args[names(new.args)]=new.args
do.call("plot",plot.args)
#error.bars(sign.lambda*log(cvobj$lam),cvobj$cvup,cvobj$cvlo,width=0.01,col="darkgrey")
  points(sign.lambda*log(cvobj$lam),cvobj$bic,pch=20,col="red")
  axis(side=3,at=sign.lambda*log(cvobj$lam),labels=paste(cvobj$df),tick=FALSE,line=0)
  invisible()
}
