#' @title
#' REMI
#' @description
#' Fit REMI-C using marginal information.
#'
#' @param stringname1 File name w/ extension name for summary statistics with columns in the order of rs name, major allele, minor allele, estimated effect size, estimated standard error squared, sample size.
#' @param stringname2 File name w/o extension name for reference panel data in the format of PLINK bfile.
#' @param nstep Number of steps for selecting tuning parameter. Default is 100.
#' @param epsilon Minimum value of tuning parameter in the ratio to the largest tuning parameter. Default is 0.05.
#' @param bandwidth The size of banded correlation to be calculated from reference panel. Default is 200.
#'@param shrinkagefactor The shrinkage amount to the estimated correlation matrix from reference panel. When shrinkagefactor is 1, there is no shrinkage. Default is 0.9.
#' @return List of model parameters including major allele (A1), solution path (beta).
#' @examples
#‘ REMI_C(summarystat_file,plink_file);
#' @details
#' \code{REMI_C} fits the REMI-C with marginal information. It requires to provide the reference panel data for SNPs in PLINK bfile format. The summary statistics in file \code{stringname1} and the reference panel is in file \code{stringname2}. This function is primarily used for comparison with REMI_R.
#' @export
REMI_C <-  function(stringname1, stringname2, nstep = 100, epsilon = 0.05, bandwidth = 200,
   shrinkagefactor = 0.9){
  nvars=as.integer(length(betah))
  this.call=match.call() 

  #REMI_run <- function(stringname1, stringname2, method = "lasso", nstep, c(6,3),
  #        epsilon , bandwidth, shrinkagefactor);
  #fit = REMI_run(stringname1, stringname2, method = "lasso", nstep, c(6,3),
  #        epsilon , bandwidth, shrinkagefactor);
  fit = REMI_run_Jiao(stringname1, stringname2, nstep, epsilon, bandwidth, shrinkagefactor)
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
