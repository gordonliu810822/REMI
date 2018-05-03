//  Created by Jin Liu on 17/11/2016.
//  Copyright © 2017年 Jin Liu. All rights reserved.
//


#ifndef rcpp_rssLassoSolver2_hpp
#define rcpp_rssLassoSolver2_hpp

#include <Rcpp.h>
#include <RcppArmadillo.h>
/*#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <string.h>
using namespace std;*/
using namespace Rcpp;
using namespace arma;

//float shrinkage(float a, float kappa);
Rcpp::List rcpp_rssMCPSolver(arma::fvec betah, arma::fvec s2, arma::fmat R, float lam, arma::fvec beta_old, arma::uvec mm, arma::uvec m, int nin);
Rcpp::List rcpp_ssLassoSolver(arma::fvec betah, arma::fvec s2, arma::fmat R, float lam, arma::fvec beta_old, arma::uvec mm, arma::uvec m, int nin);
float rcpp_ssll(arma::fvec betah, arma::fvec beta, arma::fvec s2, arma::fmat R);
Rcpp::List rcpp_bicSsMCP(arma::fvec betah, arma::fvec s2, arma::fmat R, int nstep, float epsilon, int n);
Rcpp::List rcpp_bicSsLasso(arma::fvec betah, arma::fvec s2, arma::fmat R, int nstep, float epsilon, int n);
Rcpp::List ssl(arma::fvec betah, arma::fvec s2, arma::fmat R, int nstep, float epsilon, int n);
//Rcpp::List ssLassoSolver_Jiao(arma::fvec ytilde, arma::fmat Sigma, float lam, arma::fvec beta_old, arma::uvec mm, arma::uvec m, int nin);

#endif /* rcpp_rssLassoSolver2_hpp */
