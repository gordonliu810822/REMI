//  Created by Jin Liu on 17/11/2016.
//  Copyright © 2017年 Jin Liu. All rights reserved.
//

#ifndef ssLassoSolver_Jiao_hpp
#define ssLassoSolver_Jiao_hpp

#include <RcppArmadillo.h>
#include <Rcpp.h>
//#include "rcpp_ssLassoSolver2.hpp"
/*#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <string.h>*/

//using namespace std;
using namespace Rcpp;
using namespace arma;

Rcpp::List ssLassoSolver_Jiao(arma::fvec ytilde, arma::fmat Sigma, float lam, arma::fvec beta_old, arma::uvec mm, arma::uvec m, int nin);
float ssll_Jiao(arma::fvec ytilde, arma::fvec beta, arma::fmat Sigma);
Rcpp::List bicSsLasso_Jiao(arma::fvec ytilde, arma::fmat Sigma, int nstep, float eps, int n, float sigma2);

#endif /* ssLassoSolver_Jiao_hpp */
