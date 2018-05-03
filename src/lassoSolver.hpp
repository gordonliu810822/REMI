//  Created by Jin Liu on 17/11/2016.
//  Copyright © 2017年 Jin Liu. All rights reserved.
//

#ifndef lassoSolver_hpp
#define lassoSolver_hpp

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <string.h>
using namespace std;
using namespace Rcpp;
using namespace arma;

List lasso_solver(arma::mat* X, vec Y, vec X_var, double lam, vec beta_old, uvec mm, uvec m, int nin);
Rcpp::List lassoPath(arma::mat* X, arma::vec Y, arma::vec X_var, int nstep, double eps);

#endif /* lassoSolver_hpp */
