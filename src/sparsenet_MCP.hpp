//  Created by Jin Liu on 17/11/2016.
//  Copyright © 2017年 Jin Liu. All rights reserved.
//


#ifndef sparsenet_MCP_hpp
#define sparsenet_MCP_hpp

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <string.h>
using namespace std;
using namespace Rcpp;
using namespace arma;

arma::fvec ssMCP_fixedSolver(arma::fvec betah, arma::fvec s2, arma::fmat R, float lam, float gamma, arma::fvec beta_old);
List ssMCP_pathSolver(arma::fvec betah, arma::fvec s2, arma::fmat R, arma::fvec lamseq, arma::fvec gammaseq);
List BIC_ssMCP2(arma::fvec betah, arma::fvec s2, arma::fmat R, int nstep, float epsilon, int n, arma::fvec gammaseq);

#endif /* sparsenet_MCP_hpp */
