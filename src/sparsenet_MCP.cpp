//  Created by Jin Liu on 23/04/2018.
//  Copyright © 2018年 Jin Liu. All rights reserved.

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rcpp_ssLassoSolver2.hpp"

using namespace Rcpp;
using namespace arma;

template <class T> const T& max2(const T& a, const T& b) {
	return (a < b) ? b : a;     // or: return comp(a,b)?b:a; for version (2)
}

float shrinkage2(float a, float kappa){
	float y = max2((float)(0.0), a - kappa) - max2((float)(0.0), -a - kappa);
	return y;
}

