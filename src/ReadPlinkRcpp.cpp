//  Created by Jin Liu on 23/04/2018.
//  Copyright © 2018年 Jin Liu. All rights reserved.

#include <RcppArmadillo.h>
#include <armadillo>
#include "readPlink.hpp"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

List ReadPlinkRcpp(std::string stringname, arma::uvec avbIndex, unsigned int bandwidth = 200, float shrinkagefactor = 0.9){
  //index of R vector substract 1 for the C indices
  if (avbIndex.size() > 0){
	avbIndex = avbIndex - 1;
  }
  ObjXCorr obj = ReadDataFromFile(stringname, avbIndex, bandwidth, shrinkagefactor);
  List ret;
  ret["X"] = Rcpp::wrap(obj.X);
  ret["corr"] = Rcpp::wrap(obj.corr);
  return ret;
}



//RcppExport SEXP ReadPlinkRcpp(String stringname, arma::uvec avbIndex, int bandwidth=200){

