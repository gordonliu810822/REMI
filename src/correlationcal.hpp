#ifndef correlationcal_hpp
#define correlationcal_hpp
//#include <armadillo>
#include <RcppArmadillo.h>
#include <stdio.h>
using namespace arma;
using namespace std;
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
void vec_sub(vec& v0, double m_value, double squaresum);
arma::fmat calCorr(arma::Mat<unsigned>* X0, arma::vec index, arma::uvec availIndex, int bandwidth = 200, float shrinkagefactor = 0.9);

void vec_sub2(vec& v0, double m_value);

arma::fmat calSigma(arma::Mat<unsigned>* X0, arma::vec index, arma::uvec availIndex, int bandwidth = 200, float shrinkagefactor = 0.9);
arma::vec calInnerProd(arma::Mat<unsigned>* X0, arma::vec index, arma::vec y);
//arma::fmat calCorr2(arma::Mat<unsigned>* X0, arma::vec index, arma::uvec availIndex, int bandwidth=200);
//arma::mat calCorrU(arma::umat X, arma::vec index,int badwidth = 200);

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//
#endif /* correlationcal_hpp */
