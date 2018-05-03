//
//  plinkfun.hpp
//  PlinkRead
//
//  Created by DaiMingwei on 16/10/6.
//  Copyright © 2016年 daviddai. All rights reserved.
//

#ifndef genoLoader_pred_hpp
#define genoLoader_pred_hpp

#include <RcppArmadillo.h>
#include <Rcpp.h>

#include <boost/algorithm/string.hpp>
#include "plinkfun.hpp"

using namespace std;
using namespace Rcpp;
using namespace arma;

void ReadPlinkBimFile(std::string stringname, IntegerVector A1, IntegerVector A2, CharacterVector rsname,
                      IntegerVector chr, IntegerVector bp, NumericVector morgan, int P);
CharacterVector charv_subset(CharacterVector x, uvec idx);
List match_SNPs2(std::string stringname1, CharacterVector rsname_2, uvec A1_2, uvec A2_2);
void ReadPlinkFamFile2(std::string stringname, CharacterVector FID, CharacterVector IID,
                       NumericVector pheno, int nrows, int whCol);
Rcpp::List dataLoader(std::string stringname, Rcpp::CharacterVector rsname, arma::uvec A1_r, arma::uvec A2_r,int whCol);

#endif /* genoLoader_pred_hpp */
