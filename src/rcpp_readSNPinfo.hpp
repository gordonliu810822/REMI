//  Created by Jin Liu on 17/11/2016.
//  Copyright © 2017年 Jin Liu. All rights reserved.
//

#ifndef rcpp_readSNPinfo_hpp
#define rcpp_readSNPinfo_hpp

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <string.h>
using namespace std;


void ReadSNPinfo(std::string stringname, IntegerVector A1, IntegerVector A2, CharacterVector rsname,
	IntegerVector chr, IntegerVector bp, NumericVector morgan, int N);
void Read_summarystat(std::string stringname, IntegerVector SA1, IntegerVector SA2, CharacterVector rsname,
	NumericVector betah, NumericVector s2, NumericVector sample_size, int N);
Rcpp::List preprocessn(std::string stringname1, std::string stringname2);


#endif /* rcpp_readSNPinfo_hpp */
