//
//  readPlink.cpp
//  ReadPlinkGit
//
//  Created by DaiMingwei on 16/11/8.
//  Copyright © 2016年 daviddai. All rights reserved.
//

#ifndef readPlink_hpp
#define readPlink_hpp

//#include <Rcpp.h>
//#include <armadillo>
#include <RcppArmadillo.h>
#include <iostream>
#include <stdio.h>
#include "correlationcal.hpp"
#include "plinkfun.hpp"

using namespace Rcpp;
using namespace std;
using namespace arma;

struct ObjXCorr{
    arma::fmat corr;
    arma::Mat<unsigned> X;
    
};

ObjXCorr ReadDataFromFile(string stringname, arma::uvec avbIndex, int bandwidth, float shrinkagefactor);
ObjXCorr ReadDataFromFile(string stringname, arma::uvec avbIndex);
ObjXCorr ReadDataFromFile(string stringname);
List ReadDataFromFile2(string stringname, arma::uvec avbIndex, int bandwidth, float shrinkagefactor);



#endif
