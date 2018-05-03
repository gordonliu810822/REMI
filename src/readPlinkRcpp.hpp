//
//  plinkfun.hpp
//  PlinkRead
//
//  Created by DaiMingwei on 16/10/6.
//  Copyright © 2016年 daviddai. All rights reserved.
//

#ifndef readPlinkRcpp_hpp
#define readPlinkRcpp_hpp

#include <Rcpp.h>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <string.h>
using namespace std;


Rcpp::List ReadPlinkRcpp(std::string stringname, arma::uvec avbIndex, unsigned int bandwidth = 200);

#endif /* readPlinkRcpp_hpp */
