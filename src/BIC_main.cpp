//  Created by Jin Liu on 23/04/2018.
//  Copyright © 2018年 Jin Liu. All rights reserved.

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdio.h>
#include <bitset>
#include <math.h>
#include "plinkfun.hpp"
#include "correlationcal.hpp"
#include "rcpp_ssLassoSolver2.hpp"
#include "ssLassoSolver_Jiao.hpp"
#include "readPlink.hpp"
#include "rcpp_readSNPinfo.hpp"
//#include "sparsenet_MCP.hpp"
#include<iostream>

#define MAX_LEN 20
#define Ngamma 3

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
Rcpp::List REMI_run(std::string stringname1, std::string stringname2, std::string method, int nstep, arma::fvec gamseq, float eps = 0.05, int bandwidth = 200, float shrinkagefactor = 0.9){
	//stringname1, file names (directory) for summary stat (colnames: rsname, A1, A2, betah, s2, sample_size)
	//stringname2, file names (directory) for bim bed fam files
	//summary stat input format: rsname, A1, A2, betah, s2, sample_size
	//method: "lasso" or "mcp"
	std::string lasso("lasso"), mcp("mcp");
	//cout << method.compare(lasso) << endl;
	//cout << method << endl;
	//cout << lasso << endl;
	clock_t t0 = clock();
	cout << "bandwidth:" << bandwidth << endl;

	//preprocess (includes read snpinfo, summary stat and match reference alleles for qc and change direction of effect sizes
	List pp = preprocessn(stringname1, stringname2); //List::create(betah, s2, idxin,n);
	fvec betah = pp["betah"], s2 = pp["s2"];
	uvec idxin = pp["idxin"];
	fvec sample_size = pp["sample_size"];
	CharacterVector rsname = pp["rsname"];
	uvec A1_r = pp["A1_r"];
	uvec A2_r = pp["A2_r"];

	int n = mean(sample_size);

	//load geno type and calculate R
	//cout << "Size of avbIndex: " << idxin.size() << endl;//compare sum(!is.na(N)) in R: Height_1000Q_match.R
	//cout << "min of avbIndex: " << min(idxin) << endl;
	//cout << "max of avbIndex: " << max(idxin) << endl;
	//printf("Start loading plink files: \n");
	ObjXCorr out = ReadDataFromFile(stringname2, idxin - 1, bandwidth, shrinkagefactor);
	fmat R = out.corr;

	List out2;
	cout << endl;
	cout << "Start fitting REMI: " << endl;
	clock_t t1 = clock();
	if (method.compare(lasso) == 0){
		out2 = rcpp_bicSsLasso(betah, s2, R, nstep, eps, n);
	}
	else if (method.compare(mcp) == 0){
		//out2 = BIC_ssMCP2(betah, s2, R, nstep, eps, n, gamseq);
	}
	cout << "Finished fitting REMI in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
	cout << endl;
	cout << "Total analysis done in " << (clock() - t0)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
	//List out2 = ssMCP_pathSolver(betah, s2, R, lams, gamseq);
	
	//List out2 = rcpp_bicSsMCP(betah, s2, R, nstep, eps, n, stopthresh);
	out2["idxin"] = idxin;
	out2["rsname"] = rsname;
	out2["A1_r"] = A1_r;
	out2["A2_r"] = A2_r;

	return out2;
}



// [[Rcpp::export]]
Rcpp::List REMI_run_Jiao(std::string stringname1, std::string stringname2, int nstep, 
	float eps = 0.05, int bandwidth = 200, float shrinkagefactor = 0.9){//, char* A21, char* A22){
    //stringname1, file names (directory) for summary stat (colnames: rsname, A1, A2, betah, s2, sample_size)
    //std::string lasso("lasso");
    
    clock_t t0 = clock();
    cout << "bandwidth:" << bandwidth << endl;

    List pp = preprocessn(stringname1, stringname2); //List::create(betah, s2, idxin,n);
    fvec ytilde = pp["betah"], s2 = pp["s2"];
    uvec idxin = pp["idxin"];
    fvec sample_size = pp["sample_size"];
    CharacterVector rsname = pp["rsname"];
    uvec A1_r = pp["A1_r"];
    uvec A2_r = pp["A2_r"];
    
    int n = mean(sample_size);
    float sigma2 = mean(s2);
    //ObjXCorr out = ReadDataFromFile(stringname2, idxin - 1, bandwidth, shrinkagefactor);
    //fmat Sigma = out.corr;
    List out = ReadDataFromFile2(stringname2, idxin - 1, bandwidth, shrinkagefactor);
    fmat Sigma = out["Sigma"];
    // Mat<unsigned> X = out["X"];

    List out2;
    cout << endl;
    cout << "Start fitting REMI_Jiao: " << endl;
    clock_t t1 = clock();

	
    //if (method.compare(lasso) == 0){
	out2 = bicSsLasso_Jiao(ytilde, Sigma, nstep, eps, n, sigma2);
    //}
    
    cout << "Finished fitting REMI_Jiao in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
    cout << endl;
    cout << "Total analysis done in " << (clock() - t0)*1.0 / CLOCKS_PER_SEC << " sec." << endl;

    
    out2["idxin"] = idxin;
    out2["rsname"] = rsname;
    out2["A1_r"] = A1_r;
    out2["A2_r"] = A2_r;
    out2["Sigma"] = Sigma;

    return out2;
}

/*
// [[Rcpp::export]]
Rcpp::List test_Sigma(std::string stringname1, std::string stringname2, int bandwidth, float shrinkagefactor = 0.9){
    
    //cout << method.compare(lasso) << endl;
    //cout << method << endl;
    //cout << lasso << endl;
    //clock_t t0 = clock();
    //cout << "bandwidth:" << bandwidth << endl;
    
    //preprocess (includes read snpinfo, summary stat and match reference alleles for qc and change direction of effect sizes
    Rcpp::List pp = preprocessn(stringname1, stringname2); //List::create(betah, s2, idxin,n);
    fvec betah = pp["betah"], s2 = pp["s2"];
    uvec idxin = pp["idxin"];
    fvec sample_size = pp["sample_size"];
    CharacterVector rsname = pp["rsname"];
    uvec A1_r = pp["A1_r"];
    uvec A2_r = pp["A2_r"];
    
    
    //cout << "break 1 ..." << endl;
    ObjXCorr out3 = ReadDataFromFile(stringname2, idxin - 1, bandwidth, shrinkagefactor);
    fmat R = out3.corr;
    List out = ReadDataFromFile2(stringname2, idxin - 1, bandwidth, shrinkagefactor);
    //cout << "break 2 ..." << endl;
    // cout << Sigma(0,201);
    fmat Sigma = out["Sigma"];
    Mat<unsigned> X = out["X"];
    
    List out2 = List::create(Rcpp::Named("R") = R,
                             Rcpp::Named("Sigma") = Sigma,
                             Rcpp::Named("X") = X,
                             Rcpp::Named("idxin") = idxin);
    //List out2;
    return out2;
}
*/
