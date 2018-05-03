//  Created by Jin Liu on 23/04/2018.
//  Copyright © 2018年 Jin Liu. All rights reserved.


#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <armadillo>
#include <iostream>
#include <stdio.h>
#include "read_genotype_lasso.hpp"
#include "lassoSolver.hpp"
#include "plinkfun.hpp"

using namespace std;
using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List lasso_geno(std::string stringname, std::string method, int nstep, arma::fvec gamseq, float eps){
	std::string famfile = stringname;
	famfile += ".fam"; 
	std::string bimfile = stringname;
	bimfile += ".bim";
	clock_t t0 = clock();

	int P = getLineNum(bimfile);
	int N = getLineNum(famfile);

	//mat X_std(N,P);
	arma::mat tmp(N, P);
	arma::mat* X_std = new arma::mat(tmp.memptr(), N, P, false, false);
	Rcpp::List input = ReadGenotype(stringname, X_std, N, P);

	//arma::mat X_std = input["X_std"];
	arma::vec pheno = as<vec>(input["Phenotype"]);
	//int N = (int)(pheno.n_elem);
	//int P = (int)(X_std.n_cols);

	//centering outcome
	vec Y = pheno;  
	double meanY = sum(pheno) / N;
	Y = Y - meanY;

	// calculate inner product of each SNP
	vec x_var(P);;
	for (int j = 0; j < P; j++){
		x_var(j) = sum(X_std->col(j) % X_std->col(j));
	}
	clock_t t1 = clock();
	cout << "Start fitting Lasso: " << endl;
	List out = lassoPath(X_std, Y, x_var, nstep, eps);
	cout << "Finished fitting Lasso in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
	cout << endl;
	cout << "Total analysis done in " << (clock() - t0)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
	return out;
}


// [[Rcpp::export]]
Rcpp::List lasso(arma::mat X, arma::vec Y, arma::vec X_var, int nstep, double eps)
{
	uword P = X.n_cols;
	uword N = X.n_rows;

	arma::mat* XX = new  arma::mat(X.memptr(), N, P, false, false);

	List out = lassoPath(XX, Y, X_var, nstep, eps);
	return out;
}
