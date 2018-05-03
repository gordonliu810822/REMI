//  Created by Jin Liu on 23/04/2018.
//  Copyright © 2018年 Jin Liu. All rights reserved.


#include <RcppArmadillo.h>
#include <armadillo>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

template <class T> const T& max3(const T& a, const T& b) {
	return (a < b) ? b : a;     // or: return comp(a,b)?b:a; for version (2)
}

double shrinkage3(double a, double kappa){
	double y = max3(0.0, a - kappa) - max3(0.0, -a - kappa);
	return y;
}


List lasso_solver(arma::mat* X, vec Y, vec X_var, double lam, vec beta_old, uvec mm, uvec m, int nin){
	uword p = X->n_cols;
	int  MaxIters = 10;

	double tol = 1e-04;
	int iter_outer = 0;
	int iter_inner = 0;
	int iter = 0;

	vec beta = beta_old;
	vec Yhat = *X * beta_old;

	while (iter < MaxIters){
		//loops on all variable
		double dlx = 0.0;
		//cout << "beta_old 865-th: " << beta_old[864] << endl;
		//int ncore = omp_get_num_procs();
		//int coreNum = 2 * ncore;
		//#pragma omp parallel for if( coreNum > 1) num_threads(coreNum)

		for (int j = 0; j < (int)(p); j++){
			double tmp;
			vec inner1 = X->col(j) % (Y - Yhat);
			tmp = sum(inner1) + X_var(j)*beta_old(j);
			beta(j) = shrinkage3(tmp, lam);
			beta(j) = beta(j) / X_var(j);

			double diff_beta = beta_old(j) - beta(j);
			double diffabs = fabs(diff_beta);
			if (diffabs > 0){
				Yhat = Yhat - X->col(j) * diff_beta;
				beta_old(j) = beta(j);
				//record j if j is not in active set
				if (mm(j) == 0){
					nin = nin + 1;
					mm(j) = nin;
					m(nin - 1) = j;
				}
				dlx = max3(dlx, diffabs);
			}

		}

		iter_outer = iter_outer + 1;
		if (dlx < tol){
			break;
		}

		//loops on the active set
		int iter_inner2 = 0;
		while (iter_inner2 < 100){
			dlx = 0.0;
			iter_inner = iter_inner + 1;
			iter_inner2 = iter_inner2 + 1;
			for (int k = 0; k < nin; k++){
				int j = m(k);

				double tmp;
				vec inner1 = X->col(j) % (Y - Yhat);
				tmp = sum(inner1) + X_var(j)*beta_old(j);
				beta(j) = shrinkage3(tmp, lam);
				beta(j) = beta(j) / X_var(j);

				double diff_beta = beta_old(j) - beta(j);
				double diffabs = fabs(diff_beta);
				if (diffabs > 0){
					Yhat = Yhat - X->col(j) * diff_beta;
					beta_old(j) = beta(j);
					//record j if j is not in active set
					if (mm(j) == 0){
						nin = nin + 1;
						mm(j) = nin;
						m(nin - 1) = j;
					}
					dlx = max3(dlx, diffabs);
				}

			}
			if (dlx < tol){
				break;
			}
		}
		iter = iter + 1;

	}
	//cout << "Lasso solved (outer loops: " << iter_outer << ";inner loops: " << iter_inner << endl;
	List z = List::create(Rcpp::Named("beta") = beta,
		Rcpp::Named("mm") = mm,
		Rcpp::Named("m") = m,
		Rcpp::Named("nin") = nin,
		Rcpp::Named("iter_outer") = iter_outer,
		Rcpp::Named("iter_inner") = iter_inner);
	return z;
}

Rcpp::List lassoPath(arma::mat* X, arma::vec Y, arma::vec X_var, int nstep, double eps){
	uword p = X->n_cols;
	uword n = X->n_rows;

	vec inner_xy = abs(X->t() * Y);
	double lmax = inner_xy.max();//max(abs(inner_xy));
	double lmin = lmax*eps;
	double ss = (log(lmax) - log(lmin)) / (nstep - 1);

	vec lamseq = zeros<vec>(nstep);

	for (int i = 0; i < nstep; i++){
		lamseq(i) = exp(log(lmax) - ss*(i));
	}

	uvec mm = zeros<uvec>(p), m = zeros<uvec>(p);
	int nin = 0;

	vec beta_old = zeros<vec>(p);
	vec df = -9 *ones<vec>(nstep);
	vec bic = zeros<vec>(nstep);
	vec ll = zeros<vec>(nstep);
	sp_mat beta_sparse(p, nstep);
	clock_t t1 = clock();
	vec yhat(n);
	double rss;

	for (int i = 0; i < nstep; i++){
		List out = lasso_solver(X, Y, X_var, lamseq(i), beta_old, mm, m, nin);
		nin = as<int>(out["nin"]);
		uvec tmp1 = out["mm"];
		mm = tmp1;
		uvec tmp2 = out["m"];
		m = tmp2;
		vec tmp = out["beta"];
		beta_old = tmp;

		uvec rowidx = find(abs(tmp) > 1e-10);
		vec values = tmp(rowidx);

		for (int j = 0; j < (int)(rowidx.n_elem); j++){
			beta_sparse(rowidx(j), i) = values(j);
		}
		df(i) = sum(abs(tmp) > 1e-10);
		yhat = *X * tmp;
		yhat = Y - yhat;
		rss = sum(yhat % yhat);

		ll(i) = rss;
		double alpha = log(p)/log(n);
		bic(i) = log(rss / n) + df(i)* (log(n) + 2*(1-1/2/alpha)*log(p)) / n; // EBIC

		if (i % 5 == 4){
			//cout << i + 1 << "-th lambda: df " << df(i) << "; lambda " << lamseq(i);
			cout << i + 1 << "-th lambda: df " << df(i) << "; bic " << bic(i) << "; lambda " << lamseq(i);
			cout << "; Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
		}

		if ((int)(df(i)) > 500){
			break;
		}
	}

	List output = List::create(Rcpp::Named("beta_sparse") = beta_sparse, Rcpp::Named("df") = df,
		Rcpp::Named("bic") = bic, Rcpp::Named("ll") = ll, Rcpp::Named("lamseq") = lamseq);
	return output;
}
