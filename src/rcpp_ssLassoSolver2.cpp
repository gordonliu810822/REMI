//  Created by Jin Liu on 23/04/2018.
//  Copyright © 2018年 Jin Liu. All rights reserved.


#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
#define Ngamma 3

template <class T> const T& max(const T& a, const T& b) {
	return (a < b) ? b : a;     // or: return comp(a,b)?b:a; for version (2)
}

float shrinkage(float a, float kappa){
	float y = max((float)(0.0), a - kappa) - max((float)(0.0), -a - kappa);
	return y;
}

/*template<typename T>
T mod(T a, int n)
{
	return a - floor(a / n)*n;
}*/

List rcpp_rssMCPSolver(arma::fvec betah, arma::fvec s2, arma::fmat R, float lam, float gamma, arma::fvec beta_old, arma::uvec mm, arma::uvec m, int nin){
	int p = betah.n_elem, MaxIters = 100;
	int bw = (R.n_cols - 1) / 2;

	fvec beta = beta_old;

	float tol = 1e-4;
	int iter_outer = 0;
	int iter_inner = 0;
	int iter = 0;// indmax = 0;
	fvec s = sqrt(s2);
	fvec betah_s2 = betah / s2;


	while (iter < MaxIters){
		//loops on all variable
		float dlx = 0.0;
		//cout << "beta_old 865-th: " << beta_old[864] << endl;
		//int ncore = omp_get_num_procs();
		//int coreNum = 2 * ncore;
		//#pragma omp parallel for if( coreNum > 1) num_threads(coreNum)

		for (int j = 0; j < p; j++){
			float tmp;
			fvec Rj = trans(R.row(j));
			if (j < bw){
				tmp = sum(Rj.subvec(bw - j, R.n_cols - 1) % beta.subvec(0, j + bw) / s.subvec(0, j + bw)) - beta(j) / s(j);
			}
			else if (j >= p - bw - 1){
				tmp = sum(Rj.subvec(0, p - j + bw - 1) % beta.subvec(j - bw, p - 1) / s.subvec(j - bw, p - 1)) - beta(j) / s(j);
			}
			//if (j >= bw && j < p - bw - 1){
			else{
				tmp = sum(Rj % beta.subvec(j - bw, j + bw) / s.subvec(j - bw, j + bw)) - beta(j) / s(j);
			}

			float beta_j = betah_s2(j) - tmp / s(j);
			if (fabs(beta_j) <= lam * gamma){
				beta(j) = shrinkage(beta_j, lam)*s2(j) / (1 - 1 / gamma);
			}
			else{
				beta(j) = beta_j*s2(j);
			}


			float diff_beta = beta_old(j) - beta(j);
			float diffabs = fabs(diff_beta);
			if (diffabs > 0){
				beta_old(j) = beta(j);
				//record j if j is not in active set
				if (mm(j) == 0){
					nin = nin + 1;
					mm(j) = nin;
					m(nin - 1) = j;
				}
				dlx = max(dlx, diffabs);

				//if ( fabs(dlx - diffabs) < 1e-10){
				//	indmax = j;
				//}
				//cout << "Update " << j << "-th variable. " << endl;
			}

			/*if (j == 106){
			cout << j + 1 << "-th tmp value.: " << tmp << endl;
			cout << j + 1 << "-th beta value.: " << beta(j) << endl;
			cout << j + 1 << "-th mm value.: " << mm(j) << endl;
			cout << j + 1 << "-th m value.: " << m(j) << endl;
			cout << j + 1 << "-th nin value.: " << nin << endl;
			}*/
		}
		//if ( iter == 0 ){
		//	cout << iter << "-th outer iteration; " << "max difference: " << dlx << endl;
		//}

		//cout << "max difference: " << dlx << " on the " << indmax+1 <<"-th variable"<< endl;
		//cout << "beta 865-th: " << beta[864] << endl;
		//cout << "nin value.: " << nin << endl;

		iter_outer = iter_outer + 1;
		if (dlx < tol){
			break;
		}

		//loops on the active set
		int iter_inner2 = 0;
		while (iter_inner2 < 1000){
			dlx = 0.0;
			iter_inner = iter_inner + 1;
			iter_inner2 = iter_inner2 + 1;

			//cout << "active: beta_old 865-th: " << beta_old[864] << " (outer iter:"<< iter_outer <<")"<< endl;
			for (int k = 0; k < nin; k++){
				int j = m(k);

				float tmp;
				fvec Rj = trans(R.row(j));
				if (j < bw){
					tmp = sum(Rj.subvec(bw - j, R.n_cols - 1) % beta.subvec(0, j + bw) / s.subvec(0, j + bw)) - beta(j) / s(j);
				}
				else if (j >= p - bw - 1){
					tmp = sum(Rj.subvec(0, p - j + bw - 1) % beta.subvec(j - bw, p - 1) / s.subvec(j - bw, p - 1)) - beta(j) / s(j);
				}
				//if (j >= bw && j < p - bw - 1)
				else{
					tmp = sum(Rj % beta.subvec(j - bw, j + bw) / s.subvec(j - bw, j + bw)) - beta(j) / s(j);
				}

				float beta_j = betah_s2(j) - tmp / s(j);
				if (fabs(beta_j) <= lam * gamma){
					beta(j) = shrinkage(beta_j, lam)*s2(j) / (1 - 1 / gamma);
				}
				else{
					beta(j) = beta_j*s2(j);
				}

				float diff_beta = beta_old(j) - beta(j);
				float diffabs = fabs(diff_beta);
				if (diffabs > 0){
					beta_old(j) = beta(j);
					//record j if j is not in active set
					if (mm(j) == 0){
						nin = nin + 1;
						mm(j) = nin;
						m(nin - 1) = j;
					}
					dlx = max(dlx, diffabs);
				}

			}
			//cout << "active: max difference: " << dlx << " (outer iter:" << iter_outer << ")" << endl;
			//cout << "active: beta 865-th: " << beta[864] << " (outer iter:" << iter_outer << ")" << m(294) << endl;
			if (dlx < tol){
				break;
			}
			//if ( iter_inner2 < 10 ){
			//	cout << iter_inner2 << "-th inner iteration; " << "max difference: " << dlx << endl;
			//}

		}
		iter = iter + 1;

	}
	cout << "MCP solved (outer loops: " << iter_outer << ";inner loops: " << iter_inner << endl;
	List z = List::create(Rcpp::Named("beta") = beta,
		Rcpp::Named("mm") = mm,
		Rcpp::Named("m") = m,
		Rcpp::Named("nin") = nin);
	return z;
}

List rcpp_ssLassoSolver(arma::fvec betah, arma::fvec s2, arma::fmat R, float lam, arma::fvec beta_old, arma::uvec mm, arma::uvec m, int nin){
	int p = betah.n_elem, MaxIters = 100;
	int bw = (R.n_cols - 1) / 2;

	fvec beta = beta_old;

	float tol = 1e-4;
	int iter_outer = 0;
	int iter_inner = 0;
	int iter = 0;// indmax = 0;
	fvec s = sqrt(s2);
	fvec betah_s2 = betah / s2;


	while (iter < MaxIters){
		//loops on all variable
		float dlx = 0.0;
		//cout << "beta_old 865-th: " << beta_old[864] << endl;
		//int ncore = omp_get_num_procs();
		//int coreNum = 2 * ncore;
        //#pragma omp parallel for if( coreNum > 1) num_threads(coreNum)

		for (int j = 0; j < p; j++){
			float tmp;
			fvec Rj = trans(R.row(j));
			if (j < bw){
				tmp = sum(Rj.subvec(bw - j, R.n_cols - 1) % beta.subvec(0, j + bw) / s.subvec(0, j + bw)) - beta(j) / s(j);
			}
			else if (j >= p - bw - 1){
				tmp = sum(Rj.subvec(0, p - j + bw - 1) % beta.subvec(j - bw, p - 1) / s.subvec(j - bw, p - 1)) - beta(j) / s(j);
			}
			//if (j >= bw && j < p - bw - 1){
			else{
				tmp = sum(Rj % beta.subvec(j - bw, j + bw) / s.subvec(j - bw, j + bw)) - beta(j) / s(j);
			}

			float beta_j = betah_s2(j) - tmp / s(j);
			beta(j) = shrinkage(beta_j, lam)*s2(j);

			float diff_beta = beta_old(j) - beta(j);
			float diffabs = fabs(diff_beta);
			if (diffabs > 0){
				beta_old(j) = beta(j);
				//record j if j is not in active set
				if (mm(j) == 0){
					nin = nin + 1;
					mm(j) = nin;
					m(nin-1) = j;
				}
				dlx = max(dlx, diffabs);

				//if ( fabs(dlx - diffabs) < 1e-10){
				//	indmax = j;
				//}
				//cout << "Update " << j << "-th variable. " << endl;
			}

			/*if (j == 106){
				cout << j + 1 << "-th tmp value.: " << tmp << endl;
				cout << j + 1 << "-th beta value.: " << beta(j) << endl;
				cout << j + 1 << "-th mm value.: " << mm(j) << endl;
				cout << j + 1 << "-th m value.: " << m(j) << endl;
				cout << j + 1 << "-th nin value.: " << nin << endl;
			}*/
		}
		//if ( iter == 0 ){
		//	cout << iter << "-th outer iteration; " << "max difference: " << dlx << endl;
		//}
		
		//cout << "max difference: " << dlx << " on the " << indmax+1 <<"-th variable"<< endl;
		//cout << "beta 865-th: " << beta[864] << endl;
		//cout << "nin value.: " << nin << endl;

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

			//cout << "active: beta_old 865-th: " << beta_old[864] << " (outer iter:"<< iter_outer <<")"<< endl;
			for (int k = 0; k < nin; k++){
				int j = m(k);

				float tmp;
				fvec Rj = trans(R.row(j));
				if (j < bw){
					tmp = sum(Rj.subvec(bw - j, R.n_cols - 1) % beta.subvec(0, j + bw) / s.subvec(0, j + bw)) - beta(j) / s(j);
				}
				else if (j >= p - bw - 1){
					tmp = sum(Rj.subvec(0, p - j + bw - 1) % beta.subvec(j - bw, p - 1) / s.subvec(j - bw, p - 1)) - beta(j) / s(j);
				}
				//if (j >= bw && j < p - bw - 1)
				else{
					tmp = sum(Rj % beta.subvec(j - bw, j + bw) / s.subvec(j - bw, j + bw)) - beta(j) / s(j);
				}

				float beta_j = betah_s2(j) - tmp / s(j);
				beta(j) = shrinkage(beta_j, lam)*s2(j);

				float diff_beta = beta_old(j) - beta(j);
				float diffabs = fabs(diff_beta);
				if (diffabs > 0){
					beta_old(j) = beta(j);
					//record j if j is not in active set
					if (mm(j) == 0){
						nin = nin + 1;
						mm(j) = nin;
						m(nin-1) = j;
					}
					dlx = max(dlx, diffabs);
				}

			}
			//cout << "active: max difference: " << dlx << " (outer iter:" << iter_outer << ")" << endl;
			//cout << "active: beta 865-th: " << beta[864] << " (outer iter:" << iter_outer << ")" << m(294) << endl;
			if (dlx < tol){
				break;
			}
			//if ( iter_inner2 < 10 ){
			//	cout << iter_inner2 << "-th inner iteration; " << "max difference: " << dlx << endl;
			//}

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

float rcpp_ssll(arma::fvec betah, arma::fvec beta, arma::fvec s2, arma::fmat R){
	float term1 = sum(beta % betah / s2);
	fvec beta_s = beta / sqrt(s2);

	int p = betah.n_elem, ncol = R.n_cols;
	int bw = (ncol - 1) / 2;
	fmat beta_s_mat = zeros<fmat>(p, ncol);

	//cout << "bw value.: " << bw << endl;
	for (int j = 0; j < p; j++){
		if (j < bw){
			beta_s_mat.submat(j, bw - j, j, ncol - 1) = trans(beta_s.subvec(0, j + bw));
			//cout << j + 1 << "-th p. "  << endl;
		}
		else if (j >= p - bw){
			//cout << j + 1 << "-th p. " << endl;
			beta_s_mat.submat(j, 0, j, p - j + bw - 1) = trans(beta_s.subvec(j - bw, p - 1));
			
		}
		else{
			beta_s_mat.row(j) = trans(beta_s.subvec(j - bw, j + bw));
			//cout << j + 1 << "-th p. " << endl;
		}
	}

	fmat beta_s_mat2 = repmat(beta_s, 1, ncol);
	float term3 = -accu(beta_s_mat2 % beta_s_mat % R) / 2;
	float ll = term1 + term3;
	return ll;
}

List ssl(arma::fvec betah, arma::fvec s2, arma::fmat R, int nstep, float eps, int n){
	int p = betah.n_elem;
	fvec betah_s2 = betah / s2;
	fvec abs_betah_s2 = arma::abs(betah_s2);
	float lmax = abs_betah_s2.max();
	float lmin = lmax*eps;
	float ss = (log(lmax) - log(lmin)) / (nstep - 1);
	
	fvec lams = zeros<fvec>(nstep);

	for (int i = 0; i < nstep; i++){
		lams(i) = exp(log(lmax) - ss*(i));
	}

	fvec df = zeros<fvec>(nstep), ll = zeros<fvec>(nstep);
	uvec mm = zeros<uvec>(p), m = zeros<uvec>(p);
	int nin = 0;
	fmat beta = zeros<fmat>(p, nstep);
	fvec beta_old = zeros<fvec>(p);

	for (int i = 0; i < nstep; i++){
		if(i > 0){
			beta_old = beta.col(i-1);
		}
		List out = rcpp_ssLassoSolver(betah, s2, R, lams(i), beta_old, mm, m, nin); 
		
		fvec tmp = out["beta"];
		beta.col(i) = tmp;
		df(i) = sum(beta.col(i) != 0);
		ll(i) = rcpp_ssll(betah, beta.col(i), s2, R);
	}
	List output = List::create(Rcpp::Named("lam") = lams, 
							   Rcpp::Named("df") = df, 
							   Rcpp::Named("ll") = ll, 
							   Rcpp::Named("beta") = beta);
	return output;
}

// [[Rcpp::export]]
List rcpp_bicSsLasso(arma::fvec betah, arma::fvec s2, arma::fmat R, int nstep, float eps, int n){
	int p = betah.n_elem;
	fvec betah_s2 = betah / s2;
	fvec abs_betah_s2 = arma::abs(betah_s2);
	float lmax = abs_betah_s2.max();
	float lmin = lmax*eps;
	float ss = (log(lmax) - log(lmin)) / (nstep - 1);
	
	fvec lams = zeros<fvec>(nstep);

	for (int i = 0; i < nstep; i++){
		lams(i) = exp(log(lmax) - ss*(i));
	}

	fvec bic = zeros<fvec>(nstep), df = -9 * ones<fvec>(nstep), ll = zeros<fvec>(nstep);
	fvec aic = zeros<fvec>(nstep);
	uvec mm = zeros<uvec>(p), m = zeros<uvec>(p);
	int nin = 0;

	fvec beta_old = zeros<fvec>(p);
	uvec idx;
	//int stopi = 0; 

	sp_fmat beta_sparse(p, nstep);
	clock_t t1 = clock();
	
	for (int i = 0; i < nstep; i++){		
		List out = rcpp_ssLassoSolver(betah, s2, R, lams(i), beta_old, mm, m, nin); //List::create(beta, mm, m, nin);
		
		nin = as<int>(out["nin"]);
		uvec tmp1 = out["mm"];
		mm = tmp1;
		uvec tmp2 = out["m"];
		m = tmp2;

	    fvec tmp = out["beta"];
		beta_old = tmp;

		uvec rowidx = find(abs(tmp) > 1e-10);
		fvec values = tmp(rowidx);
		
		for (int j = 0; j < (int)(rowidx.n_elem); j++){
			beta_sparse(rowidx(j), i) = values(j);
		}
		//cout << rowidx.n_elem << endl;

		ll(i) = rcpp_ssll(betah, tmp, s2, R);
		df(i) = sum(abs(tmp) > 1e-10);
		//double alpha = log(p)/log(n);
		bic(i) = -2 * ll(i) + (log(n) )*df(i); // EBIC: (log(n) + 2*(1-1/2/alpha)*log(p))
		aic(i) = -2 * ll(i) + 2*df(i);

		if (i % 5 == 4 ){
			cout << i + 1 << "-th lambda: df " << df(i) << "; bic " << bic(i) << "; lambda " << lams(i);
			cout << "; Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
		}
		

		//stopi = i;
		//if ((int)(df(i)) > 500){
		//	break;
		//}
	}
	
	/*uvec q1 = find(abs(beta) > 1e-10);
	uvec rowidx = mod(q1, p);
	uvec colidx = q1 / p;

	fvec values = beta(q1);
	umat locations(2, q1.n_elem);
	locations.rows(0, 0) = trans(rowidx);
	locations.rows(1, 1) = trans(colidx);
	sp_fmat beta_sparse(locations, values, p, (stopi + 1));*/

	List output = List::create(Rcpp::Named("bic") = bic,
		Rcpp::Named("aic") = aic,
		Rcpp::Named("lam") = lams,
		Rcpp::Named("df") = df,
		Rcpp::Named("ll") = ll,
		Rcpp::Named("beta_sparse") = beta_sparse);
							   //Rcpp::Named("stopi") = stopi);
	return output;
}


List rcpp_bicSsMCP(arma::fvec betah, arma::fvec s2, arma::fmat R, int nstep, float epsilon, int n){
	int p = betah.n_elem;
	fvec betah_s2 = betah / s2;
	fvec abs_betah_s2 = arma::abs(betah_s2);
	float lmax = abs_betah_s2.max();
	float lmin = lmax*epsilon;
	float ss = (log(lmax) - log(lmin)) / (nstep - 1);

	fvec lams = zeros<fvec>(nstep);

	for (int i = 0; i < nstep; i++){
		lams(i) = exp(log(lmax) - ss*(i));
	}
	float gamseq[Ngamma] = { 100000, 6, 2.7 };

	//float alpha[nstep];

	for (int i = 0; i < nstep; i++){

	}

	fmat bic = zeros<fmat>(nstep, Ngamma), df = -9*ones<fmat>(nstep, Ngamma), ll = zeros<fmat>(nstep, Ngamma);
	
	fcube betac = zeros<fcube>(p, nstep, Ngamma);
	fmat beta = zeros<fmat>(p, nstep);
	fvec beta_old = zeros<fvec>(p);
	//ivec stopi = zeros<ivec>(Ngamma);

	//int ncore = omp_get_num_procs();
	//int coreNum = 2 * ncore;
    //#pragma omp parallel for if( coreNum > 1) num_threads(coreNum)

	for (int j = 0; j < Ngamma; j++){
		cout << j << "-th gamma: " << gamseq[j] << endl;
		beta = betac.slice(j);
		uvec mm = zeros<uvec>(p), m = zeros<uvec>(p);
		int nin = 0;

		for (int i = 0; i < nstep; i++){
			if (i > 0){
				beta_old = beta.col(i - 1);
			}
			List out = rcpp_rssMCPSolver(betah, s2, R, lams(i), gamseq[j], beta_old, mm, m, nin); 

			nin = as<int>(out["nin"]);
			uvec tmp1 = out["mm"];
			mm = tmp1;
			uvec tmp2 = out["m"];
			m = tmp2;

			fvec tmp = out["beta"];
			beta.col(i) = tmp;
			ll(i, j) = rcpp_ssll(betah, beta.col(i), s2, R);
			df(i, j) = sum(abs(tmp) > 1e-10);
			bic(i, j) = -2 * ll(i, j) + log(n)*df(i, j);
			
			cout << i << "-th lam: df " << df(i, j) << "; bic " << bic(i, j) << "; lam " << lams(i) << "; gamma " << gamseq[j] << endl;

			//stopi(j) = i;
			//if ((int)(df(i)) > stopthresh){
			//	break;
			//}
		}
		betac.slice(j) = beta;
	}


	List output = List::create(Rcpp::Named("bic") = bic,
		Rcpp::Named("lam") = lams,
		Rcpp::Named("df") = df,
		Rcpp::Named("ll") = ll,
		Rcpp::Named("beta") = betac);
		//Rcpp::Named("stopi") = stopi);
	return output;
}
