//  Created by Jin Liu on 23/04/2018.
//  Copyright © 2018年 Jin Liu. All rights reserved.


#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

template <class T> const T& max3(const T& a, const T& b) {
	return (a < b) ? b : a;     // or: return comp(a,b)?b:a; for version (2)
}

float shrinkage3(float a, float kappa){
	float y = max3((float)(0.0), a - kappa) - max3((float)(0.0), -a - kappa);
	return y;
}


List ssLassoSolver_Jiao(arma::fvec ytilde, arma::fmat Sigma, float lam, arma::fvec beta_old, arma::uvec mm, arma::uvec m, int nin){
	int p = ytilde.n_elem, MaxIters = 100;
	int bw = (Sigma.n_cols - 1) / 2;

	fvec beta = beta_old;

	float tol = 1e-4;
	int iter_outer = 0;
	int iter_inner = 0;
	int iter = 0;// indmax = 0;

	while (iter < MaxIters){
		//loops on all variable
		float dlx = 0.0;
		//cout << "beta_old 865-th: " << beta_old[864] << endl;
		//int ncore = omp_get_num_procs();
		//int coreNum = 2 * ncore;
        //#pragma omp parallel for if( coreNum > 1) num_threads(coreNum)

		for (int j = 0; j < p; j++){
			float tmp;
			fvec Sigj = trans(Sigma.row(j));
			if (j < bw){
				fvec Ones1 = ones<fvec>(Sigma.n_cols - bw + j);
				tmp = sum(Sigj.subvec(bw - j, Sigma.n_cols - 1) % beta.subvec(0, j + bw) / Ones1) - beta(j) * Sigma(j, bw);
			}
			else if (j >= p - bw - 1){
				fvec Ones1 = ones<fvec>(p - j + bw);
				tmp = sum(Sigj.subvec(0, p - j + bw - 1) % beta.subvec(j - bw, p - 1) / Ones1) - beta(j) * Sigma(j, bw);
			}
			else{
				fvec Ones1 = ones<fvec>(2 * bw + 1);
				tmp = sum(Sigj % beta.subvec(j - bw, j + bw) / Ones1) - beta(j)* Sigma(j, bw);
			}

			float beta_j = ytilde(j) - tmp;
			beta(j) = shrinkage3(beta_j, lam)/Sigma(j,bw);

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

			//cout << "active: beta_old 865-th: " << beta_old[864] << " (outer iter:"<< iter_outer <<")"<< endl;
			for (int k = 0; k < nin; k++){
				int j = m(k);

				float tmp;
				fvec Sigj = trans(Sigma.row(j));
				if (j < bw){
					fvec Ones1 = ones<fvec>(Sigma.n_cols - bw + j);
					tmp = sum(Sigj.subvec(bw - j, Sigma.n_cols - 1) % beta.subvec(0, j + bw) / Ones1) - beta(j) * Sigma(j, bw);
				}
				else if (j >= p - bw - 1){
					fvec Ones1 = ones<fvec>(p - j + bw);
					tmp = sum(Sigj.subvec(0, p - j + bw - 1) % beta.subvec(j - bw, p - 1) / Ones1) - beta(j) * Sigma(j, bw);
				}
				else{
					fvec Ones1 = ones<fvec>(2*bw + 1);
					tmp = sum(Sigj % beta.subvec(j - bw, j + bw) / Ones1) - beta(j) * Sigma(j, bw);
				}

				float beta_j = ytilde(j) - tmp;
				beta(j) = shrinkage3(beta_j, lam)/Sigma(j,bw);

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
					dlx = max3(dlx, diffabs);
				}

			}
			if (dlx < tol){
				break;
			}
	
		}
		iter = iter + 1;

	}
	List z = List::create(Rcpp::Named("beta") = beta, 
		                  Rcpp::Named("mm") = mm, 
						  Rcpp::Named("m") = m, 
						  Rcpp::Named("nin") = nin,
						  Rcpp::Named("iter_outer") = iter_outer,
						  Rcpp::Named("iter_inner") = iter_inner);
	return z;

}

float ssll_Jiao(arma::fvec ytilde, arma::fvec beta, arma::fmat Sigma){
	fvec Ones1 = ones<fvec>(ytilde.n_elem);
	float term1 = sum(beta % ytilde / Ones1);
    // fvec beta_s = beta / sqrt(s2);
    
    int p = ytilde.n_elem, ncol = Sigma.n_cols;
    int bw = (ncol - 1) / 2;
    //fmat beta_s_mat = zeros<fmat>(p, ncol);
    fmat beta_mat = zeros<fmat>(p, ncol);
    
    //cout << "bw value.: " << bw << endl;
    for (int j = 0; j < p; j++){
        if (j < bw){
            beta_mat.submat(j, bw - j, j, ncol - 1) = trans(beta.subvec(0, j + bw));
            // beta_s_mat.submat(j, bw - j, j, ncol - 1) = trans(beta_s.subvec(0, j + bw));
            //cout << j + 1 << "-th p. "  << endl;
        }
        else if (j >= p - bw){
            //cout << j + 1 << "-th p. " << endl;
            // beta_s_mat.submat(j, 0, j, p - j + bw - 1) = trans(beta_s.subvec(j - bw, p - 1));
            beta_mat.submat(j, 0, j, p - j + bw - 1) = trans(beta.subvec(j - bw, p - 1));
            
        }
        else{
            beta_mat.row(j) = trans(beta.subvec(j - bw, j + bw));
            //beta_s_mat.row(j) = trans(beta_s.subvec(j - bw, j + bw));
            //cout << j + 1 << "-th p. " << endl;
        }
    }
    
    // beta_s_mat2 = repmat(beta_s, 1, ncol);
    fmat beta_mat2 = repmat(beta, 1, ncol);
	fmat Ones2 = ones<fmat>(Sigma.n_rows, Sigma.n_cols);
	float term3 = -accu(beta_mat2 % beta_mat % Sigma / Ones2) / 2;
    float ll = term1 + term3;
    return ll;
}

// [[Rcpp::export]]
List bicSsLasso_Jiao(arma::fvec ytilde, arma::fmat Sigma, int nstep, float eps, int n, float sigma2){
    int p = ytilde.n_elem;
    fvec abs_ytilde = arma::abs(ytilde);
    float lmax = abs_ytilde.max();
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
        // ssLassoSolver_Jiao(arma::fvec ytilde, arma::fmat Sigma, float lam, arma::fvec beta_old, arma::uvec mm, arma::uvec m, int nin)
        List out = ssLassoSolver_Jiao(ytilde, Sigma, lams(i), beta_old, mm, m, nin); //List::create(beta, mm, m, nin);
        
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
        //ssll_Jiao(arma::fvec ytilde, arma::fvec beta, arma::fmat Sigma)
        ll(i) = ssll_Jiao(ytilde, tmp, Sigma);
        df(i) = sum(abs(tmp) > 1e-10);
		//double alpha = log(p)/log(n);
        bic(i) = -2 * ll(i)*(n/sigma2) + (log(n))*df(i); // EBIC:(log(n) + 2*(1-1/2/alpha)*log(p))
        aic(i) = -2 * ll(i) + 2*df(i);
        
        if (i % 5 == 4 ){
            cout << i + 1 << "-th lambda: df " << df(i) << "; bic " << bic(i) << "; lambda " << lams(i);
            cout << "; Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
        }
        
        
        //stopi = i;
        //if ((int)(df(i)) > 500){
        //    break;
        //}
    }
    
     //uvec q1 = find(abs(beta) > 1e-10);
     //uvec rowidx = mod(q1, p);
     //uvec colidx = q1 / p;
     
     //fvec values = beta(q1);
     //umat locations(2, q1.n_elem);
     //locations.rows(0, 0) = trans(rowidx);
     //locations.rows(1, 1) = trans(colidx);
     //sp_fmat beta_sparse(locations, values, p, (stopi + 1));
    
    List output = List::create(Rcpp::Named("bic") = bic,
                               Rcpp::Named("aic") = aic,
                               Rcpp::Named("lam") = lams,
                               Rcpp::Named("df") = df,
                               Rcpp::Named("ll") = ll,
                               Rcpp::Named("beta_sparse") = beta_sparse);
    //Rcpp::Named("stopi") = stopi);
    return output;
}

