//  Created by Jin Liu on 23/04/2018.
//  Copyright © 2018年 Jin Liu. All rights reserved.

#include <RcppArmadillo.h>
#include <Rcpp.h>

#include <boost/algorithm/string.hpp>
#include "plinkfun.hpp"

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

#define MAX_LEN 20


void ReadPlinkBimFile(std::string stringname, IntegerVector A1, IntegerVector A2, CharacterVector rsname,
	IntegerVector chr, IntegerVector bp, NumericVector morgan, int P){


	FILE *stream;
	/*CharacterVector rsname(P);
	IntegerVector A1(P), A2(P);*/

	int ch, b;
	double mor;
	char s[MAX_LEN + 1], efa, nefa;

	stream = fopen(stringname.c_str(), "r");
	clock_t t1 = clock();

	//int i = 0;
	/* Put in various data. */
	for ( int i = 0; i < P; i++){
		if (i % 100000 == 0 && i != 0){
			cout << i << "-th SNP" << ",";
			cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
		}

		fscanf(stream, "%i %s %lf %i %c %c", &ch, &s[0], &mor, &b, &efa, &nefa);

		chr(i) = ch;
		rsname(i) = s;
		morgan(i) = mor;
		bp(i) = b;
		A1(i) = (int)efa;
		A2(i) = (int)nefa;

	}
}

CharacterVector charv_subset(CharacterVector x, uvec idx){
	CharacterVector v(idx.n_elem);
	for (unsigned int i = 0; i < idx.n_elem; i++){
		v(i) = x(idx(i));
	}
	return v;
}

List match_SNPs2(std::string stringname1, CharacterVector rsname_2, uvec A1_2, uvec A2_2){
	
	// stringname1: prefix for plink file 1; start working on plink file 1
	string bimfile1 = stringname1;
	bimfile1 += ".bim";
	//cout << "break 1: " << bimfile1  << endl;
	int P1 = getLineNum(bimfile1);
	cout << "## Number of SNPs (plink file 1):" << P1 << endl;

	IntegerVector A1_1(P1), A2_1(P1);
	CharacterVector rsname_1(P1);
	IntegerVector chr_1(P1), bp_1(P1);
	NumericVector morgan_1(P1);

	ReadPlinkBimFile(bimfile1, A1_1, A2_1, rsname_1, chr_1, bp_1, morgan_1, P1);
	// cout << rsname_1(0) << ";" << A1_1(0) << ";" << A2_1(0) << endl;

	
	// mathcing panel SNPs in file 1 and file 2 with correction for direction of ref allele
	// rsname in both files in the order of the first file.
	CharacterVector rs_inter = intersect(rsname_1, rsname_2);
	IntegerVector idxin = match(rsname_1, rs_inter); //index for SNPs in file 1
	CharacterVector rsname_4use = rsname_1[Rcpp::is_na(idxin) == false];
	IntegerVector chr_4use = chr_1[Rcpp::is_na(idxin) == false];
	IntegerVector bp_4use = bp_1[Rcpp::is_na(idxin) == false];

	// match snps (rsname_4use; rsname_2: file 2)
	IntegerVector idxin2 = match(rsname_4use, rsname_2);  //index for SNPs in file 2
	// match snps (rsname_4use; rsname_1: file 1)
	IntegerVector idxin1 = match(rsname_4use, rsname_1);  //index for SNPs in file 1

	/*IntegerVector chr_4use_tmp = chr_2[idxin2];
	IntegerVector bp_4use_tmp = bp_2[idxin2];

	vec idxtmp = as<vec>(chr_4use) - as<vec>(chr_4use_tmp);
	cout << "check the quality: " << sum(idxtmp) << endl;

	cout << "Size of matched SNPs: " << rsname_4use.size() << endl;*/

	// convert ascii letter to uvec and work on overlapped SNPs
	uvec idx = as<uvec>(idxin1) -1;
	uvec tmp1 = as<uvec>(A1_1);
	uvec A1_1_ = tmp1.elem(idx);
	tmp1 = as<uvec>(A2_1);
	uvec A2_1_ = tmp1.elem(idx);

	idx = as<uvec>(idxin2) -1;
	uvec tmp2 = A1_2;
	uvec A1_2_ = tmp2.elem(idx);
	tmp2 = A2_2;
	uvec A2_2_ = tmp2.elem(idx);

	//compare A1_1_ A1_2_, A2_1_ A2_2_
	//A1_1_: replace T with A,replace G with C
	idx = find(A1_1_ == 84);
	uvec idxrepl1(idx.n_elem);
	idxrepl1.fill(65);
	A1_1_.elem(idx) = idxrepl1;

	idx = find(A1_1_ == 71);
	uvec idxrepl2(idx.n_elem);
	idxrepl2.fill(67);
	A1_1_.elem(idx) = idxrepl2;

	//A1_2_: replace T with A,replace G with C
	idx = find(A1_2_ == 84);
	uvec idxrepl3(idx.n_elem);
	idxrepl3.fill(65);
	A1_2_.elem(idx) = idxrepl3;

	idx = find(A1_2_ == 71);
	uvec idxrepl4(idx.n_elem);
	idxrepl4.fill(67);
	A1_2_.elem(idx) = idxrepl4;

	//A2_1_: replace T with A,replace G with C
	idx = find(A2_1_ == 84);
	uvec idxrepl5(idx.n_elem);
	idxrepl5.fill(65);
	A2_1_.elem(idx) = idxrepl5;

	idx = find(A2_1_ == 71);
	uvec idxrepl6(idx.n_elem);
	idxrepl6.fill(67);
	A2_1_.elem(idx) = idxrepl6;

	//A2_2_: replace T with A,replace G with C
	idx = find(A2_2_ == 84);
	uvec idxrepl7(idx.n_elem);
	idxrepl7.fill(65);
	A2_2_.elem(idx) = idxrepl7;

	idx = find(A2_2_ == 71);
	uvec idxrepl8(idx.n_elem);
	idxrepl8.fill(67);
	A2_2_.elem(idx) = idxrepl8;

	// remove index;
	idx = find((A1_1_ + A2_1_) == (A1_2_ + A2_2_));
	uvec tmp3 = as<uvec>(idxin2);
	uvec idxin22 = tmp3.elem(idx);
	tmp3 = as<uvec>(idxin1);
	uvec idxin11 = tmp3.elem(idx);

	//cout << "break 1: .... " << endl;
	uvec A1_1_r = A1_1_.elem(idx), A2_1_r = A2_1_.elem(idx);
	uvec A1_2_r = A1_2_.elem(idx), A2_2_r = A2_2_.elem(idx);
	//cout << "break 2: .... " << endl;
	uvec chr_tmp = as<uvec>(chr_4use);
	uvec chr_4use_r = chr_tmp.elem(idx);
	uvec bp_tmp = as<uvec>(bp_4use);
	uvec bp_4use_r = bp_tmp.elem(idx);
	
	CharacterVector rsname_4use_r = charv_subset(rsname_4use, idx);


	vec ind(A1_1_r.n_elem);
	idx = find(A1_1_r == A1_2_r);
	ind.ones();
	ind = -ind;
	ind.elem(idx).ones();	

	cout << "Number of matched SNPs (plink file 1 and 2):" << ind.size() << endl;
	//cout << "direction (compare with trait_1000Q_match.R): " << sum(ind) << endl; //compare sum(ind) in R: Height_1000Q_match.R
	//cout << "Size of matched SNPs (remove ambiguous SNPs): " << A1_1_r.n_elem << endl;

	List output = List::create(Rcpp::Named("chr_4use_r") = chr_4use_r,
		Rcpp::Named("A1_1_r") = A1_1_r,
		Rcpp::Named("A1_2_r") = A1_2_r,
		Rcpp::Named("A2_1_r") = A2_1_r,
		Rcpp::Named("A2_2_r") = A2_2_r,
		Rcpp::Named("bp_4use_r") = bp_4use_r,
		Rcpp::Named("rsname_4use_r") = rsname_4use_r,
		Rcpp::Named("indicator") = ind,
		Rcpp::Named("idxinFile1") = idxin11,
		Rcpp::Named("idxinFile2") = idxin22,
		Rcpp::Named("ind") = ind);

	return output;
}

void ReadPlinkFamFile2(std::string stringname, CharacterVector FID, CharacterVector IID, 
	NumericVector pheno, int nrows, int whCol){

	std::ifstream myfile(stringname.c_str());
	std::string line;

	clock_t t1 = clock();

	int nrow_ind = 0;
	vector <string> tmp;

	if (myfile.is_open()){
		while (nrow_ind < nrows){
			if (nrow_ind % 1000 == 0 && nrow_ind != 0){
				cout << nrow_ind << "-th individual " << ",";
				cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
			}

			getline(myfile, line);
			boost::split(tmp, line, boost::is_any_of(" \t *"));

			FID(nrow_ind) = tmp[0];
			IID(nrow_ind) = tmp[1];
			// sex(nrow_ind) = atoi(tmp[4].c_str());
			pheno(nrow_ind) = atof(tmp[4 + whCol].c_str());

			//cout << "value: " << tmp[0] << ";" << tmp[1] << ";" << tmp[2] << ";" << tmp[3] << ";" << tmp[4] 
			//	<< ";" << tmp[5] << ";" << tmp[6] << endl;
			nrow_ind++;
		}
	}
}

Rcpp::List dataLoader(std::string stringname, Rcpp::CharacterVector rsname, arma::uvec A1_r, arma::uvec A2_r,int whCol){
	// plink file 1: stringname (testing set); from testing data
	// A1_r and A2_r: A1 and A2 in file 2; from snpinfo in ssLasso

	List tmp =  match_SNPs2(stringname, rsname, A1_r,A2_r);
	uvec idxinFile1 = tmp["idxinFile1"];
	uvec idxinFile2 = tmp["idxinFile2"];
	vec ind = tmp["indicator"];
	CharacterVector rsname_4use_r = tmp["rsname_4use_r"];
	uvec chr_4use_r = tmp["chr_4use_r"];
	uvec bp_4use_r = tmp["bp_4use_r"];
	//cout << "break 3: .... " << endl;
	uvec A1_1_r = tmp["A1_1_r"], A1_2_r = tmp["A1_2_r"], A2_1_r = tmp["A2_1_r"], A2_2_r = tmp["A2_2_r"];
	
	cout << "## Start matching SNPs in plink file. " << endl;

	string famfile1 = stringname;
	famfile1 += ".fam";
	int N1 = getLineNum(famfile1);

	IntegerVector sex_1(N1);
	NumericVector pheno_1(N1);
	CharacterVector FID_1(N1), IID_1(N1);
	// ReadPlinkFamFile(famfile2, FID_2, IID_2, sex_2, pheno_2, N2);
	ReadPlinkFamFile2(famfile1, FID_1, IID_1, pheno_1, N1, whCol);

	vec y = as<vec>(pheno_1);

	cout << "## Start loading genotype file, ";
	//read file 2 (plink bim bed fam files)
	string bimfile1 = stringname;
	bimfile1 += ".bim";
	int P1 =  getLineNum(bimfile1);
	unsigned* X1tmp = new unsigned[ N1 * P1]; 
  
	cout << "dimension: " << N1 << "X"<<P1 << endl;
	readPlink(stringname,N1, P1, X1tmp);
	//cout << "break 1 ..."<< max(idxinFile1) << endl;

	clock_t t1 = clock();

	arma::Mat<unsigned> X1(X1tmp, N1, P1, false, false);
	X1 = X1.cols(idxinFile1 -1);
	// replace NA with 0 dosage
	X1.replace(3, 0);

	cout << "Genotype loading, Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;

	List output = List::create(
		Rcpp::Named("X1") = X1,
		Rcpp::Named("y") = y,
		Rcpp::Named("rsname_4use_r") = rsname_4use_r,
		Rcpp::Named("chr_4use_r") = chr_4use_r,
		Rcpp::Named("bp_4use_r") = bp_4use_r,
		Rcpp::Named("A1_1_r") = A1_1_r,
		Rcpp::Named("A1_2_r") = A1_2_r,
		Rcpp::Named("A2_1_r") = A2_1_r,
		Rcpp::Named("A2_2_r") = A2_2_r,
		Rcpp::Named("idxinFile2") = idxinFile2,
		Rcpp::Named("ind") = ind);

	return output;
}

// [[Rcpp::export]]
arma::mat pred_REMI(std::string stringname, Rcpp::CharacterVector rsname, arma::uvec A1_r, arma::uvec A2_r,int whCol, arma::vec betah){
	// plink file 1: stringname (testing set); from testing data
	// A1_r and A2_r: A1 and A2 in file 2; from snpinfo in ssLasso (betah)
	// cout << "break 1: .... " << endl;
	cout << "## " << whCol << "-th phenotype in testing set. " << endl;
	List out = dataLoader(stringname, rsname, A1_r, A2_r, whCol);

	cout << "## Start Prediction calc, ";
	clock_t t1 = clock();
	mat X1 = out["X1"];
	vec pheno = out["y"], ind = out["ind"];
	uvec idxinFile2 = out["idxinFile2"];
	vec betah_use = betah.elem(idxinFile2-1);
	betah_use = betah_use % ind;

	vec pred = X1 * betah_use;
	mat corr = cor(pred, pheno);
	cout << "Prediction calc, Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;

	return corr;
}
