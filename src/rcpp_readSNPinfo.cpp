//  Created by Jin Liu on 23/04/2018.
//  Copyright © 2018年 Jin Liu. All rights reserved.

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdio.h>
#include <bitset>
#include <math.h>
#include "plinkfun.hpp"
#include "correlationcal.hpp"
#include "readPlink.hpp"

#include<iostream>

#define MAX_LEN 20

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


void ReadSNPinfo(std::string stringname, IntegerVector A1, IntegerVector A2, CharacterVector rsname,
	IntegerVector chr, IntegerVector bp, NumericVector morgan, int N)
{
	FILE *stream;
	int ch, b;
	double mor;
	char s[MAX_LEN + 1], efa, nefa;

	stream = fopen(stringname.c_str(), "r");
	clock_t t1 = clock();
	//int i = 0;
	/* Put in various data. */
	for (int i = 0; i < N; i++){
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

void Read_summarystat(std::string stringname, IntegerVector SA1, IntegerVector SA2, CharacterVector rsname,
	NumericVector betah, NumericVector s2, NumericVector sample_size, int N){
	
	FILE *stream;
	stream = fopen(stringname.c_str(), "r");

	char s[MAX_LEN + 1], efa, nefa;
	double bhat, se2, samsize;

	clock_t t1 = clock();
	for (int i = 0; i < N; i++){
		if (i % 100000 == 0 && i != 0){
			cout << i << "-th SNP" << ",";
			cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
		}

		fscanf(stream, "%s %c %c %lf %lf %lf", &s[0], &efa, &nefa, &bhat, &se2, &samsize);
		
		rsname(i) = s;
		//cout << s << efa << nefa << bhat << se2 << samsize << ";" << endl;
		SA1(i) = (int)efa;
		SA2(i) = (int)nefa;
		betah(i) = bhat;
		s2(i) = se2;
		sample_size(i) = samsize;
		
	}

}

Rcpp::List preprocessn(std::string stringname1, std::string stringname2){
	//Panel SNP info
	//rsname2, rsname in summary stat file

	string bimfile = stringname2;
	bimfile += ".bim";

	int N = getLineNum(bimfile);
	cout << "Number of SNPs (panel):" << N << endl;

	
	IntegerVector A1(N), A2(N);
	CharacterVector rsname1(N);
	IntegerVector chr(N), bp(N);
	NumericVector morgan(N);

	//read from SNPinfo file (pass as pointer)
	cout << endl;
	cout << "Start loading SNP info:" << endl;
	clock_t t1 = clock();
	ReadSNPinfo(bimfile, A1, A2, rsname1, chr, bp, morgan, N);
	cout << "Finish loading SNP info in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
	// summary stat SNP info, betah, s2, sample_size
	int Ns = getLineNum(stringname1);

	cout << "Number of SNPs (summarystat):" << Ns << endl;

	IntegerVector SA1(Ns), SA2(Ns);
	CharacterVector rsname2(Ns);
	NumericVector betah(Ns), s2(Ns), sample_size(Ns);

	cout << endl;
	printf("Start loading summary stat: \n");
    t1 = clock();
	Read_summarystat(stringname1, SA1, SA2, rsname2, betah, s2, sample_size, Ns);
	cout << "Finish loading summary stat in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;

	//mathcing panel SNP with summary stat and correct direction of minor allele
	CharacterVector rs_inter = intersect(rsname1, rsname2);
	IntegerVector idxin1 = match(rsname1, rs_inter);
	CharacterVector rsname_4use = rsname1[Rcpp::is_na(idxin1) == false]; // rsname in both panel and ss in the order of panel.

	//cout << "size of intersect: " << rsname_4use.size() << endl;

	/*for (int i = 0; i < 10; i++){
	cout << "rsnames_4use: " << rsname_4use(i) << idxin1(i) <<endl;
	}*/
	cout << endl;
	printf("Start matching SNPs with summary stat. \n");
	//match snps (rsname_4use; rsname2: trait)
	IntegerVector idxin2 = match(rsname_4use, rsname2);  //index for SNPs in ss

	//cout << "test" << idxin2(0) << endl;
	//cout << "test" << idxin2(2) << endl;
	//cout << "size of idxin2: " << idxin2.size() << endl;

	
	uvec idx = as<uvec>(idxin2) -1;
	fvec tmp = as<fvec>(betah);
	fvec bhat = tmp.elem(idx);
	tmp = as<fvec>(s2);
	fvec se2 = tmp.elem(idx);

	tmp = as<fvec>(sample_size);
	fvec samsize = tmp.elem(idx);
	uvec tmp2 = as<uvec>(SA1);
	uvec SA1_ = tmp2.elem(idx);
	tmp2 = as<uvec>(SA2);
	uvec SA2_ = tmp2.elem(idx);

	//match snps (rsname_4use; rsname1: panel)
	printf("Start matching SNPs with panel. \n");
	IntegerVector idxin3 = match(rsname_4use, rsname1); //index for SNPs in panel SNPs
	
	//compare direction
	cout << "Size of matched SNPs: " << rsname_4use.size() << endl; //compare keepIndx in R: Height_1000Q_match.R
	//cout << "min of idxin3: " << min(idxin3) << endl;
	//cout << "max of idxin3: " << max(idxin3) << endl;


	idx = as<uvec>(idxin3) -1;
	//mat R = tmp1.rows(idx);
	tmp2 = as<uvec>(A1);
	uvec A1_ = tmp2.elem(idx);
	tmp2 = as<uvec>(A2);
	uvec A2_ = tmp2.elem(idx);
	//ascii: (A:65; C:67; G:71; T:84) (a:97;c:99,g:103;t:116)
	//replace lower letter to upper letter in SA1_ and SA2_;
	idx = find(SA1_ > 85);
	SA1_.elem(idx) = SA1_.elem(idx) - 32;
	idx = find(SA2_ > 85);
	SA2_.elem(idx) = SA2_.elem(idx) - 32;

	//compare A1_ SA1_, A2_ SA2_
	//A1_: replace T with A,replace G with C
	idx = find(A1_ == 84);
	uvec idxrepl1(idx.n_elem);
	idxrepl1.fill(65);
	A1_.elem(idx) = idxrepl1;

	idx = find(A1_ == 71);
	uvec idxrepl2(idx.n_elem);
	idxrepl2.fill(67);
	A1_.elem(idx) = idxrepl2;

	//SA1_: replace T with A,replace G with C
	idx = find(SA1_ == 84);
	uvec idxrepl3(idx.n_elem);
	idxrepl3.fill(65);
	SA1_.elem(idx) = idxrepl3;

	idx = find(SA1_ == 71);
	uvec idxrepl4(idx.n_elem);
	idxrepl4.fill(67);
	SA1_.elem(idx) = idxrepl4;

	//A2_: replace T with A,replace G with C
	idx = find(A2_ == 84);
	uvec idxrepl5(idx.n_elem);
	idxrepl5.fill(65);
	A2_.elem(idx) = idxrepl5;

	idx = find(A2_ == 71);
	uvec idxrepl6(idx.n_elem);
	idxrepl6.fill(67);
	A2_.elem(idx) = idxrepl6;

	//SA1_: replace T with A,replace G with C
	idx = find(SA2_ == 84);
	uvec idxrepl7(idx.n_elem);
	idxrepl7.fill(65);
	SA2_.elem(idx) = idxrepl7;

	idx = find(SA2_ == 71);
	uvec idxrepl8(idx.n_elem);
	idxrepl8.fill(67);
	SA2_.elem(idx) = idxrepl8;

	//remove index;
	idx = find((A1_ + A2_) == (SA1_ + SA2_));
	uvec tmp1 = as<uvec>(idxin3);
	uvec idxin4 = tmp1.elem(idx);
	uvec A1_r = A1_.elem(idx), A2_r = A2_.elem(idx);
	uvec SA1_r = SA1_.elem(idx), SA2_r = SA2_.elem(idx);

	IntegerVector idxtmp = wrap(idx);
	CharacterVector rst = rsname_4use[idxtmp];

	//cout << "Size of idxin4: " << idxin4.size() << endl;
	//cout << "min of idxin4: " << min(idxin4) << endl;
	//cout << "max of idxin4: " << max(idxin4) << endl;

	fvec bh = bhat.elem(idx);
	fvec s2_ = se2.elem(idx);
	fvec sample_size_ = samsize.elem(idx);

	fvec ind(A1_r.n_elem);
	idx = find(A1_r == SA1_r);
	ind.ones();
	ind = -ind;
	ind.elem(idx).ones();
	cout << "direction (compare with trait_1000Q_match.R): " << sum(ind) << endl; //compare sum(ind) in R: Height_1000Q_match.R
	cout << "Size of matched SNPs (remove ambiguous SNPs): " << A1_r.n_elem << endl;

	bh = bh % ind;

	List output = List::create(Rcpp::Named("betah") = bh,
		Rcpp::Named("s2") = s2_,
		Rcpp::Named("idxin") = idxin4,
		Rcpp::Named("sample_size") = sample_size_,
		Rcpp::Named("rsname") = rst,
		Rcpp::Named("A1_r") = A1_r,
		Rcpp::Named("A2_r") = A2_r);
	return output;
}


