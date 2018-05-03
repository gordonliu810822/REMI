//  Created by Jin Liu on 23/04/2018.
//  Copyright © 2018年 Jin Liu. All rights reserved.

#include "plinkfun.hpp"
#include "correlationcal.hpp"
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <stdio.h>
#include "genoLoader_pred.hpp"
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

using namespace Rcpp;
using namespace std;
using namespace arma;

/*
// [[Rcpp::export]]
void ReadPlinkFamFile2(std::string stringname, Rcpp::CharacterVector FID, Rcpp::CharacterVector IID,
                       Rcpp::NumericVector pheno, int nrows, int whCol){

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
            //    << ";" << tmp[5] << ";" << tmp[6] << endl;
            nrow_ind++;
        }
    }
}*/

Rcpp::List ReadDataFromFile3(std::string stringname1, std::string stringname2, int whCol) {
    string famfile = stringname1;
    famfile += ".fam";
    string bimfile = stringname1;
    bimfile += ".bim";

    cout << endl;
    cout << "Start loading genotype: " << endl;

    clock_t t1 = clock();
    int N =  getLineNum(famfile);
    int P =  getLineNum(bimfile);
    //cout << bimfile << endl;
    //cout << "N: " << N << endl;
    //cout << "P: " << P << endl;
    arma::mat bimmat;
    bimmat.load(bimfile,arma::csv_ascii);
    arma::vec c1 = bimmat.col(0);
    cout << bimmat.size() << endl;
    arma::vec index(22);
    double sum2 = 0;

    for(int i=1; i <= 22; i++ ){
        index[i-1] = sum(c1 == i);
        sum2 += index[i-1];
        cout <<"number of snps in chromsome "<<i <<" =" << index[i-1] << endl;
    }
    unsigned* X = new unsigned[ N * P];

    readPlink(stringname1,N, P, X);
    cout << "Finished loading genotype " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
    arma::Mat<unsigned>* Xdata = new arma::Mat<unsigned>(X, N, P, false,false);
    cout << "Dimension of genotype: " << size(*Xdata) << endl;
    //cout << "Size of avbIndex: " << avbIndex.size() << endl;
    //cout << "min of avbIndex: " << min(avbIndex) << endl;
    //cout << "max of avbIndex: " << max(avbIndex) << endl;
    // load pheno in file 2 (fam file 2)
    string famfile2 = stringname2;
    //famfile2 += ".fam";
    int N2 = getLineNum(famfile2);

    IntegerVector sex_2(N2);
    NumericVector pheno_2(N2);
    CharacterVector FID_2(N2), IID_2(N2);
    // ReadPlinkFamFile(famfile2, FID_2, IID_2, sex_2, pheno_2, N2);
    ReadPlinkFamFile2(famfile2, FID_2, IID_2, pheno_2, N2, whCol);
    
    vec y = as<vec>(pheno_2);
    vec yc = y - mean(y);
    double sigma2 = var(yc);
    
    cout << endl;
    cout <<  "Start calculating inner product: " ;
    t1 = clock();
    arma::vec ytilde = calInnerProd(Xdata, index, yc);
    cout << "Finished calculation in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
    cout << "Length of ytilde: " << size(ytilde) << endl;
    //cout << "calCorr done." << endl;
    arma::Mat<unsigned> Xraw(X, N, P, false, false);
    
    List obj= List::create(Rcpp::Named("X") = Xraw,
                           Rcpp::Named("ytilde") = ytilde,
                           Rcpp::Named("sigma2")=sigma2);
    //obj.X = Xraw;
    //obj.Sigma = SigmaMat;


    return obj;
}


