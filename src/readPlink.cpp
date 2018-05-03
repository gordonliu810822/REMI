//  Created by Jin Liu on 23/04/2018.
//  Copyright © 2018年 Jin Liu. All rights reserved.

#include "plinkfun.hpp"
#include "readPlink.hpp"
#include "correlationcal.hpp"

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


ObjXCorr ReadDataFromFile(string stringname, arma::uvec avbIndex){
	int bandwidth = 200;
	float shrinkagefactor = 0.9;
	return ReadDataFromFile(stringname, avbIndex, bandwidth, shrinkagefactor);
}



ObjXCorr ReadDataFromFile(string stringname){
	int bandwidth = 200;
	float shrinkagefactor = 0.9;
	arma::uvec avbIndex;
	avbIndex.reset();
	return ReadDataFromFile(stringname, avbIndex, bandwidth, shrinkagefactor);

}

ObjXCorr ReadDataFromFile(string stringname, arma::uvec avbIndex, int bandwidth, float shrinkagefactor) {
  string famfile = stringname;
  famfile += ".fam";
  string bimfile = stringname;
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
  
  readPlink(stringname,N, P, X);
  cout << "Finished loading genotype " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
  arma::Mat<unsigned>* Xdata = new arma::Mat<unsigned>(X, N, P, false,false);
  cout << "Dimension of genotype: " << size(*Xdata) << endl;
  //cout << "Size of avbIndex: " << avbIndex.size() << endl;
  //cout << "min of avbIndex: " << min(avbIndex) << endl;
  //cout << "max of avbIndex: " << max(avbIndex) << endl;
  cout << endl;
  cout <<  "Start calculating correlation: " ;
  t1 = clock();
  arma::fmat corMat = calCorr(Xdata, index, avbIndex, bandwidth, shrinkagefactor);
  cout << "Finished calculation in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
  cout << "Dimension of banded correlation: " << size(corMat) << endl;
  //cout << "calCorr done." << endl;
  ObjXCorr obj;
  obj.X = *Xdata;
  obj.corr = corMat;

  
  return obj;
}


List ReadDataFromFile2(string stringname, arma::uvec avbIndex, int bandwidth, float shrinkagefactor) {
    string famfile = stringname;
    famfile += ".fam";
    string bimfile = stringname;
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
    
    readPlink(stringname,N, P, X);
    cout << "Finished loading genotype " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
    arma::Mat<unsigned>* Xdata = new arma::Mat<unsigned>(X, N, P, false,false);
    cout << "Dimension of genotype: " << size(*Xdata) << endl;
    //cout << "Size of avbIndex: " << avbIndex.size() << endl;
    //cout << "min of avbIndex: " << min(avbIndex) << endl;
    //cout << "max of avbIndex: " << max(avbIndex) << endl;
    cout << endl;
    cout <<  "Start calculating covariance: " ;
    t1 = clock();
    arma::fmat SigmaMat = calSigma(Xdata, index, avbIndex, bandwidth, shrinkagefactor);
    cout << "Finished calculation in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
    cout << "Dimension of banded covariance: " << size(SigmaMat) << endl;
    //cout << "calCorr done." << endl;
    arma::Mat<unsigned> Xraw(X, N, P, false, false);
    List obj= List::create(Rcpp::Named("X") = Xraw,
                           Rcpp::Named("Sigma") = SigmaMat);
    //obj.X = Xraw;
    //obj.Sigma = SigmaMat;
    
    
    return obj;
}

// [[Rcpp::export]]
Rcpp::List ReadPlinkFile(std::string stringname){
    string famfile = stringname;
    famfile += ".fam";
    string bimfile = stringname;
    bimfile += ".bim";
    long long N = getLineNum(famfile);
    long long P = getLineNum(bimfile);
    arma::mat bimmat;
    bimmat.load(bimfile, arma::csv_ascii);
    arma::vec c1 = bimmat.col(0);
    cout << bimmat.size() << endl;
    arma::vec index(22);
    double sum2 = 0;
    for (int i = 1; i <= 22; i++){
        index[i - 1] = sum(c1 == i);
        sum2 += index[i - 1];
        cout << "number of snps in chromsome " << i << " =" << index[i - 1] << endl;
    }
    unsigned* X = new unsigned[N * P];
    readPlink(stringname, N, P, X);
    arma::Mat<unsigned>* Xdata = new arma::Mat<unsigned>(X, N, P, false, false);
    
    Xdata->replace(3, 0);
    Mat<unsigned> Xout = *Xdata;
    
    List out2 = List::create(Rcpp::Named("X") = Xout);
    return out2;
}
