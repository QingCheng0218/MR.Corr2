#ifndef GibbsAlpGamEta_ptr_hpp
#define GibbsAlpGamEta_ptr_hpp

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>
#include <mutex>


using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]


class Options_Gibbs3{
public:
  // Constructor definition
  // The complier deciedes which constructor to be called depending on 
  // the number of argument present with the object
  Options_Gibbs3(){
    this -> agm = 0;
    this -> bgm = 0;
    this -> aal = 0;
    this -> bal = 0;
    this -> a = 100;
    this -> b = 1;
    this -> maxIter = 4000;
    this -> thin = 10;
    this -> burnin = 1000;
  }
  
  Options_Gibbs3(double agm, double bgm, double aal, double bal, 
                 double a, double b, uword maxIter, uword thin, uword burnin){
    
    this -> agm = agm;
    this -> bgm = bgm;
    this -> aal = aal;
    this -> bal = bal;
    this -> a = a;
    this -> b = b;
    this -> maxIter = maxIter;
    this -> thin = thin;
    this -> burnin = burnin;
    
  }
  
  double agm;
  double bgm;
  double aal;
  double bal;
  double a;
  double b;
  uword maxIter;
  uword thin;
  uword burnin;
  
};
struct ObjGibbs3{
  ivec Eta;
  vec Etares;
  vec Beta0res;
  vec Beta1res;
  vec Sgga2Res;
  vec Sgal2Res;
  imat EtaAll;
};


class paraBlock_GamAlpEta{
public:
  int current_idx=0;
  // uword Ngene_active;
  int n_thread = 1;
  umat block_inf;
  
  uword nblocks;
  
  int conspar;
  ivec Eta;
  vec beta0sig, beta1sig, beta0mean, beta1mean, Mu2, MuA2;
  double beta0, beta1, logW, sgga2, sgal2;
  
  field<vec> F4se1, F4se2, F4Gammah, F4gammah, F4ginvsg2, F4GinvsG2;
  field<vec> F4mu, F4muA, F4Rinsgmu,  F4RinsGmu, F4RinsGmuA;
  field<mat> F4Rins, F4Rins2, F4Rblock, F4insGRinsG, F4insgRinsg;
  // field<mat> F4sgRinvsg, F4sGRinvsG, F4sgRsg, F4sGRsG;
  
  
  paraBlock_GamAlpEta(uword &nblocks, field<vec> &F4se1, field<vec> &F4se2, 
                      field<vec> &F4gammah,
                      field<vec> &F4Gammah, field<vec> &F4ginvsg2,
                      field<vec> &F4GinvsG2,  field<vec> &F4mu, field<vec> &F4muA, field<vec> &F4Rinsgmu, field<vec> &F4RinsGmu,
                      field<vec> &F4RinsGmuA, field<mat> &F4insgRinsg, field<mat> &F4insGRinsG,
                      field<mat> &F4Rins, field<mat> &F4Rins2, field<mat> &F4Rblock,
                      // field<mat> &F4sgRinvsg,field<mat> &F4sGRinvsG, field<mat> &F4sgRsg,field<mat> &F4sGRsG,
                      arma::vec &beta0sig, arma::vec &beta1sig, arma::vec &beta0mean, arma::vec &beta1mean, arma::ivec Eta,
                      arma::vec &Mu2, arma::vec &MuA2,
                      double &beta0, double &beta1, double &logW, double &sgga2, double &sgal2, const int &conspar){
    this -> nblocks = nblocks;
    this -> F4se1 = F4se1;
    this -> F4se2 = F4se2;
    this -> F4gammah = F4gammah;
    this -> F4Gammah = F4Gammah;
    this -> F4Rins = F4Rins;
    this -> F4Rins2 = F4Rins2;
    this -> F4Rblock = F4Rblock;
    // this -> F4sgRinvsg = F4sgRinvsg;
    // this -> F4sGRinvsG = F4sGRinvsG;
    // this -> F4sgRsg = F4sgRsg;
    // this -> F4sGRsG = F4sGRsG;
    this -> F4ginvsg2 = F4ginvsg2;
    this -> F4GinvsG2 = F4GinvsG2;
    this -> F4insgRinsg = F4insgRinsg;
    this -> F4insGRinsG = F4insGRinsG;
    this -> F4mu = F4mu;
    this -> F4muA = F4muA;
    this -> F4Rinsgmu = F4Rinsgmu;
    this -> F4RinsGmu = F4RinsGmu;
    this -> F4RinsGmuA = F4RinsGmuA;
    this -> beta0sig = beta0sig;
    this -> beta0mean = beta0mean;
    this -> beta1sig = beta1sig;
    this -> beta1mean = beta1mean;
    this -> Mu2 = Mu2;
    this -> MuA2 = MuA2;
    this -> Eta = Eta;
    this -> beta0 = beta0;
    this -> beta1 = beta1;
    this -> logW = logW;
    this -> sgga2 = sgga2;
    this -> sgal2 = sgal2;
    this -> conspar = conspar;
  }
  
  
  int  next_GamAlpEta();
  void loop_by_block_gibbs_GamAlpEta(int i);
  void update_by_thread_GamAlpEta(int thread_id);
  
};

#endif /* GibbsAlpGamEta_ptr_hpp */
