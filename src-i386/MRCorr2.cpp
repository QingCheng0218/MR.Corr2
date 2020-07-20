#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>
#include <mutex>
#include "pdsoft.hpp"
#include "GibbsAlpGamEta_ptr.hpp"
#include "ReadGeneFile.hpp"
#include "CalCorr.hpp"
#include "data_loader.hpp"
#include <ctime>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]


ObjGibbs3 gibbs3group_blockpar(arma::vec &gammah, arma::vec &Gammah, arma::vec &se1, arma::vec &se2, 
                               arma::mat R, arma::umat block_inf,  int coreNum, const int&conspar, Options_Gibbs3* opts)
{
  
  // ----------------------------------------------------------------------
  // check number of input arguments
  double agm = opts -> agm;
  double bgm = opts -> bgm;
  double aal = opts -> aal;
  double bal = opts -> bal;
  double a = opts -> a;
  double b = opts -> b;
  uword maxIter = opts -> maxIter;
  uword thin = opts -> thin;
  uword burnin = opts -> burnin;
  
  // ----------------------------------------------------------------------
  // initial values
  int p = Gammah.n_elem;
  double sgga2 = 1; double sgal2 = 1; double beta0 = 0.01; double w = 0.1;
  double beta1 = 0.01; 
  int numsave = maxIter / thin;
  vec Beta0res = ones(numsave, 1);
  vec Beta1res = ones(numsave, 1);
  vec Etares = ones(numsave, 1);
  
  
  vec mu = 0.01*ones(p, 1);
  vec muA = 0.01*ones(p, 1);
  
  uword nblocks = block_inf.n_rows;
  ivec Eta;
  if(conspar==0||conspar==2){
    Eta = zeros<ivec>(nblocks, 1);
  }else if(conspar==1){
    Eta = ones<ivec>(nblocks, 1);
  }
  ivec Eta1 = zeros<ivec>(nblocks, 1);
  // ----------------------------------------------------------------------
  // cout << endl;
  // cout << "Preprocessing Process:" << endl;
  // clock_t t1 = clock();
  
  ivec NB = zeros<ivec>(nblocks, 1);
  
  field<vec> F4se1(nblocks, 1), F4se2(nblocks, 1), F4gammah(nblocks, 1), F4Gammah(nblocks, 1), F4mu(nblocks, 1), F4muA(nblocks, 1);
  field<vec> F4sg2(nblocks, 1), F4sG2(nblocks, 1), F4GinvsG2(nblocks, 1), F4ginvsg2(nblocks, 1);
  field<vec> F4diaginsGRinsG(nblocks, 1), F4diaginsgRinsg(nblocks, 1), F4RinsGmu(nblocks, 1), F4Rinsgmu(nblocks, 1), F4RinsGmuA(nblocks, 1);
  field<mat> F4insGRinsG(nblocks, 1), F4insgRinsg(nblocks, 1), F4insgRinsG(nblocks, 1), F4Rins(nblocks, 1), F4Rins2(nblocks, 1);
  field<mat> F4Rblock(nblocks, 1);
  // field<mat> F4sgRinvsg(nblocks, 1), F4sGRinvsG(nblocks, 1), F4sgRsg(nblocks, 1), F4sGRsG(nblocks, 1);
  
  // cout<<"nblocks:"<< nblocks<< endl;
  for (int nn = 0; nn < (int)(nblocks); nn = nn+1){
    
    NB[nn] = block_inf(nn, 1) - block_inf(nn, 0) + 1;
    // cout << block_inf(nn, 0) << block_inf(nn, 1)<< endl;
    vec se1_block = se1.subvec(block_inf(nn, 0), block_inf(nn, 1));
    vec se2_block = se2.subvec(block_inf(nn, 0), block_inf(nn, 1));
    vec sg2_block = pow(se1_block, 2);
    vec sG2_block = pow(se2_block, 2);
    // cout<< "nn:" << nn << "check error 1" << endl;
    
    vec mu_block = mu.subvec(block_inf(nn, 0), block_inf(nn, 1));
    vec muA_block = muA.subvec(block_inf(nn, 0), block_inf(nn, 1));
    
    vec bh1_block = gammah.subvec(block_inf(nn, 0), block_inf(nn, 1));
    vec bh2_block = Gammah.subvec(block_inf(nn, 0), block_inf(nn, 1));
    
    mat R_block = R.submat(block_inf(nn, 0), block_inf(nn, 0), block_inf(nn, 1), block_inf(nn, 1));
    F4Rblock(nn, 0) = R_block;
    
    F4mu(nn, 0) = mu_block;
    F4muA(nn, 0) = muA_block;
    F4se1(nn, 0) = se1_block;
    F4se2(nn, 0) = se2_block;
    
    F4sg2(nn, 0) = sg2_block;
    F4sG2(nn, 0) = sG2_block;
    
    F4gammah(nn, 0) = bh1_block;
    F4Gammah(nn, 0) = bh2_block;
    F4GinvsG2(nn, 0) = bh2_block / sG2_block;
    F4ginvsg2(nn, 0) = bh1_block / sg2_block;
    
    // cout<< "nn:" << nn << "check error 2" << endl;
    
    F4insGRinsG(nn, 0) = diagmat(1. / se2_block)*R_block*diagmat(1. / se2_block);
    F4insgRinsg(nn, 0) = diagmat(1. / se1_block)*R_block*diagmat(1. / se1_block);
    
    F4Rins(nn, 0) = R_block*diagmat(1 / se1_block);
    F4Rins2(nn, 0) = R_block*diagmat(1 / se2_block);
    
    F4RinsGmu(nn, 0) = R_block*diagmat(1 / se2_block)*mu_block;
    F4Rinsgmu(nn, 0) = R_block*diagmat(1 / se1_block)*mu_block;
    F4RinsGmuA(nn, 0) = R_block*diagmat(1. / se2_block)*muA_block;
    
    
    // F4sgRinvsg(nn, 0) = diagmat(se1_block)*R_block*diagmat(1. / se1_block);
    // F4sGRinvsG(nn, 0) = diagmat(se2_block)*R_block*diagmat(1. / se2_block);
    // F4sgRsg(nn, 0) = diagmat(se1_block)*R_block*diagmat(se1_block);
    // F4sGRsG(nn, 0) = diagmat(se2_block)*R_block*diagmat(se2_block);
    // cout<< "nn:" << nn << "check error 3" << endl;
  }
  
  // cout << "Finish the Preprocessing Process in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
  
  int l = 0;
  
  // cout << endl;
  // cout << "Start gibbs sampling Process:" << endl;
  // t1 = clock();
  // 
  
  vec beta0mean = zeros(nblocks, 1); // mean for each block.
  vec beta0sig = zeros(nblocks, 1); // sigma for each block.
  vec beta1mean = zeros(nblocks, 1); // mean for each block.
  vec beta1sig = zeros(nblocks, 1); // sigma for each block.
  vec Mu2 = zeros(nblocks, 1); // sgga2 for each block.
  vec MuA2 = zeros(nblocks, 1); // sgal2 for each block.
  
  for(int iter = 0; iter < (int_fast32_t)(maxIter + burnin); iter ++){
    double logW = log(w /( 1 - w));
    
    // ------------------------------------------------------------------
    // set parallel computation for gamma, alpha, eta
    paraBlock_GamAlpEta parobj_GamAlpEta(nblocks, F4se1, F4se2, F4gammah, F4Gammah, F4ginvsg2,
                                         F4GinvsG2,  F4mu, F4muA, F4Rinsgmu, F4RinsGmu, 
                                         F4RinsGmuA, F4insgRinsg, F4insGRinsG,
                                         F4Rins, F4Rins2, F4Rblock, 
                                         // F4sgRinvsg, F4sGRinvsG, F4sgRsg, F4sGRsG,
                                         beta0sig, beta1sig, beta0mean, beta1mean, Eta,
                                         Mu2, MuA2,
                                         beta0, beta1, logW, sgga2, sgal2, conspar);
    
    
    const int n_thread = coreNum;
    std::vector<std::thread> threads(n_thread);
    for(int i_thread = 0; i_thread < n_thread; i_thread++){
      threads[i_thread] = std::thread(&paraBlock_GamAlpEta::update_by_thread_GamAlpEta, &parobj_GamAlpEta, i_thread);
    }
    
    
    for(int i = 0; i < n_thread; i++){
      threads[i].join();
    }
    
    // save the parallel result
    beta0sig =  parobj_GamAlpEta.beta0sig;
    beta1sig =  parobj_GamAlpEta.beta1sig;
    
    beta0mean =  parobj_GamAlpEta.beta0mean;
    beta1mean =  parobj_GamAlpEta.beta1mean;
    
    Mu2 = parobj_GamAlpEta.Mu2;
    MuA2 = parobj_GamAlpEta.MuA2;
    
    F4mu = parobj_GamAlpEta.F4mu;
    F4muA = parobj_GamAlpEta.F4muA;
    F4Rinsgmu = parobj_GamAlpEta.F4Rinsgmu;
    F4RinsGmu = parobj_GamAlpEta.F4RinsGmu;
    F4RinsGmuA = parobj_GamAlpEta.F4RinsGmuA;
    
    Eta = parobj_GamAlpEta.Eta;
    // ------------------------------------------------------------------
    
    // ------------------ //
    // update omega
    // ------------------ //
    Eta1 = 1 - Eta;
    double wpa1, wpa2;
    wpa1 = a + sum(Eta1);
    wpa2 = b + Eta1.n_elem - sum(Eta1);
    w = 1 - R::rbeta(wpa1, wpa2);
    
    // cout<<"Eta1:"<< Eta1.subvec(0, 5).t() << "Eta:" << Eta.subvec(0, 5).t()<< endl;
    // cout<<"sum(Eta1):"<< sum(Eta1) << "sum(Eta):" << sum(Eta) << "Eta.n_elem:" << Eta.n_elem << endl;
    // cout<<"wpa1:"<< wpa1 << "wpa2:" << wpa2 << endl;
    // ------------------ //
    // update beta0
    // ------------------ //
    double sig2b0, mub0;
    if(sum(Eta)==0){
      beta0 = 0;
    }else{
      sig2b0 = 1. / sum(beta0sig);
      mub0 = sum(beta0mean) * sig2b0;
      beta0 = mub0 + randn()*sqrt(sig2b0);
      // beta0 = mub0;
    }
    // cout<<"Mb0:"<< beta0 << "--sigb0:" << sqrt(sig2b0) << endl;
    // ------------------ //
    // update beta1
    // ------------------ //
    double sig2b1, mub1;
    if(sum(Eta)==(int)nblocks){
      beta1 = 0;
    }else{
      sig2b1 = 1. / sum(beta1sig);
      mub1 = sum(beta1mean) * sig2b1;
      beta1 = mub1 + randn()*sqrt(sig2b1);
      // beta1 = mub1;
      // cout <<"check error 2" << endl;
    }
    // cout<<"Mb1:"<< beta1 << "--sigb1:" << sqrt(sig2b1) << endl;
    // ------------------ //
    // update sgga2
    // ------------------ //
    double tagm, tbgm, taal, tbal;
    tagm = agm + p / 2;
    tbgm = sum(Mu2) / 2 + bgm;
    sgga2 =  1 / randg<double>(distr_param(tagm, 1/tbgm));
    // sgga2 = tbgm;
    // ------------------ //
    // update sgal2
    // ------------------ //
    taal = aal + p / 2;
    tbal = sum(MuA2) / 2 + bal;
    sgal2 =  1 / randg<double>(distr_param(taal, 1/tbal));
    // sgal2 = tbal;
    
    // cout<<"sgga2:" << sgga2 <<"--tagm--"<< tagm << endl;
    // cout<<"sgal2:" << sgal2 <<"--taal--"<< taal << endl;
    // cout<<"beta0--"<<beta0<<"--beta1--"<<beta1 << endl;
    // cout<<"--------------------------------------------"<<endl;
    // cout<<endl;
    if(iter >= (int)burnin){
      if((iter - burnin) % thin ==0){
        // cout << "Eta:" << sum(Eta)<< endl;
        Beta0res[l] = beta0;
        Beta1res[l] = beta1;
        Etares[l] = 1 - (float)sum(Eta)/(float)Eta.n_elem;
        l += 1;
      }
    }
    
    // if(iter%500==0){
    //   cout << "the iteration number is:" << iter << endl;
    // }
    
    
  }
  

  ObjGibbs3 obj;
  obj.Eta = Eta1;
  obj.Etares = Etares;
  obj.Beta0res = Beta0res;
  obj.Beta1res = Beta1res;
  
  return obj;
}

//[[Rcpp::export]]
List MRCorr2Sim(arma::vec &gammah, arma::vec &Gammah, arma::vec &se1, arma::vec &se2,
                           arma::mat R, arma::umat block_inf, int coreNum,  SEXP opts = R_NilValue)
{
  Options_Gibbs3* lp_opt = NULL;
  if (!Rf_isNull(opts)){
    Rcpp::List opt(opts);
    lp_opt = new Options_Gibbs3(opt["agm"], opt["bgm"], opt["aal"], opt["bal"], opt["a"], opt["b"], opt["maxIter"], opt["thin"], opt["burnin"]);
  }
  if (Rf_isNull(opts)){
    lp_opt = new Options_Gibbs3();
  }
  const int conspar = 0;
  
  ObjGibbs3 obj = gibbs3group_blockpar(gammah, Gammah, se1, se2, R, block_inf, coreNum, conspar, lp_opt);
  List output = List::create(
    Rcpp::Named("Eta") = Rcpp::wrap(obj.Eta),
    Rcpp::Named("Etares") = Rcpp::wrap(obj.Etares),
    Rcpp::Named("Beta0res") = Rcpp::wrap(obj.Beta0res),
    Rcpp::Named("Beta1res") = Rcpp::wrap(obj.Beta1res)
  );
  return output;
}

ObjGibbs3 gibbs3group_blockparReal(arma::vec &gammah, arma::vec &Gammah, arma::vec &se1, arma::vec &se2, 
                                   arma::field<mat> F4Rblock, arma::umat block_inf,  int coreNum, const int&conspar, Options_Gibbs3* opts)
{
  
  // ----------------------------------------------------------------------
  // check number of input arguments
  double agm = opts -> agm;
  double bgm = opts -> bgm;
  double aal = opts -> aal;
  double bal = opts -> bal;
  double a = opts -> a;
  double b = opts -> b;
  uword maxIter = opts -> maxIter;
  uword thin = opts -> thin;
  uword burnin = opts -> burnin;
  
  // ----------------------------------------------------------------------
  // initial values
  int p = Gammah.n_elem;
  double sgga2 = 1; double sgal2 = 1; double beta0 = 0.01; double w = 0.1;
  double beta1 = 0.01; 
  int numsave = maxIter / thin;
  vec Beta0res = ones(numsave, 1);
  vec Beta1res = ones(numsave, 1);
  vec Etares = ones(numsave, 1);
  vec Sgga2Res = ones(numsave, 1);
  vec Sgal2Res = ones(numsave, 1);
  
  vec mu = 0.01*ones(p, 1);
  vec muA = 0.01*ones(p, 1);
  
  uword nblocks = block_inf.n_rows;
  
  imat EtaAll = ones<imat>(nblocks, numsave);
  ivec Eta;
  if(conspar==0||conspar==2){
    Eta = zeros<ivec>(nblocks, 1);
  }else if(conspar==1){
    Eta = ones<ivec>(nblocks, 1);
  }
  ivec Eta1 = zeros<ivec>(nblocks, 1);
  
  
  // ----------------------------------------------------------------------
  cout << endl;
  cout << "Preprocessing Process:" << endl;
  clock_t t1 = clock();
  
  ivec NB = zeros<ivec>(nblocks, 1);
  
  field<vec> F4se1(nblocks, 1), F4se2(nblocks, 1), F4gammah(nblocks, 1), F4Gammah(nblocks, 1), F4mu(nblocks, 1), F4muA(nblocks, 1);
  field<vec> F4sg2(nblocks, 1), F4sG2(nblocks, 1), F4GinvsG2(nblocks, 1), F4ginvsg2(nblocks, 1);
  field<vec> F4diaginsGRinsG(nblocks, 1), F4diaginsgRinsg(nblocks, 1), F4RinsGmu(nblocks, 1), F4Rinsgmu(nblocks, 1), F4RinsGmuA(nblocks, 1);
  field<mat> F4insGRinsG(nblocks, 1), F4insgRinsg(nblocks, 1), F4insgRinsG(nblocks, 1), F4Rins(nblocks, 1), F4Rins2(nblocks, 1);
  field<mat> F4sGRsG(nblocks, 1), F4BB(nblocks, 1);
  
  
  for (int nn = 0; nn < (int)(nblocks); nn = nn+1){
    
    NB[nn] = block_inf(nn, 1) - block_inf(nn, 0) + 1;
    
    vec se1_block = se1.subvec(block_inf(nn, 0), block_inf(nn, 1));
    vec se2_block = se2.subvec(block_inf(nn, 0), block_inf(nn, 1));
    vec sg2_block = pow(se1_block, 2);
    vec sG2_block = pow(se2_block, 2);
    
    vec mu_block = mu.subvec(block_inf(nn, 0), block_inf(nn, 1));
    vec muA_block = muA.subvec(block_inf(nn, 0), block_inf(nn, 1));
    
    vec bh1_block = gammah.subvec(block_inf(nn, 0), block_inf(nn, 1));
    vec bh2_block = Gammah.subvec(block_inf(nn, 0), block_inf(nn, 1));
    
    mat R_block =  symmatu(F4Rblock(nn, 0));
    
    F4mu(nn, 0) = mu_block;
    F4muA(nn, 0) = muA_block;
    F4se1(nn, 0) = se1_block;
    F4se2(nn, 0) = se2_block;
    
    F4sg2(nn, 0) = sg2_block;
    F4sG2(nn, 0) = sG2_block;
    
    F4Gammah(nn, 0) = bh2_block;
    F4gammah(nn, 0) = bh1_block;
    F4GinvsG2(nn, 0) = bh2_block / sG2_block;
    F4ginvsg2(nn, 0) = bh1_block / sg2_block;
    
    F4insGRinsG(nn, 0) = diagmat(1. / se2_block)*R_block*diagmat(1. / se2_block);
    F4insgRinsg(nn, 0) = diagmat(1. / se1_block)*R_block*diagmat(1. / se1_block);
    
    F4Rins(nn, 0) = R_block*diagmat(1 / se1_block);
    F4Rins2(nn, 0) = R_block*diagmat(1 / se2_block);
    
    F4RinsGmu(nn, 0) = R_block*diagmat(1 / se2_block)*mu_block;
    F4Rinsgmu(nn, 0) = R_block*diagmat(1 / se1_block)*mu_block;
    F4RinsGmuA(nn, 0) = R_block*diagmat(1. / se2_block)*muA_block;
    
    F4sGRsG(nn, 0) = diagmat(se2_block)*R_block*diagmat(se2_block);
    F4BB(nn, 0) = diagmat(se2_block)*R_block*diagmat(1. / sG2_block)*R_block*diagmat(se2_block);
  }
  
  cout << "Finish the Preprocessing Process in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
  
  int l = 0;
  
  cout << endl;
  cout << "Start gibbs sampling Process:" << endl;
  t1 = clock();
  
  
  vec beta0mean = zeros(nblocks, 1); // mean for each block.
  vec beta0sig = zeros(nblocks, 1); // sigma for each block.
  vec beta1mean = zeros(nblocks, 1); // mean for each block.
  vec beta1sig = zeros(nblocks, 1); // sigma for each block.
  vec Mu2 = zeros(nblocks, 1); // sgga2 for each block.
  vec MuA2 = zeros(nblocks, 1); // sgal2 for each block.
  
  for(int iter = 0; iter < (int_fast32_t)(maxIter + burnin); iter ++){
    double logW = log(w /( 1 - w));
    
    // ------------------------------------------------------------------
    // set parallel computation for gamma, alpha, eta
    paraBlock_GamAlpEta parobj_GamAlpEta(nblocks, F4se1, F4se2, F4gammah, F4Gammah, F4ginvsg2,
                                         F4GinvsG2,  F4mu, F4muA, F4Rinsgmu, F4RinsGmu, 
                                         F4RinsGmuA, F4insgRinsg, F4insGRinsG, F4Rins, F4Rins2, F4Rblock, 
                                         beta0sig, beta1sig, beta0mean, beta1mean, Eta,
                                         Mu2, MuA2, beta0, beta1, logW, sgga2, sgal2, conspar);
    
    const int n_thread = coreNum;
    std::vector<std::thread> threads(n_thread);
    for(int i_thread = 0; i_thread < n_thread; i_thread++){
      threads[i_thread] = std::thread(&paraBlock_GamAlpEta::update_by_thread_GamAlpEta, &parobj_GamAlpEta, i_thread);
    }
    
    
    for(int i = 0; i < n_thread; i++){
      threads[i].join();
    }
    
    // save the parallel result
    beta0sig =  parobj_GamAlpEta.beta0sig;
    beta1sig =  parobj_GamAlpEta.beta1sig;
    
    beta0mean =  parobj_GamAlpEta.beta0mean;
    beta1mean =  parobj_GamAlpEta.beta1mean;
    
    Mu2 = parobj_GamAlpEta.Mu2;
    MuA2 = parobj_GamAlpEta.MuA2;
    
    F4mu = parobj_GamAlpEta.F4mu;
    F4muA = parobj_GamAlpEta.F4muA;
    F4Rinsgmu = parobj_GamAlpEta.F4Rinsgmu;
    F4RinsGmu = parobj_GamAlpEta.F4RinsGmu;
    F4RinsGmuA = parobj_GamAlpEta.F4RinsGmuA;
    
    Eta = parobj_GamAlpEta.Eta;
    // ------------------------------------------------------------------
    
    
    // ------------------ //
    // update omega
    // ------------------ //
    Eta1 = 1 - Eta;
    double wpa1, wpa2;
    wpa1 = a + sum(Eta1);
    wpa2 = b + Eta1.n_elem - sum(Eta1);
    w = 1 - R::rbeta(wpa1, wpa2);
    // if(iter > 300){
    //   cout<<"sum(Eta)::"<<sum(Eta) << "wpa1:"<< wpa1 << "--wpa2:" << wpa2 << "--w:" << w << endl;
    //   
    // }
    // 
    // cout<<"iter:"<<iter << "--sum(Eta)::"<<sum(Eta) << "wpa1:"<< wpa1 << "--wpa2:" << wpa2 << "--w:" << w << endl;
    
    // ------------------ //
    // update beta0
    // ------------------ //
    double sig2b0, mub0;
    if(sum(Eta)==0){
      beta0 = 0;
    }else{
      sig2b0 = 1. / sum(beta0sig);
      mub0 = sum(beta0mean) * sig2b0;
      beta0 = mub0 + randn()*sqrt(sig2b0);
      // beta0 = mub0;
    }
    // cout<<"Mb0:"<< beta0 << "--sigb0:" << sqrt(sig2b0) << endl;
    // ------------------ //
    // update beta1
    // ------------------ //
    double sig2b1, mub1;
    if(sum(Eta)==(int)nblocks){
      beta1 = 0;
    }else{
      sig2b1 = 1. / sum(beta1sig);
      mub1 = sum(beta1mean) * sig2b1;
      beta1 = mub1 + randn()*sqrt(sig2b1);
      // beta1 = mub1;
      // cout <<"check error 2" << endl;
    }
    // cout<<"Mb1:"<< beta1 << "--sigb1:" << sqrt(sig2b1) << endl;
    // ------------------ //
    // update sgga2
    // ------------------ //
    double tagm, tbgm, taal, tbal;
    tagm = agm + p / 2;
    tbgm = sum(Mu2) / 2 + bgm;
    sgga2 =  1 / randg<double>(distr_param(tagm, 1/tbgm));
    // sgga2 = tbgm;
    // ------------------ //
    // update sgal2
    // ------------------ //
    taal = aal + p / 2;
    tbal = sum(MuA2) / 2 + bal;
    sgal2 =  1 / randg<double>(distr_param(taal, 1/tbal));
    // sgal2 = tbal;
    
    // cout<<"sgga2:" << sgga2 <<"--tagm--"<< tagm << endl;
    // cout<<"sgal2:" << sgal2 <<"--taal--"<< taal << endl;
    // cout<<"beta0--"<<beta0<<"--beta1--"<<beta1 << endl;
    // cout<<"--------------------------------------------"<<endl;
    // cout<<endl;
    if(iter >= (int)burnin){
      if((iter - burnin) % thin ==0){
        // cout << "Eta:" << sum(Eta)<< endl;
        Beta0res[l] = beta0;
        Beta1res[l] = beta1;
        Etares[l] = 1 - (float)sum(Eta)/(float)Eta.n_elem;
        // EtaAll.col(l) = Eta;
        Sgga2Res[l] = sgga2;
        Sgal2Res[l] = sgal2;
        l += 1;
      }
    }
    
    if(iter%500==0){
      cout << "the iteration number is:" << iter << endl;
    }
    
    
  }
  
  
  cout << "Finish the gibbs sampling Process in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
  ObjGibbs3 obj;
  // 
  
  obj.Eta = Eta1;
  obj.Etares = Etares;
  obj.Beta0res = Beta0res;
  obj.Beta1res = Beta1res;
  obj.Sgga2Res = Sgga2Res;
  obj.Sgal2Res = Sgal2Res;
  
  return obj;
}

// [[Rcpp::export]]
Rcpp::List MRCorr2Real(arma::ivec &bp, arma::ivec &chr, arma::uvec &avbIndex, arma::uvec &idx4panel, std::string &block_file,
                                     std::string stringname3, double ld_r2_thresh,
                                     arma::vec &bh1, arma::vec &bh2, arma::vec &se1, arma::vec &se2,
                                     double lam,  int coreNum,   SEXP opts = R_NilValue){
  
  
  
  List Rblockres =  Cal_blockR(bp, chr, avbIndex, idx4panel, block_file,
                               stringname3,  ld_r2_thresh, coreNum, lam);
  arma::field<mat> F4Rblock = Rblockres["F4Rblock"];
  arma::umat block_inf = Rblockres["block_inf"];
  uword nblocks = Rblockres["nblocks"];
  arma::ivec Indpid = Rblockres["Indpid"];
  
  
  Options_Gibbs3* lp_opt = NULL;
  if (!Rf_isNull(opts)){
    Rcpp::List opt(opts);
    lp_opt = new Options_Gibbs3(opt["agm"], opt["bgm"], opt["aal"], opt["bal"], opt["a"], opt["b"], opt["maxIter"], opt["thin"], opt["burnin"]);
  }
  if (Rf_isNull(opts)){
    lp_opt = new Options_Gibbs3();
  }
  
  const int conspar = 0;
  ObjGibbs3 obj = gibbs3group_blockparReal(bh1,  bh2, se1, se2, F4Rblock, block_inf, coreNum, conspar, lp_opt);
  
  // gibbs3group_blockparReal(arma::vec &gammah, arma::vec &Gammah, arma::vec &se1, arma::vec &se2, 
  //                          arma::field<mat> F4Rblock, arma::umat block_inf,  int coreNum, Options_Gibbs3* opts)
  
  List output = List::create(
    Rcpp::Named("nblocks") = Rcpp::wrap(nblocks),
    Rcpp::Named("Indpid") = Rcpp::wrap(Indpid),
    Rcpp::Named("Eta") = Rcpp::wrap(obj.Eta),
    Rcpp::Named("Etares") = Rcpp::wrap(obj.Etares),
    // Rcpp::Named("EtaAll") = Rcpp::wrap(obj.EtaAll),
    Rcpp::Named("Beta0res") = Rcpp::wrap(obj.Beta0res),
    Rcpp::Named("Beta1res") = Rcpp::wrap(obj.Beta1res),
    Rcpp::Named("Sgga2Res") = Rcpp::wrap(obj.Sgga2Res),
    Rcpp::Named("Sgal2Res") = Rcpp::wrap(obj.Sgal2Res)
  );
  return output;
}


mat cal_blockcor(arma::umat& X){
  uword size1 = X.n_rows;
  uword size2 = X.n_cols;
  vec meanX(size2);
  vec sqrtsum(size2);
  mat Xnorm = conv_to<mat>::from(X);
  for (int i = 0; i < (int)(size2); i++) { //calculate the mean of the vector and sqrt sum
    meanX[i] = sum(Xnorm.col(i))*1.0 / size1;
    Xnorm.col(i) -= meanX[i];
    vec v_i = Xnorm.col(i);
    mat pd = v_i.t() * v_i;
    sqrtsum[i] = sqrt(pd.at(0));
    if (sqrtsum[i] > 0){
      Xnorm.col(i) /= sqrtsum[i];
    }
  }
  arma::mat corr(size2, size2);
  arma::mat eyeI(size2, size2);
  eyeI.eye();
  mat cor_ij(1, 1);
  corr.eye();
  for (int i = 0; i < (int)(size2); i++){
    for (int j = i + 1; j < (int)(size2); j++) {
      cor_ij = ((Xnorm.col(i)).t() * (Xnorm.col(j)));
      double value = cor_ij.at(0);
      corr(i, j) = value;
      corr(j, i) = value;
    }
  }
  //
  // corr *= 0.9;
  // corr += 0.1*eyeI;
  
  return corr;
  
}

// [[Rcpp::export]]
mat Cal_block_SimR(umat block_inf, arma::umat &X, double lam){
  
  
  int p = X.n_cols;
  int nblocks = block_inf.n_rows;
  mat R = zeros(p, p);
  for(int j = 0; j < (int)(nblocks); j++){
    umat subX = X.cols(block_inf(j, 0), block_inf(j, 1));
    
    arma::mat corr0 = cal_blockcor(subX);
    arma::mat corr;
    int p1 = corr0.n_rows;
    arma::mat LAM = zeros(p1, p1);
    LAM.fill(lam);
    if(lam > 0.5)
    {
      corr = corr0;
      arma::mat eyeI(p1, p1);
      eyeI.eye();
      corr = corr0;
      corr *= lam;
      corr += (1 - lam)*eyeI;
    }
    else
    {
      pdsoftObj out = pdsoft(corr0, LAM);
      corr = out.theta;
    }
    
    R.submat(block_inf(j, 0), block_inf(j, 0), block_inf(j, 1), block_inf(j, 1)) = corr;
    
  }
  
  
  return R;
}


void fastLm(const arma::vec & y, const arma::mat & X, double &coefb, double &stdb) {
  
  int n = X.n_rows, k = X.n_cols;
  
  arma::colvec coef = arma::solve(X, y); 
  arma::colvec resid = y - X*coef; 
  
  double sig2 = arma::as_scalar(arma::trans(resid)*resid/(n-k));
  arma::colvec stderrest = 
    arma::sqrt(sig2 * arma::diagvec( arma::inv(arma::trans(X)*X)) );
  
  coefb = coef[1];
  stdb = stderrest[1];
}

// [[Rcpp::export]]
List fastSigLm(const arma::vec & y, const arma::mat & X) {
  
  // int n = X.n_rows, k = X.n_cols;
  int p = X.n_cols;int n = X.n_rows;
  arma::mat xx = zeros(p, 2);
  double coefb = 0;
  double stdb = 0;
  vec coef = zeros(p, 1);
  vec std = zeros(p, 1);
  
  for( int j = 0; j < p; j = j + 1 )
  {
    xx = join_rows(ones(n, 1), X.col(j));
    fastLm(y, xx, coefb, stdb);
    coef[j] = coefb;
    std[j] = stdb;
  }
  
  return List::create(Named("coef") = coef,
                      Named("std") = std);
}



