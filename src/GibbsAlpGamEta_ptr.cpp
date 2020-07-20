#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>
#include <mutex>
#include "GibbsAlpGamEta_ptr.hpp"
#include <random>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]


void paraBlock_GamAlpEta::loop_by_block_gibbs_GamAlpEta(int i){
  
  int eta = Eta[i];
  vec se1 = F4se1(i, 0);
  vec se2 = F4se2(i, 0);
  mat Rins = F4Rins(i, 0);
  mat Rins2 = F4Rins2(i, 0);
  mat Rblock = F4Rblock(i, 0);
  mat insGRinsG = F4insGRinsG(i, 0);
  mat insgRinsg = F4insgRinsg(i, 0);
  vec ginvsg2 = F4ginvsg2(i, 0);
  vec GinvsG2 = F4GinvsG2(i, 0);
  vec mu = F4mu(i, 0);
  vec muA = F4muA(i, 0);
  
  
  vec Rinsgmu = F4Rinsgmu(i, 0);
  vec RinsGmu = F4RinsGmu(i, 0);
  vec RinsGmuA = F4RinsGmuA(i, 0);
  
  vec diaginsgRinsg = diagvec(insgRinsg);
  vec diaginsGRinsG = diagvec(insGRinsG);
  
  
  vec invse1 = 1. / se1;
  vec invse2 = 1. / se2;
  double invsgga2 = 1. / sgga2;
  double invsgal2 = 1. / sgal2;
  vec Rdiag = diagvec(Rblock);
  
  int p = Rblock.n_cols;
  
  // ------------------ //
  // update gamma
  // ------------------ //
  vec v21, v20 = zeros(p, 1);
  v21 = 1. / (diaginsgRinsg + pow(beta0, 2)*diaginsGRinsG + invsgga2);
  v20 = 1. / (diaginsgRinsg + pow(beta1, 2)*diaginsGRinsG + invsgga2);
  
  if(eta==1){
    for(int j = 0; j < p; j++){
      vec tmp1, tmp2;
      double RinSmujj1, RinSmujj2, mu1;
      tmp1 = Rinsgmu - Rins.col(j)*mu[j];
      tmp2 = RinsGmu - Rins2.col(j)*mu[j];
      RinSmujj1 = Rinsgmu[j] - Rdiag[j]*mu[j] / se1[j];
      RinSmujj2 = RinsGmu[j] - Rdiag[j]*mu[j] / se2[j];
      
      mu1 = (ginvsg2[j] + beta0*GinvsG2[j] - RinSmujj1 / se1[j] - beta0*beta0/se2[j] * RinSmujj2)*v21[j];
      mu[j] = mu1 + randn()*sqrt(v21[j]);
      // mu[j] = mu1;
      Rinsgmu = tmp1 + Rins.col(j)*mu[j];
      RinsGmu = tmp2 + Rins2.col(j)*mu[j];
      
      
    }
  }else{
    for(int j = 0; j < p; j++){
      vec tmp1, tmp2;
      double RinSmujj1, RinSmujj2, mu0;
      tmp1 = Rinsgmu - Rins.col(j)*mu[j];
      tmp2 = RinsGmu - Rins2.col(j)*mu[j];
      RinSmujj1 = Rinsgmu[j] - Rdiag[j]*mu[j] / se1[j];
      RinSmujj2 = RinsGmu[j] - Rdiag[j]*mu[j] / se2[j];
      
      mu0 = (ginvsg2[j] + beta1*GinvsG2[j] - RinSmujj1 / se1[j] - beta1*beta1/se2[j] * RinSmujj2 - beta1/se2[j]*RinsGmuA[j])*v20[j];
      
      mu[j] = mu0 + randn()*sqrt(v20[j]);
      // mu[j] = mu0;
      Rinsgmu = tmp1 + Rins.col(j)*mu[j];
      RinsGmu = tmp2 + Rins2.col(j)*mu[j];
      
    }
  }
  
  
  // ------------------ //
  // update alpha
  // ------------------ //
  vec v2A = zeros(p, 1);
  v2A = 1./(diaginsGRinsG + invsgal2);
  
  if(eta==1){
    for(int k = 0; k < p; k++){
      muA[k] = randn()*sqrt(sgal2);
      // muA[k] = 0.01;
    }
    RinsGmuA = Rblock*diagmat(invse2)*muA; // remember this update !!
    
  }else{
    for(int k = 0; k < p; k++){
      vec tmp3;
      double RinSmuAkk, muA0;
      tmp3 = RinsGmuA - Rins2.col(k)*muA[k];
      RinSmuAkk = RinsGmuA[k] - Rdiag[k]*muA[k]/se2[k];
      muA0 = (GinvsG2[k] - beta1/se2[k]*RinsGmu[k] - 1/se2[k]*RinSmuAkk)*v2A[k];
      muA[k] = muA0 + randn()*sqrt(v2A[k]);
      // muA[k] = muA0;
      RinsGmuA = tmp3 + Rins2.col(k)*muA[k];
    }
    
  }
  
  // ------------------ //
  // update eta
  // ------------------ //
  if(conspar==0){
    
    vec gammah = F4gammah(i, 0);
    vec Gammah = F4Gammah(i, 0);
    vec ginvsg2 = F4ginvsg2(i, 0);
    vec GinvsG2 = F4GinvsG2(i, 0);
    double term1, term2, prob0, prob;
    int lm = GinvsG2.n_elem;
    mat A1, A0, invA1, invA0;
    vec am1, am0;
    A1 = invsgga2*diagmat(ones(lm, 1)) + insgRinsg + beta0*beta0*insGRinsG;
    invA1 = inv_sympd(A1);
    A0 = join_cols(join_rows(invsgga2*diagmat(ones(lm, 1)) + insgRinsg + beta1*beta1*insGRinsG, beta1*insGRinsG),
                   join_rows(beta1*insGRinsG, invsgal2*diagmat(ones(lm, 1)) + insGRinsG));
    invA0 = inv_sympd(A0);
    am1 = ginvsg2 + beta0*GinvsG2;
    am0 = join_cols(ginvsg2 + beta1*GinvsG2, GinvsG2);
    term1 = 0.5*as_scalar(am1.t()*invA1*am1) - sum(log(diagvec(chol(A1))));
    term2 = - 0.5*as_scalar(am0.t()*invA0*am0) + sum(log(diagvec(chol(A0))))  + 0.5*lm*log(sgal2);
    prob0 = logW + term1 + term2;
    prob = 1. / (1 + exp(-prob0));
    eta = R::rbinom(1, prob);
    Eta[i] = eta;
    // 
    //     // The second caculation of likelihood(The original caculation)
    // 
    //     vec phim = join_cols(gammah, Gammah);
    //     mat sgRinvsg = F4sgRinvsg(i, 0);
    //     mat sGRinvsG = F4sGRinvsG(i, 0);
    //     mat sGRsG = F4sGRsG(i, 0);
    //     mat sgRsg = F4sgRsg(i, 0);
    //     mat Bgm2 = sgRinvsg*sgRinvsg.t();
    //     mat BgGm = sgRinvsg*sGRinvsG.t();
    //     mat BGgm = BgGm.t();
    //     mat BGm2 = sGRinvsG*sGRinvsG.t();
    // 
    //     mat Am = join_cols(join_rows(sgRsg, zeros(lm, lm)), join_rows(zeros(lm, lm), sGRsG));
    //     mat Am1 = sgga2*join_cols(join_rows(Bgm2, beta0*BgGm),
    //                               join_rows(beta0*BGgm, beta0*beta0*BGm2));
    //     mat Am0 = join_cols(join_rows(sgga2*Bgm2, beta1*sgga2*BgGm),
    //                         join_rows(beta1*sgga2*BGgm, (beta1*beta1*sgga2 + sgal2)*BGm2));
    // 
    //     mat Sig1 = Am + Am1;
    //     mat invSig1 = inv_sympd(Sig1);
    //     mat Sig0 = Am + Am0;
    //     mat invSig0 = inv_sympd(Sig0);
    // 
    //     double lik1, lik0, prob2;
    // 
    //     lik1 = - 0.5*as_scalar(phim.t()*invSig1*phim)- sum(log(diagvec(chol(Sig1))));
    // 
    //     lik0 = - 0.5*as_scalar(phim.t()*invSig0*phim) - sum(log(diagvec(chol(Sig0))));
    //     prob2 = logW + lik1 - lik0;
    // 
    //     cout <<"i:" << i << "prob0:" << prob0 << "prob2:" << prob2 <<endl;
    
    
    
  }
  
  // -----------------------------------------------------------------------
  F4mu(i, 0) = mu;
  F4muA(i, 0) = muA;
  F4Rinsgmu(i, 0) = Rinsgmu;
  F4RinsGmu(i, 0) = RinsGmu;
  F4RinsGmuA(i, 0) = RinsGmuA;
  // -----------------------------------------------------------------------
  // for beta0 and beta1 iteration
  vec muEta = mu*eta;
  vec muEta1 = mu*(1 - eta);
  vec muAEta1 = muA*(1 - eta);
  
  
  beta0sig[i] = as_scalar(muEta.t()*insGRinsG*muEta);
  beta0mean[i] = as_scalar(GinvsG2.t()*muEta);
  
  beta1sig[i] = as_scalar(muEta1.t()*insGRinsG*muEta1);
  beta1mean[i] = as_scalar(GinvsG2.t()*muEta1 - muEta1.t()*insGRinsG*muAEta1);
  
  // // -----------------------------------------------------------------------
  // for sgga2 and sgal2 iteration
  Mu2[i] = sum(mu%mu);
  MuA2[i] = sum(muA%muA);
  // -----------------------------------------------------
  se1.reset();
  se2.reset();
  Rins.reset();
  Rins2.reset();
  insGRinsG.reset();
  insgRinsg.reset();
  Rblock.reset();
  ginvsg2.reset();
  GinvsG2.reset();
  
}

std::mutex _mtx0;
int paraBlock_GamAlpEta::next_GamAlpEta(){
  std::lock_guard<std::mutex> lockGuard(_mtx0);
  if(current_idx >= (int)nblocks){
    return -1;
  }
  current_idx++;
  return current_idx-1;
}

void paraBlock_GamAlpEta::update_by_thread_GamAlpEta(int thread_id){
  while(true){
    int idx = next_GamAlpEta();
    if(idx == -1){
      break;
    }
    loop_by_block_gibbs_GamAlpEta(idx);
  }
}
