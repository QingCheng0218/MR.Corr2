#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>
#include <mutex>
#include "GibbsAlpGamEta_ptr.hpp"


using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

ObjGibbs3 MRcorrObj(arma::vec &gammah, arma::vec &Gammah, arma::vec &se1, arma::vec &se2, 
                    Options_Gibbs3* opts)
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
  vec Sgal2Res = ones(numsave, 1);
  vec Sgga2Res = ones(numsave, 1);
  imat EtaAll = ones<imat>(p, numsave);
  
  vec mu = 0.01*ones(p, 1);
  vec muA = 0.01*ones(p, 1);
  
  
  ivec Eta = zeros<ivec>(p, 1);
  ivec Eta1 = zeros<ivec>(p, 1);
  // ----------------------------------------------------------------------
  vec sG2 = se2%se2;
  vec sg2 = se1%se1;
  vec invsG2 = 1. / sG2;
  vec invsg2 = 1. / sg2;
  
  // cout <<"check error 1" << endl;
  vec GinvsG2 = Gammah / sG2;
  vec ginvsg2 = gammah / sg2;
  
  double b02 = beta0*beta0;
  double b12 = beta1*beta1;
  double logw = log(w / (1 - w));
  int l = 0;
  for(int iter = 0; iter < (int_fast32_t)(maxIter + burnin); iter ++){
    double invsgga2 = 1. / sgga2;
    double invsgal2 = 1. / sgal2;
    

    
    vec v21, v20 = zeros(p, 1);
    v21 = 1. / (invsg2 + b02*invsG2 + invsgga2);
    v20 = 1. / (invsg2 + b12*invsG2 + invsgga2);
   
    vec v2A = zeros(p, 1);
    v2A = 1./(invsG2 + invsgal2);
    

    vec Bm1, bm1, Bm1invsgal2, invBm1;
    vec Bm011, Bm022, Bm012;
    mat bm0;
    
    Bm1 = invsgga2 + invsg2 + b02*invsG2;
    bm1 = ginvsg2 + beta0*GinvsG2;
    Bm1invsgal2 = Bm1*invsgal2;
    invBm1 = 1. / Bm1;
    bm0 = join_rows(ginvsg2 + beta1*GinvsG2, GinvsG2);
    Bm011 = invsgga2 + invsg2 + b12*invsG2;
    Bm022 = invsgal2 + invsG2;
    Bm012 = beta1*invsG2;
    
    mat Bm0j = zeros(2, 2);
    mat invBm0j = zeros(2, 2);
    double tmp, lik1, lik0, prob0, prob;
    
    for(int j = 0; j < p; j++){
      // ------------------ //
      // update gamma
      // ------------------ //
      if(Eta[j]==1){
        double mu1;
        mu1 = (ginvsg2[j] + beta0*GinvsG2[j])*v21[j];
        mu[j] = mu1 + randn()*sqrt(v21[j]);
        // mu[j] = mu1;
      }else{
        double mu0;
        mu0 = (ginvsg2[j] + beta1*GinvsG2[j] - beta1*muA[j]*invsG2[j])*v20[j];
        mu[j] = mu0 + randn()*sqrt(v20[j]);
        // mu[j] = mu0;
      }
      
      // ------------------ //
      // update alpha
      // ------------------ //
      if(Eta[j]==1){
        muA[j] = randn()*sqrt(sgal2);
      }else{
        double muA0;
        muA0 = (GinvsG2[j] - beta1*invsG2[j]*mu[j])*v2A[j];
        muA[j] = muA0 + randn()*sqrt(v2A[j]);
        // muA[k] = muA0;
      }
      
      // ------------------ //
      // update eta
      // ------------------ //
      rowvec bm0j = bm0.row(j);
      // cout<<"check error 11 " <<bm0j<< endl;
      // cout << Bm0j << endl;
      Bm0j(0, 0) = Bm011[j];
      Bm0j(0, 1) = Bm012[j];
      Bm0j(1, 0) = Bm012[j];
      Bm0j(1, 1) = Bm022[j];
      // cout<<"check error 1 " << endl;
      invBm0j = inv_sympd(Bm0j);
      tmp = as_scalar(bm0j*invBm0j*bm0j.t());
      // 
      lik1 = 0.5*invBm1[j]*bm1[j]*bm1[j] - 0.5*log(Bm1invsgal2[j]);
      lik0 = 0.5*tmp - sum(log(diagvec(chol(Bm0j))));
      // cout<<"check error 2 " <<  "lik1:" << 0.5*invBm1[j]*bm1[j]*bm1[j]  << endl;
      prob0 = logw + lik1 - lik0;
      prob = 1. / ( 1 + exp(-prob0));
      Eta[j] = R::rbinom(1, prob);
    }
    // cout << "mu:" << mu.subvec(0, 4).t() << "sum:"<< sum(mu) << endl; 
    // cout << "muA:" << muA.subvec(0, 4).t() << "sum:"<< sum(muA) << endl; 
    // ------------------ //
    // update beta0
    // ------------------ //
    vec muEta, muEta1, muAEta1;
    muEta = mu%Eta;
    muEta1 = mu%(1 - Eta);
    muAEta1 = muA%(1 - Eta);
    
    double sig2b0, mub0;
    if(sum(Eta)==0){
      beta0 = 0;
    }else{
      sig2b0 = 1. / sum(muEta%invsG2%muEta);
      mub0 = sum(GinvsG2%muEta)*sig2b0;
      beta0 = mub0 + randn()*sqrt(sig2b0);
      // beta0 = mub0;
    }
    
    // ------------------ //
    // update beta1
    // ------------------ //
    double sig2b1, mub1;
    if(sum(Eta)==p){
      beta1 = 0;
    }else{
      sig2b1 = 1. / sum(muEta1%invsG2%muEta1);
      mub1 = (sum(GinvsG2%muEta1) - sum(muEta1%invsG2%muAEta1))*sig2b1;
      beta1 = mub1 + randn()*sqrt(sig2b1);
      // beta1 = mub1;
    }
    b02 = beta0*beta0;
    b12 = beta1*beta1;
    
    // cout << "sig2b0:"<< sig2b0 <<"beta0:" <<beta0 << "sig2b1:"<<sig2b1 << "beta1:"<< beta1 << endl;
    // ------------------ //
    // update sgga2
    // ------------------ //
    double tagm, tbgm, taal, tbal;
    tagm = agm + p / 2;
    tbgm = as_scalar(mu.t()*mu) / 2 + bgm;
    sgga2 =  1 / randg<double>(distr_param(tagm, 1/tbgm));
    // sgga2 = tbgm;
    // ------------------ //
    // update sgal2
    // ------------------ //
    taal = aal + p / 2;
    tbal = as_scalar(muA.t()*muA) / 2 + bal;
    sgal2 =  1 / randg<double>(distr_param(taal, 1/tbal));
    // sgal2 = tbal;
    
    // ------------------ //
    // update omega
    // ------------------ //
    // double wpa1, wpa2, w, logw;
    Eta1 = 1 - Eta;
    double wpa1, wpa2;
    wpa1 = a + sum(Eta1);
    wpa2 = b + Eta1.n_elem - sum(Eta1);
    w = 1 - R::rbeta(wpa1, wpa2);
    
    // w = 0.1;
    logw = log(w / (1 - w));
    
    
    
    if(iter >= (int)burnin){
      if((iter - burnin) % thin ==0){
        Etares[l] = 1 - (float)sum(Eta)/(float)Eta.n_elem;
        // cout << "Eta:" << sum(eta)<< endl;
        Sgga2Res[l] = sgga2;
        Sgal2Res[l] = sgal2;
        Beta0res[l] = beta0;
        Beta1res[l] = beta1;
        EtaAll.col(l) = Eta1;
        l += 1;
      }
    }
    
  }
  
  ObjGibbs3 obj;
  
  obj.Eta = Eta1;
  obj.Etares = Etares;
  obj.EtaAll = EtaAll;
  obj.Beta0res = Beta0res;
  obj.Beta1res = Beta1res;
  obj.Sgga2Res = Sgga2Res;
  obj.Sgal2Res = Sgal2Res;
  
  return obj;
}


//[[Rcpp::export]]
List MRcorr(arma::vec &gammah, arma::vec &Gammah, arma::vec &se1, arma::vec &se2,
                           SEXP opts = R_NilValue)
{
  Options_Gibbs3* lp_opt = NULL;
  if (!Rf_isNull(opts)){
    Rcpp::List opt(opts);
    lp_opt = new Options_Gibbs3(opt["agm"], opt["bgm"], opt["aal"], opt["bal"], opt["a"], opt["b"], opt["maxIter"], opt["thin"], opt["burnin"]);
  }
  if (Rf_isNull(opts)){
    lp_opt = new Options_Gibbs3();
  }
  
  ObjGibbs3 obj = MRcorrObj(gammah, Gammah, se1, se2, lp_opt);
  List output = List::create(
    Rcpp::Named("Eta") = Rcpp::wrap(obj.Eta),
    Rcpp::Named("Etares") = Rcpp::wrap(obj.Etares),
    Rcpp::Named("EtaAll") = Rcpp::wrap(obj.EtaAll),
    Rcpp::Named("Beta0res") = Rcpp::wrap(obj.Beta0res),
    Rcpp::Named("Beta1res") = Rcpp::wrap(obj.Beta1res),
    Rcpp::Named("Sgga2Res") = Rcpp::wrap(obj.Sgga2Res),
    Rcpp::Named("Sgal2Res") = Rcpp::wrap(obj.Sgal2Res)
    
  );
  return output;
}


// ObjGibbs3 MRcorrObj(arma::vec &gammah, arma::vec &Gammah, arma::vec &se1, arma::vec &se2, 
//                     Options_Gibbs3* opts)
// {
//   
//   // ----------------------------------------------------------------------
//   // check number of input arguments
//   double agm = opts -> agm;
//   double bgm = opts -> bgm;
//   double aal = opts -> aal;
//   double bal = opts -> bal;
//   double a = opts -> a;
//   double b = opts -> b;
//   uword maxIter = opts -> maxIter;
//   uword thin = opts -> thin;
//   uword burnin = opts -> burnin;
//   
//   // ----------------------------------------------------------------------
//   // initial values
//   int p = Gammah.n_elem;
//   double sgga2 = 1; double sgal2 = 1; double beta0 = 0.01; double w = 0.1;
//   double beta1 = 0.01; 
//   int numsave = maxIter / thin;
//   vec Beta0res = ones(numsave, 1);
//   vec Beta1res = ones(numsave, 1);
//   vec Etares = ones(numsave, 1);
//   vec Sgal2Res = ones(numsave, 1);
//   vec Sgga2Res = ones(numsave, 1);
//   imat EtaAll = ones<imat>(p, numsave);
//   
//   vec mu = 0.01*ones(p, 1);
//   vec muA = 0.01*ones(p, 1);
//   
//   
//   ivec Eta = zeros<ivec>(p, 1);
//   // ----------------------------------------------------------------------
//   vec sG2 = se2%se2;
//   vec sg2 = se1%se1;
//   vec invsG2 = 1. / sG2;
//   vec invsg2 = 1. / sg2;
//   
//   // cout <<"check error 1" << endl;
//   vec GinvsG2 = Gammah / sG2;
//   vec ginvsg2 = gammah / sg2;
//   
//   double b02 = beta0*beta0;
//   double b12 = beta1*beta1;
//   double logw = log(w / (1 - w));
//   int l = 0;
//   for(int iter = 0; iter < (int_fast32_t)(maxIter + burnin); iter ++){
//     double invsgga2 = 1. / sgga2;
//     double invsgal2 = 1. / sgal2;
//     
//     // ------------------ //
//     // update gamma
//     // ------------------ //
//     
//     vec v21, v20 = zeros(p, 1);
//     v21 = 1. / (invsg2 + b02*invsG2 + invsgga2);
//     v20 = 1. / (invsg2 + b12*invsG2 + invsgga2);
//     
//     for(int j = 0; j < p; j++){
//       if(Eta[j]==1){
//         double mu1;
//         mu1 = (ginvsg2[j] + beta0*GinvsG2[j])*v21[j];
//         mu[j] = mu1 + randn()*sqrt(v21[j]);
//         // mu[j] = mu1;
//       }else{
//         double mu0;
//         mu0 = (ginvsg2[j] + beta1*GinvsG2[j] - beta1*muA[j]*invsG2[j])*v20[j];
//         mu[j] = mu0 + randn()*sqrt(v20[j]);
//         // mu[j] = mu0;
//       }
//     }
//     // cout << "mu:" << mu.subvec(0, 4).t() << "sum:"<< sum(mu) << endl; 
//     // ------------------ //
//     // update alpha
//     // ------------------ //
//     vec v2A = zeros(p, 1);
//     v2A = 1./(invsG2 + invsgal2);
//     if(sum(Eta)==p){
//       muA = zeros(p, 1);
//     }else{
//       for(int k = 0; k < p; k++){
//         if(Eta[k]==1){
//           muA[k] = randn()*sqrt(sgal2);
//         }else{
//           double muA0;
//           muA0 = (GinvsG2[k] - beta1*invsG2[k]*mu[k])*v2A[k];
//           muA[k] = muA0 + randn()*sqrt(v2A[k]);
//           // muA[k] = muA0;
//         }
//       }
//     }
//     // cout << "muA:" << muA.subvec(0, 4).t() << "sum:"<< sum(muA) << endl; 
//     // ------------------ //
//     // update beta0
//     // ------------------ //
//     vec muEta, muEta1, muAEta1;
//     muEta = mu%Eta;
//     muEta1 = mu%(1 - Eta);
//     muAEta1 = muA%(1 - Eta);
//     
//     double sig2b0, mub0;
//     if(sum(Eta)==0){
//       beta0 = 0;
//     }else{
//       sig2b0 = 1. / sum(muEta%invsG2%muEta);
//       mub0 = sum(GinvsG2%muEta)*sig2b0;
//       beta0 = mub0 + randn()*sqrt(sig2b0);
//       // beta0 = mub0;
//     }
//     
//     // ------------------ //
//     // update beta1
//     // ------------------ //
//     double sig2b1, mub1;
//     if(sum(Eta)==p){
//       beta1 = 0;
//     }else{
//       sig2b1 = 1. / sum(muEta1%invsG2%muEta1);
//       mub1 = (sum(GinvsG2%muEta1) - sum(muEta1%invsG2%muAEta1))*sig2b1;
//       beta1 = mub1 + randn()*sqrt(sig2b1);
//       // beta1 = mub1;
//     }
//     b02 = beta0*beta0;
//     b12 = beta1*beta1;
//     
//     // cout << "sig2b0:"<< sig2b0 <<"beta0:" <<beta0 << "sig2b1:"<<sig2b1 << "beta1:"<< beta1 << endl;
//     // ------------------ //
//     // update sgga2
//     // ------------------ //
//     double tagm, tbgm, taal, tbal;
//     tagm = agm + p / 2;
//     tbgm = as_scalar(mu.t()*mu) / 2 + bgm;
//     sgga2 =  1 / randg<double>(distr_param(tagm, 1/tbgm));
//     // sgga2 = tbgm;
//     // ------------------ //
//     // update sgal2
//     // ------------------ //
//     taal = aal + p / 2;
//     tbal = as_scalar(muA.t()*muA) / 2 + bal;
//     sgal2 =  1 / randg<double>(distr_param(taal, 1/tbal));
//     // sgal2 = tbal;
//     
//     // ------------------ //
//     // update omega
//     // ------------------ //
//     // double wpa1, wpa2, w, logw;
//     double wpa1, wpa2;
//     wpa1 = a + sum(Eta);
//     wpa2 = b + p - sum(Eta);
//     // cout<<"wpa1:"<< wpa1 << "wpa2:" << wpa2 << endl;
//     w = R::rbeta(wpa1, wpa2);
//     // w = 0.1;
//     logw = log(w / (1 - w));
//     
//     
//     // ------------------ //
//     // update eta
//     // ------------------ //
//     vec Bm1, bm1, Bm1invsgal2, invBm1;
//     vec Bm011, Bm022, Bm012;
//     mat bm0;
//     
//     Bm1 = invsgga2 + invsg2 + b02*invsG2;
//     bm1 = ginvsg2 + beta0*GinvsG2;
//     Bm1invsgal2 = Bm1*invsgal2;
//     invBm1 = 1. / Bm1;
//     bm0 = join_rows(ginvsg2 + beta1*GinvsG2, GinvsG2);
//     Bm011 = invsgga2 + invsg2 + b12*invsG2;
//     Bm022 = invsgal2 + invsG2;
//     Bm012 = beta1*invsG2;
//     
//     mat Bm0j = zeros(2, 2);
//     mat invBm0j = zeros(2, 2);
//     double tmp, lik1, lik0, prob0, prob;
//     // vec psum = zeros(p, 1);
//     // vec bm0j;
//     // cout<<"check error 0 " << endl;
//     // cout <<bm0.row(0)<< endl;
//     for(int j = 0; j < p; j++){
//       rowvec bm0j = bm0.row(j);
//       // cout<<"check error 11 " <<bm0j<< endl;
//       // cout << Bm0j << endl;
//       Bm0j(0, 0) = Bm011[j];
//       Bm0j(0, 1) = Bm012[j];
//       Bm0j(1, 0) = Bm012[j];
//       Bm0j(1, 1) = Bm022[j];
//       // cout<<"check error 1 " << endl;
//       invBm0j = inv_sympd(Bm0j);
//       tmp = as_scalar(bm0j*invBm0j*bm0j.t());
//       // 
//       lik1 = 0.5*invBm1[j]*bm1[j]*bm1[j] - 0.5*log(Bm1invsgal2[j]);
//       lik0 = 0.5*tmp - sum(log(diagvec(chol(Bm0j))));
//       // cout<<"check error 2 " <<  "lik1:" << 0.5*invBm1[j]*bm1[j]*bm1[j]  << endl;
//       prob0 = logw + lik1 - lik0;
//       prob = 1. / ( 1 + exp(-prob0));
//       Eta[j] = R::rbinom(1, prob);
//       // cout<<"prob:" << prob0<< endl;
//     }
//     
//     if(iter >= (int)burnin){
//       if((iter - burnin) % thin ==0){
//         Etares[l] = (float)sum(Eta)/(float)Eta.n_elem;
//         // cout << "Eta:" << sum(eta)<< endl;
//         Sgga2Res[l] = sgga2;
//         Sgal2Res[l] = sgal2;
//         Beta0res[l] = beta0;
//         Beta1res[l] = beta1;
//         EtaAll.col(l) = Eta;
//         l += 1;
//       }
//     }
//     
//   }
//   
//   ObjGibbs3 obj;
//   
//   obj.Eta = Eta;
//   obj.Etares = Etares;
//   obj.EtaAll = EtaAll;
//   obj.Beta0res = Beta0res;
//   obj.Beta1res = Beta1res;
//   obj.Sgga2Res = Sgga2Res;
//   obj.Sgal2Res = Sgal2Res;
//   
//   return obj;
// }