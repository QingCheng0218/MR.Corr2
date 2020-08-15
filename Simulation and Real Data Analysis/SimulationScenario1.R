
rm(list = ls());
library("mvtnorm");
library(MR.Corr2);
library("mvtnorm");
library("dplyr")
library(cause);
library(dplyr);
library(MR.LDP);
library(MRMix);
library(mr.raps);
library(gsmr);
# -------------------------------------------------------
rho <- 0.4; h2z <- 0.1; L <-50; b0 <- 0.1; rho_ag <- 0.2;
# -------------------------------------------------------
Alrate <- 0.1;
h2y = 0.1;
lam = 0.055;
Garate = 1;
M = 10;
p = M*L
n1 <- 20000;
n2 <- 20000;
n3 <- 500;
m <- p;
coreNum = 20;

filename <- paste0("r", 10*rho,  "hz", h2z*100, "hy", h2y*10,
                   "L", L, "b", b0,  "AG",rho_ag,
                   "Ga", Garate, "Al", Alrate, "cause.Rdata");

if(b0!=0){nrep = 100;}else{nrep = 1000}
bhat = pvalue = matrix(0, nrow = nrep, ncol = 6); # save IND Result for other methods.
bhat1 = pvalue1 = matrix(0, nrow = nrep, ncol = 4); # save Noprune Result for other methods.
nblocks = L;
block_inf <- cbind(seq(1, p, M), seq(M, p, M));
block_inf1 <- block_inf - 1;
p = M*L;
maf = runif(p,0.05,0.5);
x = genRawGeno(maf, L, M, rho, n1 + n2 + n3);
# -------------------------------------------------------------------
x1 = x[1:n1,];
x2 = x[(n1+1):(n1+n2),];
x12 = x[1:(n1+n2),];
x3 = x[(n1+n2+1):(n1+n2+n3),];
R = Cal_block_SimR(block_inf1, x3, lam)
#-------------------------------------------------------------------------#
# functions used in cause
ld_prune_cormat <- function(R, snp_names, p_vals,  p_val_thresh, r2_thresh){
  stopifnot(nrow(R) == length(snp_names))
  stopifnot(length(snp_names) == length(p_vals))
  
  ix <- which(p_vals < p_val_thresh)
  
  if(length(ix) == 0) return(c())
  
  snp_names <- snp_names[ix]
  R <- R[ix, ix, drop=FALSE]
  p_vals <- p_vals[ix]
  o_c <- order(p_vals, decreasing=FALSE)
  keep <- c()
  while(length(o_c) > 0){
    keep <- c(keep, snp_names[o_c[1]])
    myld <- R[o_c, o_c, drop=FALSE]
    remove_ix <- which(myld[,1]^2 > r2_thresh)
    o_c <- o_c[-remove_ix]
  }
  return(keep)
}
new_cause_data1 <- function(x = data.frame()){
  stopifnot(inherits(x, "data.frame"))
  # x <- validate_cause_data(x)
  #stopifnot(all(c("snp", "beta_hat_1", "seb1", "beta_hat_2", "seb2") %in% names(x)))
  structure(x, class = c("cause_data", "data.frame"))
}
cause_sims <- function(dat, param_ests, sigma_g, qalpha=1, qbeta=10,
                       no_ld = FALSE,
                       thresh = 1e-3){
  
  if(no_ld) dat <- process_dat_nold(dat)
  
  vars <- filter(dat, ld_prune == TRUE & p_value < thresh) %>% with(., snp)
  X <- new_cause_data(dat)
  
  #get sigma_g from data
  if(missing(sigma_g)) sigma_g <- cause:::eta_gamma_prior(X, vars)
  if(is.na(sigma_g)) sigma_g <- cause:::eta_gamma_prior(X, vars)
  
  res <- cause::cause(X=X, variants = vars, param_ests = param_ests, sigma_g = sigma_g,
                      qalpha = qalpha, qbeta = qbeta, force=TRUE)
  return(res)
}

#this is a helper function for switching from with ld to no ld data
process_dat_nold <- function(dat){
  dat <- dat %>%  select(-beta_hat_1, -beta_hat_2, -p_value, -ld_prune) %>%
    rename(beta_hat_1 = beta_hat_1_nold,
           beta_hat_2 = beta_hat_2_nold,
           p_value = p_value_nold) %>%
    mutate(ld_prune = TRUE)
  return(dat)
}


#-------------------------------------------------------------------------#

for(irep in 1:nrep){
  
  Stm <- proc.time();
  q = 50
  u = matrix(rnorm( (n1+n2) * q),ncol=q);
  
  # ------------------------------------------------------------------------
  sigma2g = 1;
  sigma2a = 1;
  
  S_ag = matrix(c(sigma2a, rho_ag, rho_ag, sigma2g), nrow=2);
  AG = rmvnorm(p, mean=rep(0, 2), sigma = S_ag,method="chol")
  alpha = AG[,1]; gamma = AG[,2];
  
  
  if(Garate!=1){
    gano = floor(p*(1 - Garate));
    indxGA = sample(1:p,gano);
    gamma[indxGA] = 0;
  }
  
  # ------------------------------------------------------------------------
  Su = matrix(c(1,0.8,0.8,1),nrow=2)
  bu = rmvnorm(q,mean=rep(0,2), sigma = Su,method="chol")
  by = bu[,1]; bz = bu[,2];
  uby = u%*%by; ubz = u%*%bz;
  uby = uby/sqrt(as.numeric(var(uby)/0.6));
  ubz = ubz/sqrt(as.numeric(var(ubz)/0.2));
  
  x12g = x12%*%gamma;
  
  if(b0!=0){
    h2ga = (h2y *( 1 + b0^2))/(b0^2 * (1 - h2y));
    gamma0 = gamma/sqrt(as.numeric(var(x12g)/h2ga));
    x12g = x12%*%gamma0;
  }
  
  
  yall = x12g + uby + rnorm(n1+n2)*as.numeric(sqrt(1-var(uby)));
  # ------------------------------------------------------------------------
  # The direct effects on Z
  h2yb = var(b0*yall);
  h2al = (h2z + h2z*h2yb)/(1 - h2z);
  
  if(h2z==0){
    alpha = rep(0, p);
    x12a = x12%*%alpha;
  }else{
    if(Alrate!=1){
      alno = floor(L*Alrate);
      ind = sample(1:L,alno);
      if(length(ind)==1){
        indxAL = block_inf[ind, 1]:block_inf[ind, 2]
        alpha[-indxAL] = 0;
      }else{
        indxAL = NULL;
        for(i in 1:length(ind)){
          tmp = block_inf[ind[i], 1]:block_inf[ind[i], 2]
          indxAL = append(indxAL, tmp);
        }
        alpha[-indxAL] = 0;
        x12a = x12%*%alpha;
        alpha0 = alpha/sqrt(as.numeric(var(x12a)/(h2al)));
        x12a = x12%*%alpha0;
      }
    }
  }
  
  # ------------------------------------------------------------------------
  resz = ubz + rnorm(n1+n2)*as.numeric(sqrt(1-var(ubz)));
  zall = b0*yall  + x12a +  resz;
  # H2a.res[irep] <- var(x12a)/var(zall);
  # H2g.res[irep] <- var(b0*x12g)/var(zall);
  
  y = yall[1:n1];
  z = zall[(n1+1):(n1+n2)];
  
  x1 = x12[1:n1, ];
  x2 = x12[(n1+1):(n1+n2), ]
  
  
  # create summary statistics
  gammaSS = fastSigLm(y, x1);
  GammaSS = fastSigLm(z, x2);
  gammah = gammaSS$coef;
  se1 = gammaSS$std;
  Gammah = GammaSS$coef;
  se2 = GammaSS$std;
  #-------------------------------------------------------------------------#
  # MR-Corr2
  # -----------------------------------------------------
  maxIter = 4000;
  burnin = 1000;
  a1 = 1;
  b1 = L
  opt = list(agm = 0.001, bgm = 0.001, aal = 0.001, bal = 0.001,
             a = a1, b = b1, maxIter = maxIter, thin = 10, burnin = burnin);
  
  
  SimRes = MRCorr2Sim(gammah, Gammah, se1, se2, R,  block_inf1,coreNum, opt);
  
  # EstEta[, irep] = SimRes$Eta;
  beta_hat = mean(SimRes$Beta0res);
  beta_se = sd(SimRes$Beta0res);
  bhat[irep, 1]  = beta_hat;
  pvalue[irep, 1] = 2*(1 - pnorm(abs(beta_hat/beta_se)));
  #-------------------------------------------------------------------------#
  # MR-LDP
  # -----------------------------------------------------
  # initial value.
  beta0 <- 0;
  sgga2 <- 1;
  agm <- bgm <- aal <- bal <- 0.001;
  gamma <- rep(0.01, p);
  alpha <- rep(0.01, p);
  sgga2 <- 1;
  sgal2 <- 1;
  nblocks = L;
  # block_inf <- cbind(seq(0, p-M, M), seq(M-1, p, M));
  epsStopLogLik = 1e-7;
  maxIter = 4000;
  
  
  pxresult2 = MRLDP_SimPXvb(gammah, Gammah, se1, se2, gamma, alpha, beta0, sgga2, sgal2, R,
                            0, epsStopLogLik, maxIter, model = 2);
  tstat20 = pxresult2$tstat;
  bhat[irep, 2] = pxresult2$beta0;
  
  if(b0==0){
    pxresult21 = MRLDP_SimPXvb(gammah, Gammah, se1, se2, gamma, alpha,  beta0, sgga2, sgal2, R,
                               1, epsStopLogLik, maxIter, model = 2);
    tstat21 = pxresult21$tstat;
    tstat2 = 2*(tstat20 - tstat21);
    pval2 = pchisq(tstat2, 1, lower.tail = F);
    pvalue[irep, 2] = pval2;
  }
  
 
  # -------------------------------------------------------------------------------------
  # CAUSE method:
  ld_prune_pval_thresh = 1e-3;
  r2_thresh = 0.1;
  
  snp = 1:p;
  p_value = 2*pnorm(-abs(gammah/se1));
  keep <- sapply(seq(L), function(i){
    strt <- block_inf[i, 1]
    stp <- block_inf[i, 2]
    # print(strt);
    # print(stp);
    R1 = R[strt:stp, strt:stp];
    ld_prune_cormat(R1, snp[strt:stp], p_value[strt:stp],  ld_prune_pval_thresh, r2_thresh)
    # print(ld_prune_cormat(R, df$snp[strt:stp], df$p_value[strt:stp],  ld_prune_pval_thresh, r2_thresh))
  }) %>% unlist()
  
  ld_prune = case_when(!snp %in% keep ~ FALSE, TRUE ~ TRUE);
  
  df <- data.frame(snp = snp,
                   beta_hat_1 = as.vector(gammah),
                   beta_hat_2 = as.vector(Gammah),
                   seb1 = as.vector(se1),
                   seb2 = as.vector(se2),
                   ld_prune = ld_prune,
                   p_value = p_value);
  
  X1 = data.frame(snp = snp,
                  beta_hat_1 = as.vector(gammah),
                  beta_hat_2 = as.vector(Gammah),
                  seb1 = as.vector(se1),
                  seb2 = as.vector(se2))
  
  X1 = new_cause_data1(X1);
  param = est_cause_params(X1, X1$snp, null_wt = 10, max_candidates = Inf);
  cause_res <- cause_sims(df, param, no_ld = FALSE)
  qs = summary(cause_res)$quants
  bhat[irep, 3] = qs[[2]][1,1];
  pvalue[irep, 3] <- pnorm(summary(cause_res)$z);
  
  # --------------------------------------------
  # gsmr
  # --------------------------------------------
  p_val_thresh  = 5e-8
  id4gsmr = which(p_value < p_val_thresh & ld_prune=="TRUE");
  Rgsmr = R[id4gsmr, id4gsmr];
  Rgsmr <- data.frame(Rgsmr)
  names(Rgsmr) <- seq_along(id4gsmr)
  rownames(Rgsmr) <- seq_along(id4gsmr)
  dat <- df[id4gsmr,] %>%
    mutate(p_value2 = 2*pnorm(-abs(beta_hat_2/seb2)))
  
  resgsmr <-  try(with(dat, gsmr(beta_hat_1, seb1, p_value,
                                 beta_hat_2, seb2, bzy_pval = p_value2,
                                 ldrho=Rgsmr, snpid=seq_along(id4gsmr),
                                 n_ref = 1, nsnps_thresh=1, gwas_thresh=p_val_thresh)))
  
  if(class(resgsmr)!="try-error") {
    bhat[irep, 4] <- resgsmr$bxy;
    pvalue[irep, 4] <- resgsmr$bxy_pval;
  }
  # --------------------------------------------
  # MRMix:
  # --------------------------------------------
  
  estMRMix = MRMix(gammah[id4gsmr], Gammah[id4gsmr], se1[id4gsmr], se2[id4gsmr], theta_temp_vec = seq(-0.1,0.2,by=0.001));
  
  bhat[irep, 5] = estMRMix$theta;
  pvalue[irep, 5] = estMRMix$pvalue_theta;

  # --------------------------------------------
  # RAPS:
  # --------------------------------------------
  
  id4ld = which(ld_prune == "TRUE")
  # # The RAPs method over.dispersion = TRUE for systematic pleiotropy.
  if(h2z==0){
    raps = mr.raps(gammah[id4ld], Gammah[id4ld], se1[id4ld], se2[id4ld]);
  }else if(h2z!=0){
    raps = mr.raps(gammah[id4ld], Gammah[id4ld], se1[id4ld], se2[id4ld], over.dispersion = TRUE);
  }
  
  bhat[irep, 6] <- raps$beta.hat;
  pvalue[irep, 6] = raps$beta.p.value;
  # -------------------------------------------
  # Noprune Result for CASUE, GSMR, MRMix, RAPS
  # -------------------------------------------
  # CAUSE method:
  
  snp = 1:p;
  p_value = 2*pnorm(-abs(gammah/se1));
  ld_prune = rep("TRUE", p)
  
  df <- data.frame(snp = snp,
                   beta_hat_1 = as.vector(gammah),
                   beta_hat_2 = as.vector(Gammah),
                   seb1 = as.vector(se1),
                   seb2 = as.vector(se2),
                   ld_prune = ld_prune,
                   p_value = p_value);
  X1 = data.frame(snp = snp,
                  beta_hat_1 = as.vector(gammah),
                  beta_hat_2 = as.vector(Gammah),
                  seb1 = as.vector(se1),
                  seb2 = as.vector(se2))
  
  X1 = new_cause_data1(X1);
  param = est_cause_params(X1, X1$snp, null_wt = 10, max_candidates = Inf);
  cause_res <- cause_sims(df, param, no_ld = FALSE, thresh = 1)
  qs = summary(cause_res)$quants
  bhat1[irep, 1] = qs[[2]][1,1];
  pvalue1[irep, 1] <- pnorm(summary(cause_res)$z);
  
  # --------------------------------------------
  # gsmr
  # --------------------------------------------
  # p_val_thresh  = 5e-8
  # id4gsmr = which(p_value < p_val_thresh & ld_prune=="TRUE");
  id4gsmr = 1:p;
  Rgsmr = R[id4gsmr, id4gsmr];
  Rgsmr <- data.frame(Rgsmr)
  names(Rgsmr) <- seq_along(id4gsmr)
  rownames(Rgsmr) <- seq_along(id4gsmr)
  dat <- df[id4gsmr,] %>%
    mutate(p_value2 = 2*pnorm(-abs(beta_hat_2/seb2)))
  
  resgsmr <-  try(with(dat, gsmr(beta_hat_1, seb1, p_value,
                                 beta_hat_2, seb2, bzy_pval = p_value2,
                                 ldrho=Rgsmr, snpid=seq_along(id4gsmr),
                                 n_ref = 1, nsnps_thresh=1, gwas_thresh=1)))
  
  if(class(resgsmr)!="try-error") {
    bhat1[irep, 2] <- resgsmr$bxy;
    pvalue1[irep, 2] <- resgsmr$bxy_pval;
  }
  
  # --------------------------------------------
  # MRMix:
  # --------------------------------------------
  
  estMRMix = MRMix(gammah[id4gsmr], Gammah[id4gsmr], se1[id4gsmr], se2[id4gsmr], theta_temp_vec = seq(-0.1,0.2,by=0.001));
  
  bhat1[irep, 3] = estMRMix$theta;
  pvalue1[irep, 3] = estMRMix$pvalue_theta
  # --------------------------------------------
  # RAPS:
  # --------------------------------------------
  
  # # The RAPs method over.dispersion = TRUE for systematic pleiotropy.
  if(h2z==0){
    raps = mr.raps(gammah, Gammah, se1, se2);
  }else if(h2z!=0){
    raps = mr.raps(gammah, Gammah, se1, se2, over.dispersion = TRUE);
  }
  
  bhat1[irep, 4] <- raps$beta.hat;
  pvalue1[irep, 4] = raps$beta.p.value;
  
  
  
  
  save(bhat, pvalue, file = filename);
  
}
