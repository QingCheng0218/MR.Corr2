rm(list = ls());
library(MR.LDP);
library(stringr);
library(mr.raps);
library(cause);
library(gsmr);
library(MRMix);
library(MR.Corr2);


lam = 0.9;
coreNum = 24;
maxIter = 5000;
burnin = 5000;
saveInd = maxIter/thin;
thin = 10;
ld_r2_thresh = 0.001;
stringname3 = "all_chr_1000G";
block_file = "fourier_ls-all.bed";

filename = paste0("HheightLam", lam*100, "bu", burnin, "max", maxIter, "ld_r2_thresh", ld_r2_thresh, ".Rdata");


filescreen <- "UKbiobank_height_used.txt";
fileexposure <- "GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N_used.txt";
fileoutcome <- "GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_WOMEN_N_used.txt";

# ------------------------------------------------------------------------------


PVA = sort(c(10^seq(-4, -7, -1), 5*10^seq(-4, -8, -1)));
BRES = list();
Parares = matrix(0, nrow = 7, ncol = length(PVA));
Bhat0Res = Bhat1Res = matrix(0, nrow = saveInd, ncol = length(PVA));
Etasave = matrix(0, nrow = saveInd, ncol = length(PVA));
# ---------------------------------------------------------------------
pva_cutoff = 2;
scrres = matchscreen(filescreen, fileexposure, fileoutcome, stringname3, pva_cutoff);

bh1 = as.numeric(scrres$bh1);
bh2 = as.numeric(scrres$bh2);
s12 = as.numeric(scrres$s12);
s22 = as.numeric(scrres$s22);
chr = as.numeric(scrres$chr);
bp = scrres$bp;
rsname = scrres$rsname
avbIndex = scrres$idxin;
idx4panel = scrres$idx4panel;

idx4panel_in = idx4panel;
screendata = read.table(filescreen, header = T);

# ------------------------------------------------------------------------------
# CAUSE estimate parameters
print(Sys.time())
cat("Start to estimate parameters:", "\n");
X1 = read.table(fileexposure, header = T);
X2 = read.table(fileoutcome, header = T);
X <- gwas_merge(X1, X2, snp_name_cols = c("SNP", "SNP"),
                beta_hat_cols = c("beta", "beta"),
                se_cols = c("se", "se"),
                A1_cols = c("A1", "A1"),
                A2_cols = c("A2", "A2"));

size1 = min(dim(X)[1], 1000000);
varlist <- with(X, sample(snp, size=size1, replace=FALSE))
params <- est_cause_params(X, varlist);

new_cause_data1 <- function(x = data.frame()){
  stopifnot(inherits(x, "data.frame"))
  # x <- validate_cause_data(x)
  #stopifnot(all(c("snp", "beta_hat_1", "seb1", "beta_hat_2", "seb2") %in% names(x)))
  structure(x, class = c("cause_data", "data.frame"))
};
# ------------------------------------------------------------------------------
# Select instrumental SNPs using screening data under different thresholds
idx = match(rsname, screendata$SNP);
Xdata0 = screendata[idx, ];

for(irep in 1:length(PVA)){
  BHres = matrix(0, nrow = 7, ncol = 4);
  Stm <- proc.time();
  
  pvalue = PVA[irep];
  idx1 = which(Xdata0$pvalue < pvalue);
  Xdata = Xdata0[idx1, ];
  idx0 = match(Xdata$SNP, rsname);
  # ------------------------------------------------------------------------------
  bh1_in = bh1[idx0];
  bh2_in = bh2[idx0];
  s12_in = s12[idx0];
  s22_in = s22[idx0];
  bp_in = bp[idx0];
  chr_in = chr[idx0];
  avbIndex_in = avbIndex[idx0];
  rsname_in = rsname[idx0];
  
  p = length(bh1_in);
  # ------------------------------------------------------------------------------
  #####---- MR-Corr2----#####
  # --------------------------
  # block information
  infres = Cal_blockinf(bp_in, chr_in, block_file);
  nblocks = infres$nblocks;
  
  opt = list(agm = 0.001, bgm = 0.001, aal = 0.001, bal = 0.001,
             a = 1, b = nblocks, maxIter = maxIter, thin = thin, burnin = burnin);
  
  result0 = MRCorr2Real(bp_in, chr_in, avbIndex_in-1, idx4panel_in, block_file, stringname3,
                        ld_r2_thresh, bh1_in, bh2_in, s12_in, s22_in, lam, coreNum, opt)
  
  Bhat0 = result0$Beta0res;
  Bhat1 = result0$Beta1res;
  Etasave[, irep] = result0$Etares;
  
  bhat1 = mean(Bhat0);
  sd1 = sd(Bhat0);
  tstat1 = bhat1 / sd1;
  pval1 = 2*(1 - pnorm(abs(tstat1)))
  BHres[1, ] = c(bhat1,  bhat1-1.96*sd1, bhat1+1.96*sd1, pval1);
  
  Etarate0 = length(which(result0$Eta==1))/result0$nblocks;
  Eta = result0$Etares;
  NB = result0$nblocks;
  Indpid = result0$Indpid + 1;
  p1 = length(Indpid);
  
  sgga2hat = median(result0$Sgga2Res);
  sgal2hat = median(result0$Sgal2Res);
  para = c(p, NB, p1, Etarate0, sgga2hat, sgal2hat);
  Parares[, irep] = para;
  Bhat0Res[, irep] = Bhat0;
  Bhat1Res[, irep] = Bhat1;
  # ------------------------------------------------------------------------------
  #####---- MR-LDP----#####
  # --------------------------
  print(Sys.time())
  cat("Start MR-LDP procedure:", "\n");
  
  alpha = rep(0.01, p);
  gamma = rep(0.01, p);
  sgga2 =  0.01;
  sgal2 =  0.01;
  beta0 = 0;
  maxIter1 = 10000
  coreNum = 24;
  epsStopLogLik = 1e-7;
  lam1 = 0; #for lambda = 0.9 shrinkage.
  RealMRLDP_Hb = MRLDP_RealPXvb_block(bp_in, chr_in, avbIndex_in-1, idx4panel_in, block_file,
                                      stringname3, bh1_in, bh2_in, s12_in, s22_in,
                                      gamma, alpha, beta0, sgga2, sgal2, coreNum,
                                      lam1, 0, epsStopLogLik, maxIter1, model = 2);
  
  
  RealMRLDP_H0 = MRLDP_RealPXvb_block(bp_in, chr_in, avbIndex_in-1,idx4panel_in,  block_file,
                                      stringname3, bh1_in, bh2_in, s12_in, s22_in,
                                      gamma, alpha, beta0, sgga2, sgal2, coreNum,
                                      lam1, 1, epsStopLogLik, maxIter1, model = 2);
  beta0_MRLDP = RealMRLDP_Hb$beta0;
  Tstat_MRLDP = 2*(RealMRLDP_Hb$tstat - RealMRLDP_H0$tstat);
  MRLDP_se = abs(RealMRLDP_Hb$beta0/sqrt(Tstat_MRLDP));
  pval2 = pchisq(Tstat_MRLDP, 1, lower.tail = F);
  BHres[2, ] = c(beta0_MRLDP, beta0_MRLDP-1.96*MRLDP_se, beta0_MRLDP+1.96*MRLDP_se, pval2);
  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  #####---- CAUSE----#####
  # --------------------------
  Rres = Cal_block_Rmatrix(bp_in, chr_in, avbIndex_in-1,idx4panel_in, block_file,
                           stringname3, ld_r2_thresh, coreNum, lam);
  IndR = Rres$IndR;
  Indpid = Rres$Indpid + 1;
  
  top_vars = rsname_in[Indpid];
  XX = data.frame(
    snp = rsname_in,
    beta_hat_1 = bh1_in,
    seb1 = s12_in,
    beta_hat_2 = bh2_in,
    seb2 = s22_in,
    A1 = rep("A", length(bh1_in)),
    A2 = rep("G", length(bh1_in))
  )
  XX1 = new_cause_data1(XX);
  
  resCAUSE <- cause(X=XX1, variants = top_vars, param_ests = params);
  
  qs = summary(resCAUSE)$quants
  
  BHres[3, ] = c(qs[[2]][1,1], qs[[2]][2,1], qs[[2]][3,1], summary(resCAUSE)$p);
  # ------------------------------------------------------------------------------
  #####---- GSMR----#####
  # --------------------------
  print(Sys.time())
  cat("Start GSMR procedure:", "\n");
  
  
  dat = data.frame(
    snp = rsname_in[Indpid],
    beta_hat_1 = bh1_in[Indpid],
    seb1 = s12_in[Indpid],
    beta_hat_2 = bh2_in[Indpid],
    seb2 = s22_in[Indpid],
    p_value = 2*pnorm(-abs(bh1_in[Indpid]/s12_in[Indpid])),
    p_value2 = 2*pnorm(-abs(bh2_in[Indpid]/s22_in[Indpid]))
  )
  idxsave = which(dat$beta_hat_1!=0);
  
  if(length(idxsave)!=length(dat$beta_hat_1)){
    cat("irep:", irep, "with threshold:", PVA[irep], "has 0 bh1!", "\n" );
  }
  dat1 = data.frame(
    snp = dat$snp[idxsave],
    beta_hat_1 = dat$beta_hat_1[idxsave],
    beta_hat_2 = dat$beta_hat_2[idxsave],
    seb1 = dat$seb1[idxsave],
    seb2 = dat$seb2[idxsave],
    p_value = dat$p_value[idxsave],
    p_value2 = dat$p_value2[idxsave]
  )
  p_val_thresh  = 1;
  IndR <- data.frame(IndR[idxsave, idxsave]);
  colnames(IndR) = rownames(IndR) = 1:nrow(IndR);
  resgsmr <-  try(with(dat1, gsmr(beta_hat_1, seb1, p_value,
                                  beta_hat_2, seb2, bzy_pval = p_value2,
                                  ldrho=IndR, snpid=1:nrow(IndR),
                                  n_ref = 1, nsnps_thresh=1, gwas_thresh=p_val_thresh)))
  
  if(class(resgsmr)!="try-error") {
    BHres[4, ] = c(resgsmr$bxy, resgsmr$bxy - 1.96*resgsmr$bxy_se, resgsmr$bxy + 1.96*resgsmr$bxy_se, resgsmr$bxy_pval)
  }
  # ------------------------------------------------------------------------------
  #####---- MRMix----#####
  # --------------------------
  print(Sys.time())
  cat("Start MRMix procedure:", "\n");
  
  # theta_temp_vec = seq(-0.1,0.2,by=0.001)
  estMRMix = MRMix(dat$beta_hat_1, dat$beta_hat_2, dat$seb1, dat$seb2);
  BHres[5, ] = c(estMRMix$theta, estMRMix$theta - 1.96*estMRMix$SE_theta,
                 estMRMix$theta + 1.96*estMRMix$SE_theta, estMRMix$pvalue_theta);
  # ------------------------------------------------------------------------------
  #####---- RAPS----#####
  # --------------------------
  print(Sys.time())
  cat("Start RAPS procedure:", "\n");
  # # The RAPs method over.dispersion = TRUE for systematic pleiotropy.
  # raps = mr.raps(dat$beta_hat_1, dat$beta_hat_2, dat$seb1, dat$seb2, over.dispersion = TRUE);
  # raps <-  try(mr.raps(dat$beta_hat_1, dat$beta_hat_2, dat$seb1, dat$seb2, over.dispersion = TRUE));
  raps <-  try(mr.raps(dat1$beta_hat_1, dat1$beta_hat_2, dat1$seb1, dat1$seb2, over.dispersion = TRUE));
  
  if(class(raps)!="try-error") {
    BHres[6, ] = c(raps$beta.hat, raps$beta.hat - 1.96*raps$beta.se, raps$beta.hat + 1.96*raps$beta.se, raps$beta.p.value)
    
  }
  
  # ------------------------------------------------------------------------------
  #####---- MR-Corr----#####
  # --------------------------
  p1 = length(Indpid);
  opt = list(agm = 0.001, bgm = 0.001, aal = 0.001, bal = 0.001,
             a = 1, b = p1, maxIter = maxIter, thin = 10, burnin = burnin);
  
  resgibbs = MRcorr(dat$beta_hat_1, dat$beta_hat_2, dat$seb1, dat$seb2, opt);
  
  Bhat0 = resgibbs$Beta0res;
  
  sgga2hat = median(resgibbs$Sgga2Res);
  sgal2hat = median(resgibbs$Sgal2Res);
  
  bhat1 = mean(Bhat0);
  sd1 = sd(Bhat0);
  tstat1 = bhat1 / sd1;
  pval1 = 2*(1 - pnorm(abs(tstat1)));
  
  
  BHres[7, ] = c(bhat1, bhat1-1.96*sd1, bhat1+1.96*sd1, pval1);
  # ------------------------------------------------------------------------------
  BRES[[irep]] = BHres;
  
  Ftm <- proc.time();
  S.time <- Ftm[3] - Stm[3];
  print(Sys.time())
  cat( irep,"/",length(PVA)," time: ",S.time,"\n",sep="");
  
  save(BRES, Parares, Bhat0Res,  Bhat1Res, Etasave, file = filename);
}
