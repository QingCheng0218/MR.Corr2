rm(list = ls());
library(MR.LDP);
library(stringr);
library(mr.raps);
library(cause);
library(gsmr);
library(MRMix);
library(MR.Corr2)
NAMEpa = c("ADHD", "MDD", "BIP", "SCZ",  "ASD");
namexy = NAMEpa[1];
QCindex = 1;
pva_cutoff = 1e-4;
QCindex = 1;
lam = 0.9;
ld_r2_thresh = 0.1;
maxIter = 5000;
burnin = 5000;
coreNum = 24;
thin = 10;
# ---------------------------------------------------------------------
filename = paste0(namexy, "Outcomelam", lam, "QC", QCindex, "pva" ,pva_cutoff,".Rdata")
stringname3 = "all_chr_1000G";
block_file = "fourier_ls-all.bed";

fileAllScr = c("PGC.ADHD2012onlyZscore.txt", "MDD_UKB_SummaryStat.txt", "UKB.BIPmatch1KG.txt",
               "SCZ_PGC_EastAsian_SummaryStat.txt", "PGC.ASD2015match1KG.txt");

fileAllXY = c("ADHD_eur_jun2017SummaryStat.txt", "MDD_PGC_SummaryStat.txt", "BIP_PGC_BD1_SummaryStat.txt",
              "SCZ_PGC_EUR_SummaryStat.txt", "ASD_PGC_EUR_SummaryStat.txt");
nameScr = fileAllScr[which(str_detect(fileAllScr, namexy))];
nameX = fileAllXY[which(str_detect(fileAllXY, namexy))];

filescreen = nameScr;
fileexposure = nameX;
fileAllXY0 = fileAllXY[-which(str_detect(fileAllXY, namexy))];


NumSave = maxIter / thin;
BRES = list();
Parares = matrix(0, nrow = 6, ncol = length(fileAllXY0));
Bhat0Res = Bhat1Res = matrix(0, nrow = NumSave, ncol = length(fileAllXY0));
Outcome = rep(0, length(fileAllXY0));

for(yy in 1:length(fileAllXY0)){
  
  Stm <- proc.time();
  
  nameY = fileAllXY0[yy];
  Outcome[yy] = nameY;
  fileoutcome = paste0(pwdXY, nameY);
  
  scrres = matchscreen(filescreen, fileexposure,  fileoutcome, stringname3,  pva_cutoff)
  bh1 = as.numeric(scrres$bh1);
  bh2 = as.numeric(scrres$bh2);
  s12 = as.numeric(scrres$s12);
  s22 = as.numeric(scrres$s22);
  chr = as.numeric(scrres$chr);
  bp = scrres$bp;
  rsname = scrres$rsname
  avbIndex = scrres$idxin;
  idx4panel = scrres$idx4panel;
  p = length(avbIndex);
  
  if(QCindex){
    QCresult = summaryQC(mhcstart, mhcend, bh1, bh2, s12, s22, bp,
                         chr, rsname, avbIndex, idx4panel, Inf, Inf)
    bh1new = QCresult$bh1new;
    bh2new = QCresult$bh2new;
    s12new = QCresult$s12new;
    s22new = QCresult$s22new;
    bpnew = QCresult$bpnew;
    chrnew = QCresult$chrnew;
    avbIndexnew = QCresult$avbIndexnew;
    idx4panelnew = QCresult$idx4panel
    rsnamenew = QCresult$rsnamenew;
  }else{
    bh1new = bh1;
    bh2new = bh2;
    s12new = s12;
    s22new = s22;
    bpnew = bp;
    chrnew = chr;
    rsnamenew = rsname;
    idx4panelnew = idx4panel;
    avbIndexnew = avbIndex;
  }
  p = length(avbIndexnew);

  BHres = matrix(0, nrow = 6, ncol = 4);
  # ------------------------------------------------------------------------------
  #####---- MR-LDCP----#####
  # --------------------------
  print(Sys.time())
  cat("Start MR-LDCP procedure:", "\n");
  # block information
  infres = Cal_blockinf(bpnew, chrnew, block_file);
  nblocks = infres$nblocks;
  a1 = 1;
  b1 = nblocks;
  
  opt = list(agm = 0.001, bgm = 0.001, aal = 0.001, bal = 0.001,
             a = a1, b = b1, maxIter = maxIter, thin = 10, burnin = burnin);
  
  result0 = MRCorr2Real(bpnew, chrnew, avbIndexnew-1,idx4panelnew, block_file,
                        stringname3, ld_r2_thresh, bh1new, bh2new, s12new, s22new, lam, coreNum, opt)
  
  
  Bhat0 = result0$Beta0res;
  bhat1 = mean(Bhat0);
  sd1 = sd(Bhat0);
  tstat1 = bhat1 / sd1;
  pval1 = 2*(1 - pnorm(abs(tstat1)))
  BHres[1, ] = c(bhat1, bhat1-1.96*sd1, bhat1+1.96*sd1, pval1);
  
  Etarate0 = length(which(result0$Eta==1))/result0$nblocks;
  Eta = result0$Etares;
  NB = result0$nblocks;
  Indpid = result0$Indpid + 1;
  p1 = length(Indpid);
  
  sgga2hat = median(result0$Sgga2Res);
  sgal2hat = median(result0$Sgal2Res);
  para = c(p, NB, p1, Etarate0, sgga2hat, sgal2hat);
  Parares[, yy] = para;
  Bhat0Res[, yy] = Bhat0;
  Bhat1Res[, yy] = result0$Beta1res;
  # ------------------------------------------------------------------------------
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
  epsStopLogLik = 1e-7
  lam1 = 0; #for lambda = 0.9 shrinkage.
  RealMRLDP_Hb = MRLDP_RealPXvb_block(bpnew, chrnew, avbIndexnew-1, idx4panelnew, block_file,
                                      stringname3, bh1new, bh2new, s12new, s22new,
                                      gamma, alpha, beta0, sgga2, sgal2, coreNum,
                                      lam1, 0, epsStopLogLik, maxIter1, model = 2);
  
  
  RealMRLDP_H0 = MRLDP_RealPXvb_block(bpnew, chrnew, avbIndexnew-1,idx4panelnew,  block_file,
                                      stringname3, bh1new, bh2new, s12new, s22new,
                                      gamma, alpha, beta0, sgga2, sgal2, coreNum,
                                      lam1, 1, epsStopLogLik, maxIter1, model = 2);
  beta0_MRLDP = RealMRLDP_Hb$beta0;
  Tstat_MRLDP = 2*(RealMRLDP_Hb$tstat - RealMRLDP_H0$tstat);
  MRLDP_se = abs(RealMRLDP_Hb$beta0/sqrt(Tstat_MRLDP));
  pval2 = pchisq(Tstat_MRLDP, 1, lower.tail = F);
  BHres[2, ] = c(beta0_MRLDP, beta0_MRLDP-1.96*MRLDP_se, beta0_MRLDP+1.96*MRLDP_se, pval2);
  
  # ------------------------------------------------------------------------------
  #####---- CAUSE----#####
  # --------------------------
  print(Sys.time())
  cat("Start CAUSE procedure:", "\n");

  scrres4para = matchscreen(fileexposure, fileexposure,  fileoutcome, stringname3,  1)
  X = data.frame(
    snp = scrres4para$rsname,
    beta_hat_1 = scrres4para$bh1,
    seb1 = scrres4para$s12,
    beta_hat_2 = scrres4para$bh2,
    seb2 = scrres4para$s22,
    A1 = rep("A", length(scrres4para$bh1)),
    A2 = rep("G", length(scrres4para$bh1))
  );
  new_cause_data1 <- function(x = data.frame()){
    stopifnot(inherits(x, "data.frame"))
    # x <- validate_cause_data(x)
    #stopifnot(all(c("snp", "beta_hat_1", "seb1", "beta_hat_2", "seb2") %in% names(x)))
    structure(x, class = c("cause_data", "data.frame"))
  };
  X1 = new_cause_data1(X);


  varlist <- with(X1, sample(snp, size=dim(X1)[1], replace=FALSE))
  params <- est_cause_params(X1, varlist);

  top_vars = rsnamenew[Indpid];
  XX = data.frame(
    snp = rsnamenew,
    beta_hat_1 = bh1new,
    seb1 = s12new,
    beta_hat_2 = bh2new,
    seb2 = s22new,
    A1 = rep("A", length(bh1new)),
    A2 = rep("G", length(bh1new))
  )
  XX1 = new_cause_data1(XX);

  resCAUSE <- cause(X=XX1, variants = top_vars, param_ests = params);
  if(class(resCAUSE)!="try-error") {
    qs = summary(resCAUSE)$quants
    BHres[3, ] = c(qs[[2]][1,1], qs[[2]][2,1], qs[[2]][3,1], summary(resCAUSE)$p);
  }

  # ------------------------------------------------------------------------------
  #####---- GSMR----#####
  # --------------------------
  print(Sys.time())
  cat("Start GSMR procedure:", "\n");
  
  Rres = Cal_block_Rmatrix(bpnew, chrnew, avbIndexnew-1,idx4panelnew, block_file,
                           stringname3, ld_r2_thresh, coreNum, lam);
  IndR = Rres$IndR;
  Indpid = Rres$Indpid + 1;
  dat = data.frame(
    snp = rsnamenew[Indpid],
    beta_hat_1 = bh1new[Indpid],
    seb1 = s12new[Indpid],
    beta_hat_2 = bh2new[Indpid],
    seb2 = s22new[Indpid],
    p_value = 2*pnorm(-abs(bh1new[Indpid]/s12new[Indpid])),
    p_value2 = 2*pnorm(-abs(bh2new[Indpid]/s22new[Indpid]))
  )
  
  p_val_thresh  = 1;
  IndR <- data.frame(IndR);
  names(IndR) <- seq_along(Indpid);
  rownames(IndR) <- seq_along(Indpid);
  resgsmr <-  try(with(dat, gsmr(beta_hat_1, seb1, p_value,
                                 beta_hat_2, seb2, bzy_pval = p_value2,
                                 ldrho=IndR, snpid=seq_along(Indpid),
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
  raps <-  try(mr.raps(dat$beta_hat_1, dat$beta_hat_2, dat$seb1, dat$seb2, over.dispersion = TRUE));
  if(class(raps)!="try-error") {
    BHres[6, ] = c(raps$beta.hat, raps$beta.hat - 1.96*raps$beta.se, raps$beta.hat + 1.96*raps$beta.se, raps$beta.p.value)
    
  }
  
  # ------------------------------------------------------------------------------
  BRES[[yy]] = BHres;
  
  Ftm <- proc.time();
  S.time <- Ftm[3] - Stm[3];
  print(Sys.time())
  cat( yy,"/", length(fileAllXY0), " time: ",S.time,"\n",sep="");
  
  save(BRES, Parares, Bhat0Res, Outcome, Bhat1Res, file = filename);
  
}