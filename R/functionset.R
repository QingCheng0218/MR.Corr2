genRawGeno <- function(maf, L, M, rho, n){
  SIGMA = matrix(nrow=M,ncol=M)
  for (i in 1:M){
    for (j in 1:M){
      SIGMA[i,j] = rho^(abs(i-j));
    }
  }
  
  nsnp = L*M;
  X = NULL;
  for ( l in 1:L ){
    
    index = (M*(l-1)+1): (M*l);
    AAprob = maf[index]^2.;
    Aaprob = 2*maf[index]*(1-maf[index]);
    quanti = matrix(c(1-Aaprob-AAprob, 1- AAprob),M,2);
    Xt = rmvnorm(n, mean=rep(0,M), sigma=SIGMA, method="chol")
    Xt2 = matrix(0,n,M);
    for (j in 1:M){
      cutoff = qnorm(quanti[j,]);
      Xt2[Xt[,j] < cutoff[1],j] = 0;
      Xt2[Xt[,j] >= cutoff[1] & Xt[,j] < cutoff[2],j] = 1;  ## attention
      Xt2[Xt[,j] >= cutoff[2],j] = 2;
    }
    X <- cbind(X,Xt2);
  }
  return(X)
}

mhcstart = 28477797;
mhcend = 33448354;
summaryQC = function(mhcstart, mhcend, bh1, bh2, s12, s22, bp, chr,
                     rsname, avbIndex, idx4panel, xbound, ybound){
  # remove SNPs in MHC region:
  idxchr6 = which(chr==6);
  idxcut = idxchr6[which(bp[idxchr6]>=mhcstart & bp[idxchr6]<=mhcend)];
  pmhc = length(idxcut);
  
  if(pmhc!=0){
    bh1Rmhc = bh1[-idxcut];
    bh2Rmhc = bh2[-idxcut];
    s12Rmhc = s12[-idxcut];
    s22Rmhc = s22[-idxcut];
    bpRmhc = bp[-idxcut];
    chrRmhc = chr[-idxcut];
    rsnameRmhc = rsname[-idxcut];
    avbIndexRmhc = avbIndex[-idxcut];
    
    tmp0 = 1:length(bh1);
    tmp = tmp0[-idxcut];
    if(length(idx4panel)!=0){
      idx4panelRmhc = match(avbIndex[intersect((idx4panel + 1), tmp)], avbIndexRmhc) -1
    }else{
      idx4panelRmhc = idx4panel;
    }
    
  }else{
    bh1Rmhc = bh1;
    bh2Rmhc = bh2;
    s12Rmhc = s12;
    s22Rmhc = s22;
    bpRmhc = bp;
    chrRmhc = chr;
    rsnameRmhc = rsname;
    avbIndexRmhc = avbIndex;
    idx4panelRmhc = idx4panel;
  }
  
  
  # remove SNPs(exposure) with chi-square >80
  idx = which((bh1Rmhc/s12Rmhc)^2>xbound);
  px = length(idx);
  if(px!=0){
    bh1Rmhc_x = bh1Rmhc[-idx];
    bh2Rmhc_x = bh2Rmhc[-idx];
    s12Rmhc_x = s12Rmhc[-idx];
    s22Rmhc_x = s22Rmhc[-idx];
    bpRmhc_x = bpRmhc[-idx];
    chrRmhc_x = chrRmhc[-idx];
    rsnameRmhc_x = rsnameRmhc[-idx];
    avbIndexRmhc_x = avbIndexRmhc[-idx];
    # idx4panelRmhc_x = idx4panelRmhc[-idx];
    
    
    tmp0 = 1:length(bh1Rmhc);
    tmp = tmp0[-idx];
    if(length(idx4panel)!=0){
      idx4panelRmhc_x = match(avbIndexRmhc[intersect((idx4panelRmhc + 1), tmp)], avbIndexRmhc_x) -1;
    }else{
      idx4panelRmhc_x = idx4panel;
    }
    
    
  }else{
    bh1Rmhc_x = bh1Rmhc;
    bh2Rmhc_x = bh2Rmhc;
    s12Rmhc_x = s12Rmhc;
    s22Rmhc_x = s22Rmhc;
    bpRmhc_x = bpRmhc;
    chrRmhc_x = chrRmhc;
    rsnameRmhc_x = rsnameRmhc;
    avbIndexRmhc_x = avbIndexRmhc;
    idx4panelRmhc_x = idx4panelRmhc;
  }
  
  # remove SNPs(outcome) with chi-square >80
  idy = which((bh2Rmhc_x/s22Rmhc_x)^2>ybound);
  py = length(idy);
  if(py!=0){
    bh1Rmhc_xy = bh1Rmhc_x[-idy];
    bh2Rmhc_xy = bh2Rmhc_x[-idy];
    s12Rmhc_xy = s12Rmhc_x[-idy];
    s22Rmhc_xy = s22Rmhc_x[-idy];
    bpRmhc_xy = bpRmhc_x[-idy];
    chrRmhc_xy = chrRmhc_x[-idy];
    rsnameRmhc_xy = rsnameRmhc_x[-idy];
    avbIndexRmhc_xy = avbIndexRmhc_x[-idy];
    # idx4panelRmhc_xy = idx4panelRmhc_x[-idy];
    
    tmp0 = 1:length(bh1Rmhc_x);
    tmp = tmp0[-idx];
    
    if(length(idx4panel)!=0){
      idx4panelRmhc_xy = match(avbIndexRmhc_x[intersect((idx4panelRmhc_x + 1), tmp)], avbIndexRmhc_xy) -1;
    }else{
      idx4panelRmhc_xy = idx4panel;
    }
    
  }else{
    bh1Rmhc_xy = bh1Rmhc_x;
    bh2Rmhc_xy = bh2Rmhc_x;
    s12Rmhc_xy = s12Rmhc_x;
    s22Rmhc_xy = s22Rmhc_x;
    bpRmhc_xy = bpRmhc_x;
    chrRmhc_xy = chrRmhc_x;
    rsnameRmhc_xy = rsnameRmhc_x;
    avbIndexRmhc_xy = avbIndexRmhc_x;
    idx4panelRmhc_xy = idx4panelRmhc_x;
  }
  return(list(bh1new = bh1Rmhc_xy, bh2new = bh2Rmhc_xy, s12new = s12Rmhc_xy, s22new = s22Rmhc_xy,
              bpnew = bpRmhc_xy, chrnew = chrRmhc_xy, rsnamenew = rsnameRmhc_xy,
              avbIndexnew = avbIndexRmhc_xy, idx4panelnew = idx4panelRmhc_xy, pmhc = pmhc, px = px, py = py))
}

traceplot <- function(bhatpoint){
  y <- bhatpoint
  x <- 1:length(bhatpoint);
  da <- cbind(x, y);
  dat <- data.frame(da);
  p1 <- ggplot(data = dat, aes(x= x, y = y))+  geom_line()  +
    labs( x = paste0("GibbsSampleIndex"), y =  expression(hat(beta[0])));
  p1 = p1 + theme(axis.title.x = element_text(size=10,face = "bold"),
                  axis.text.x = element_text(size=12,face = "bold"),
                  axis.title.y = element_text(size=10,face = "bold"),
                  axis.text.y = element_text(size=12,face = "bold"));
  return(p1);
}
