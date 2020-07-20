//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "pdsoft.hpp"

using namespace std;
using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH, bigmemory)]]
// [[Rcpp::plugins(cpp11)]]


double **make_mat(int nrow, int ncol);
void delete_mat( double **mat);
double *make_vec(int len);
void delete_vec(double *vec);         

void pdsc(double * S, double *Sigma, double * Omega, double * tosoft, int * pin, 
          double * lam, 
          double * tauin,
          double * tolin, 
          int * maxitin, 
          double * tolout, 
          int * maxitout, 
          int * totalout)
{
  int p=*pin;
  double tau=*tauin;  
  int i,j,jj,k, numit_in, numit_out;
  double diff_in, diff_out, tmp, newbeta, betadiff, tau_over_Sigma_jj, newtmp, val;
  
  double * S_jj, *Omega_jj, *Sigma_jj, * S_0j, *lam_0j, *Sigma_0j, *Omega_0j,
  *Omega_j0, *Sigma_j0,
  * beta_k, *beta_jj, *tosoft_jj, *tosoft_k, 
  * Omega_k_jj, *Omega_jj_jj,  *S_jj_j,
  * lam_jj_j, *Omega_0_jj,
  * Omega_jj_j, *Omega_j_jj, *Sigma_j_jj,       
  * Omega_kj;
  
  
  
  numit_out=0;
  diff_out=*tolout+1;
  while( (diff_out > *tolout) && (numit_out < *maxitout) )
  {
    numit_out++;
    diff_out = 0;
    
    S_jj=S;
    S_0j=S;
    lam_0j=lam; 
    Omega_jj=Omega;
    Sigma_jj=Sigma;
    Sigma_0j=Sigma;
    Omega_0j=Omega;
    Omega_j0=Omega;
    Sigma_j0=Sigma;
    
    for(j=0; j < p; j++)
    {
      // update jth diagonal entry of Sigma
      tmp = *S_jj + tau* *Omega_jj; 
      diff_out += fabs(tmp - *Sigma_jj);
      *Sigma_jj = tmp;
      tau_over_Sigma_jj=tau / *Sigma_jj; 
      
      numit_in=0;
      diff_in=*tolin+1;      
      while( (diff_in > *tolin) && (numit_in < *maxitin) )
      {
        numit_in++;
        diff_in = 0;
        
        //first run:
        
        if(numit_in==1)
        {
          tosoft_jj=tosoft;    
          Omega_0_jj=Omega;
          Omega_jj_jj=Omega;
          S_jj_j=S_0j;
          beta_jj=Sigma_0j;
          lam_jj_j=lam_0j;
          for(jj=0; jj < p; jj++)
          {        
            if(jj != j)
            {
              tmp = 0;
              beta_k=Sigma_0j;
              Omega_k_jj=Omega_0_jj;
              for(k=0; k<p; k++)
              {
                if(k != jj && k != j)
                {
                  tmp+= *Omega_k_jj * *beta_k;  
                }
                beta_k++;
                Omega_k_jj++;
              }
              
              tmp = *S_jj_j - tau_over_Sigma_jj *tmp;
              *tosoft_jj=tmp;        
              newtmp=0;
              val = fabs(tmp) - *lam_jj_j;
              if( val > 0)
                if(tmp > 0)
                  newtmp = val;
                else
                  newtmp = -val;
                tmp=newtmp;  
                
                
                newbeta = tmp/( 1 + tau_over_Sigma_jj * *Omega_jj_jj);
                diff_in += fabs(newbeta - *beta_jj);
                *beta_jj = newbeta;
            } 
            lam_jj_j++;
            beta_jj++;
            S_jj_j++;
            Omega_jj_jj+=p+1;
            Omega_0_jj+=p;
            tosoft_jj++;
          }
        } else
        {
          /*  Update tosoft only if beta_k changes */    
          tosoft_jj=tosoft;
          Omega_jj_jj=Omega;
          S_jj_j=S_0j;
          beta_jj=Sigma_0j;
          lam_jj_j=lam_0j;     
          Omega_0_jj=Omega;
          for(jj=0; jj < p; jj++)
          {        
            if(jj != j)
            {              
              newtmp=0;
              val = fabs(*tosoft_jj) - *lam_jj_j;
              if( val > 0)
                if(*tosoft_jj > 0)
                  newtmp = val;
                else
                  newtmp = -val;
                tmp=newtmp;  
                
                newbeta = tmp/( 1 + tau_over_Sigma_jj * *Omega_jj_jj);
                
                
                if(newbeta != *beta_jj ) // then update tosoft
                {
                  
                  betadiff=*beta_jj - newbeta;
                  tosoft_k=tosoft;
                  Omega_k_jj=Omega_0_jj;
                  for(k=0; k<p; k++)
                  {
                    if(k != jj && k != j)
                    {
                      *tosoft_k +=  tau_over_Sigma_jj* *Omega_k_jj * betadiff;  
                    }
                    tosoft_k++;
                    Omega_k_jj++;
                  }
                  diff_in+=fabs(betadiff);
                  *beta_jj=newbeta;
                } 
            }
            beta_jj++;
            Omega_0_jj+=p;            
            tosoft_jj++;
            lam_jj_j++;
            Omega_jj_jj+=p+1;
          } // end the jj for loop
        } /* End: update tosoft only if beta_k changes */
          diff_out+=diff_in;
      } // end inner while
      
      /* the lasso regression for the jth column is complete */   
      
      // update the jth row/column of Omega (no diagonal entry)
      // update the jth row of Sigma  for symmetry
      Omega_0_jj=Omega;
      Omega_jj_j=Omega_0j;
      Omega_j_jj=Omega_j0;
      Sigma_j_jj=Sigma_j0;
      beta_jj=Sigma_0j;
      for(jj=0; jj <p; jj++)
      {
        if(jj != j)
        {
          tmp=0;
          beta_k=Sigma_0j;
          Omega_k_jj=Omega_0_jj;
          for(k=0; k < p; k++)
          {
            if(k != j)
              tmp+= *Omega_k_jj * *beta_k;
            beta_k++;
            Omega_k_jj++;
          }
          tmp=-tmp/ *Sigma_jj;
          *Omega_jj_j= tmp;
          *Omega_j_jj = tmp;
          *Sigma_j_jj = *beta_jj;
        } 
        Omega_0_jj+=p;        
        beta_jj++;        
        Omega_jj_j++;
        Omega_j_jj+=p;
        Sigma_j_jj+=p;        
      }
      
      // update Omega jj
      tmp=0;
      beta_k=Sigma_0j;
      Omega_kj=Omega_0j;
      for(k=0; k < p; k++)
      {
        if(k != j)
          tmp = tmp + *beta_k * *Omega_kj;
        Omega_kj++;
        beta_k++;
      }
      tmp = 1 - tmp;
      tmp = tmp/ *Sigma_jj;
      *Omega_jj = tmp;
      
      
      S_jj+=p+1;
      S_0j+=p;
      Omega_jj+=p+1;
      Sigma_jj+=p+1;
      Sigma_0j+=p;
      lam_0j+=p;
      Omega_0j+=p;
      Omega_j0++;
      Sigma_j0++;
      
    }// end for loop over j   
  } // end outer while
  totalout[0] = numit_out;
}



void bchol(double * xin, int * nin, int * pin, int * kin, double * bcov)
{
  int p=*pin;
  int n=*nin;
  int k=*kin;
  int i,ii,iii,j,kk,q;
  double xtx,xty,rss,tmpsum,sum,tmp,tmp2,d,a,b;
  double **L=make_mat(p,p); 
  /* L is the covariance cholesky factor, 
  it's diagonal is D */
  double **x=make_mat(n,p);
  double **r=make_mat(n,p); //matrix of residuals
  
  // Read in input 
  for(j=0; j < p; j++)
  {
    for ( i=0; i < n; i++)
    {
      x[i][j]=xin[j*n+i];
      r[i][j]=xin[j*n+i];
    }
  }
  
  // get first residual variance
  rss=0;
  for(ii=0;ii < n; ii++)
  {
    rss+=x[ii][0]*x[ii][0];
  }  
  L[0][0] = rss/n;
  
  
  // go thru rows of cholesky factor (starting with the second row)  
  for( i=1; i <p; i++)
  {
    //get coefficients in the ith row of L
    for(j=1; (j<=k) && ((i-j) >= 0); j++)
    {
      xtx=0;
      xty=0;
      for(ii=0; ii < n; ii++)
      {
        xtx += r[ii][i-j] *r[ii][i-j];
        xty += r[ii][i-j] *x[ii][i];
      }
      L[i][i-j] = xty/xtx; 
    }
    //get residuals
    rss=0;
    for(ii=0; ii < n; ii++)
    {
      sum=0;
      for(iii=1; (iii <= k) && ((i-iii) >= 0); iii++)
      {
        sum += r[ii][i-iii] * L[i][i-iii];
      }
      r[ii][i] = x[ii][i]-sum;
      rss += r[ii][i]*r[ii][i];	   
    }
    L[i][i]=rss/n;
  }
  
  for(i=1; i <=p; i++)
  {
    for(j = i; (i-j) <= k && j >= 1; j--)
    {
      tmpsum=0;
      for(kk=0; (kk <= k) && (kk <= (k-i+j)) &&((j-kk) >= 1); kk++)
      {
        //for cov i,j
        q = j-kk;
        d = L[q-1][q-1];
        a = L[i-1][q-1];
        b = L[j-1][q-1];
        if ( i == q)
        {
          a=1;
        }
        if( j == q)
        {
          b=1;
        }
        tmpsum += a * d * b;
      }
      bcov[(j-1)*p+(i-1)] = tmpsum;
      bcov[(i-1)*p+(j-1)] = tmpsum;
    }
  }
  delete_mat(L);
  delete_mat(x);
  delete_mat(r);
}


double **make_mat(int nrow, int ncol)
{
  double ** mat;
  int k;
  
  mat = (double **) malloc(nrow*sizeof(double*));
  mat[0]=(double*) malloc(nrow*ncol*sizeof(double));
  
  for(k=1; k < nrow; k++)
    mat[k] = mat[k-1] + ncol;
  return mat;
}
void delete_mat( double **mat)
{
  free(mat[0]);
  free(mat);  
}


double *make_vec(int len)
{
  double * vec;
  vec = (double *) malloc(len*sizeof(double));
  return vec;
}

void delete_vec(double *vec)
{
  free(vec);  
}



pdsoftObj  pdsoft(mat s, mat lam,  double tau, string init, bool standard, double tolin, double tolout, int maxitin, int maxitout, bool quiet){
  
  int p= s.n_rows;
  
  vec oe_val = vec();
  mat oe_vec = mat();
  mat S = mat();
  mat sigma = mat();
  mat omega = mat();
  vec dhat = vec();
  vec dhat_inv = vec();
  mat theta = mat();
  mat theta_inv = mat();
  mat s0 = mat();
  mat i0 = mat();
  
  // output object
  pdsoftObj out;
    
  if(sum(sum(lam))==0) init="dense";
  
  if(standard)
  {
    // correlation matrix S
    dhat = sqrt(diagvec(s));
    dhat_inv = 1/dhat;
    S= repmat(dhat_inv, 1, p) % s % repmat(dhat_inv.t(), p, 1);
    
    if(init=="diag")
    {
      s0=eye(p, p);
      i0=eye(p, p);
    } else if(init=="soft")
    {
      mat tmp=abs(S) - lam;
      tmp=tmp%(tmp > 0);
      mat S_soft=sign(S)%tmp;
      S_soft.diag()=diagvec(S);
      eig_sym(oe_val, oe_vec, S_soft);
      vec evs=oe_val/2 + sqrt(oe_val%oe_val + 4*tau) /2;
      s0 = (oe_vec%repmat(evs.t(), p, 1))*oe_vec.t();
      //s0=tcrossprod(oe$vec*rep(evs, each=p), oe$vec)
      i0 = (oe_vec%repmat((1/evs).t(), p, 1))*oe_vec.t();
      //i0=tcrossprod(oe$vec*rep(1/evs, each=p), oe$vec)
    } else if( init=="dense")
    {
      eig_sym(oe_val, oe_vec, S);
      vec evs=oe_val/2 + sqrt( oe_val%oe_val + 4*tau) /2;
      s0 = (oe_vec%repmat(evs.t(), p, 1))*oe_vec.t();
      i0 = (oe_vec%repmat((1/evs).t(), p, 1))*oe_vec.t();
    } else
    {
      if( s0.is_empty() & (!i0.is_empty() ))
        s0=inv(i0);
        if( (!s0.is_empty()) & i0.is_empty() )
          i0=inv(s0);
          if(s0.is_empty() & i0.is_empty() )
            throw std::range_error("Error: must specify either s0 or i0 with user initialization");
    }
  }

  // if(is.null(dim(lam)[1]))
  // {
  //   lam = matrix(lam, nrow=p, ncol=p)
  // }
  
  if( sum(sum(lam)) > 0 )
  {
    mat Soff=S;
    Soff.diag().fill(0);
    double tolmult=sum(sum(abs(Soff)))/2;
    S.reshape(p*p, 1);
    i0.reshape(p*p, 1);
    s0.reshape(p*p, 1);
    lam.reshape(p*p, 1);
    double b = tau;
    tolin = tolin*(tolmult/p);
    tolout = tolout*tolmult;
    int totalout = 1;
    vec tosoft=zeros(p);
    
    pdsc(S.begin(), s0.begin(), i0.begin(), tosoft.begin(), &p, lam.begin(), &b, &tolin, &maxitin, &tolout, &maxitout, &totalout);
    
    i0.reshape(p, p);
    s0.reshape(p, p);
    
    sigma = s0;
    omega = i0;
    
    if(!quiet)
    {
      cout << "Total outer iterations = " << totalout << endl;
    }
  } else // lam == 0
  {
    sigma=s0;
    omega=i0;
  }
  
  
  if(standard)
  {
    theta=sigma;
    theta_inv=omega;
    sigma = repmat(dhat, 1, p) % theta % repmat(dhat.t(), p, 1);
    omega = repmat(dhat_inv, 1, p) % theta_inv % repmat(dhat_inv.t(), p, 1);
  } else
  {
    theta.reset();
    theta_inv.reset();
  }
  
  out.omega = omega;
  out.sigma = sigma;
  out.theta = theta;
  out.theta_inv = theta_inv;

  return out;
}


 
