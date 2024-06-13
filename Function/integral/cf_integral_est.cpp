#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;


// [[Rcpp::export]]
SEXP cf_integral_est(
    double betaD,
    arma::vec time,
    arma::vec event,
    arma::vec IV,
    arma::vec cf_IV_c,
    arma::mat Covariates,
    arma::mat D_status,
    arma::vec stime,
    arma::mat cf_beta,
    arma::mat cf_dLam,
    int max_iter,
    double tol,
    double eta,
    double contraction
)
{
  int n = time.size();
  int k = stime.size();
  int p = Covariates.n_cols;
  int step = 0;
  bool Convergence;
  double Hessian, fn, new_betaD, pk, fn_abs, diff, Asyvar, hatsigma, hatphi, ave_hatsigma;
  arma::vec hatsigma_v(n), res(n), Covbeta(n), SY(k), dSY(k), dPhi(n);
  arma::mat int_D(n, k), dNt(n, k), Yt(n, k), int_expbetaD(n, k), int_dexpbetaD(n, k);
  // double H_com_val, bktrkg;
  
  for (int iter = 0; iter < max_iter; iter++)
  {
    fn = 0.0;
    Hessian = 0.0;
    res.zeros(); Covbeta.zeros(), dPhi.zeros();
    SY.zeros(); dSY.zeros(); int_D.zeros(); int_expbetaD.zeros(), int_dexpbetaD.zeros();
    
    // arma::arma_print(beta);
    // if(iter == 0) 
    // {
    //     arma::arma_print(betaD);
    //     arma::arma_print(beta);
    // }
    for (int i = 0; i < n; i++)
    {
      for (int ii = 0; ii < p; ii++)
      {
        Covbeta[i] += Covariates(i, ii) * cf_beta[ii];
      }
    }
    
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < k; j++)
      {
        if(j == 0)
          int_D(i, j) = IV[i] * stime[j];
        else 
          int_D(i, j) = int_D(i, j-1) + D_status(i, j-1) * (stime[j] - stime[j-1]);
        
        
        dNt(i, j) = (time[i] == stime[j]) ? 1:0;
        dNt(i, j) *= event[i];
        Yt(i, j) = (time[i] >= stime[j]) ? 1:0;
        SY[j] = SY[j] + exp(betaD * int_D(i, j)) * Yt(i, j);
        dSY[j] = dSY[j] + int_D(i, j) * exp(betaD * int_D(i, j)) * Yt(i, j);
      }
      
    }
    
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < k; j++)
      {
        // Integral computation for exp(D\theta) and D * exp(D\theta)
        if (j == 0)
        {
          if (betaD != 0)
          {
            int_expbetaD(i, j) = (IV[i] > 0.5 ) ? (exp(betaD * int_D(i, j))-1)/betaD:stime[j];
            int_dexpbetaD(i, j) = (IV[i] > 0.5) ? (int_D(i, j) * exp(betaD * int_D(i, j))/betaD - int_expbetaD(i, j)/betaD):0;
          }
          else
          { 
            int_expbetaD(i, j) = stime[j];
            int_dexpbetaD(i, j) = int_D(i, j) * int_D(i, j) / 2;
          }
          res[i] = res[i] + exp(int_D(i, j) * betaD) * (dNt(i, j)- Yt(i, j) * cf_dLam(i,j)) - int_expbetaD(i, j) * (IV[i] * betaD + Covbeta[i]);

              
          dPhi[i] = dPhi[i]+ (dNt(i, j) - Yt(i, j) * cf_dLam(i,j)) * int_D(i, j) * exp(betaD * int_D(i, j)) -
            int_dexpbetaD(i, j) * (IV[j] * betaD + Covbeta[i]) -
            int_expbetaD(i, j) * IV[j];
        }
        else
        {
          if (betaD != 0)
          {
            int_expbetaD(i, j) = (D_status(i, j-1) > 0.5) ? (exp(betaD * int_D(i, j)) - exp(betaD * int_D(i, j-1)))/betaD:(stime[j] - stime[j-1]) * exp(betaD * int_D(i, j));
            int_dexpbetaD(i, j) = (D_status(i, j-1) > 0.5) ? ((int_D(i, j) * exp(betaD * int_D(i, j)) - int_D(i, j-1) * exp(betaD * int_D(i, j-1)))/betaD - int_expbetaD(i, j)/betaD) : (int_D(i, j)*exp(betaD * int_D(i, j))*(stime[j] - stime[j-1]));
          }
          else 
          {
            int_expbetaD(i, j) = stime[j] - stime[j-1];
            int_dexpbetaD(i, j) = (stime[j] - stime[j-1]) * int_D(i, j-1) + D_status(i, j-1) * (stime[j] - stime[j-1]) * (stime[j] - stime[j-1])/2;
            
          }
          res[i] = res[i] + exp(int_D(i, j) * betaD) * (dNt(i, j) - Yt(i, j)*cf_dLam(i, j)) - Yt(i, j-1) * int_expbetaD(i, j) * (D_status(i, j-1) * betaD + Covbeta[i]);
              
          dPhi[i] = dPhi[i] + (dNt(i, j) - Yt(i, j) * cf_dLam(i, j)) * int_D(i, j) * exp(betaD * int_D(i, j)) -
            Yt(i, j-1) * int_dexpbetaD(i, j) * (D_status(i, j-1) * betaD + Covbeta[i]) -
            Yt(i, j-1) * int_expbetaD(i, j) * D_status(i, j-1);
        }
        
        
        
        
        
        // if (j == 1 && i ==0) arma::arma_print(res);
      }
      fn = fn + cf_IV_c[i] * res[i];
      Hessian = Hessian + cf_IV_c[i] * dPhi[i];
    }
    
    pk = -fn/Hessian;
    new_betaD = betaD + pk;
    fn_abs = abs(fn);
    // arma::arma_print(fn_abs);
    // arma::arma_print(arma::sum(fn_abs));
    diff = abs(new_betaD - betaD);
    if(fn_abs < tol ||diff < (tol * tol))
    {
      step = iter;
      Convergence = true;
      break;
    }
    betaD = new_betaD;
    // prev_fn = fn;
    if (iter >= (max_iter - 1))
    {
      step = iter;
      Convergence = false;
    }
  }
  
  Asyvar = 0.0;
  hatsigma = 0.0;
  hatphi = 0.0;
  ave_hatsigma = 0.0;
  hatsigma_v.zeros();
  
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < k; j++)
    {
      if (j == 0)
      {
        ave_hatsigma = ave_hatsigma + (cf_IV_c[i]*exp(int_D(i, j) * betaD) *(dNt(i, j) - Yt(i, j) * cf_dLam(i, j)) - 
          cf_IV_c[i] * int_expbetaD(i, j)* 1*(betaD* IV[i] + Covbeta[i]))/n;
      } else
      {
        ave_hatsigma = ave_hatsigma + (cf_IV_c[i]*exp(int_D(i, j) * betaD) *(dNt(i, j)- Yt(i, j) * cf_dLam(i, j)) - 
          cf_IV_c[i] * int_expbetaD(i, j)* Yt(i, j-1)*(betaD*D_status(i, j-1) + Covbeta[i]))/n;
      }
    }
  }
  
  
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < k; j++)
    {
      // hatsigma = hatsigma + ((IV_c[i] * IV_c[i])*exp(2*int_D(i, j)*betaD)*dNt(i, j));
      
      if(j==0)
      {
        hatsigma_v[i] = hatsigma_v[i] + (cf_IV_c[i]*exp(int_D(i, j) * betaD) *(dNt(i, j)- Yt(i, j) * cf_dLam(i,j)) - 
          cf_IV_c[i] * int_expbetaD(i, j)* 1*(betaD*IV[i] + Covbeta[i]));
        hatphi = hatphi + (cf_IV_c[i] * int_D(i, j)* exp(int_D(i,j)*betaD)*(dNt(i, j) - Yt(i, j)*cf_dLam(i,j)) -
          cf_IV_c[i] * int_dexpbetaD(i, j)*1*(betaD*IV[i] + Covbeta[i]) -
          cf_IV_c[i] * int_expbetaD(i, j)*1*IV[i]);
      } else 
      {
        hatsigma_v[i] = hatsigma_v[i] + cf_IV_c[i]*exp(int_D(i, j) * betaD) *(dNt(i, j)- Yt(i, j) * cf_dLam(i,j)) - 
          cf_IV_c[i] * int_expbetaD(i, j)* Yt(i, j-1)*(betaD*D_status(i, j-1) + Covbeta[i]);
        hatphi = hatphi + (cf_IV_c[i] * int_D(i, j)* exp(int_D(i,j)*betaD)*(dNt(i, j)- Yt(i, j)*cf_dLam(i,j)) -
          cf_IV_c[i] * int_dexpbetaD(i, j)*Yt(i, j-1)*(betaD*D_status(i, j-1) + Covbeta[i]) -
          cf_IV_c[i] * int_expbetaD(i, j)*Yt(i, j-1)*D_status(i, j-1));
      }
      
    }
    hatsigma = hatsigma + (hatsigma_v[i] - ave_hatsigma)*(hatsigma_v[i] - ave_hatsigma);
  }
  Asyvar = hatsigma / (hatphi*hatphi);
  
  
  
  return Rcpp::List::create(
    Named("x") = betaD,
    Named("Convergence") = Convergence,
    Named("iter") = step,
    Named("var") = Asyvar
  );
}
