#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;


// [[Rcpp::export]]
SEXP cf_integral(
    double betaD,
    arma::vec time,
    arma::vec event,
    arma::vec IV,
    arma::vec cv_IV_c,
    arma::mat Covariates,
    arma::mat D_status,
    arma::vec stime,
    arma::mat cv_beta,
    arma::mat cv_dLam,
    int max_iter,
    double tol,
    double eta,
    double contraction
)
{
  int n = time.size();
  int k = stime.size();
  int p = Covariates.n_cols;
  int step;
  bool Convergence;
  double Hessian, fn, new_betaD, pk, fn_abs, diff;
  arma::vec res(n), Covbeta(n), int_cexpbetaD(n), int_cdexpbetaD_1(n), int_cdexpbetaD_2(n), SY(k), dSY(k), new_parameters(p+1), dPhi(n);
  arma::mat int_D(n, k), dNt(n, k), Yt(n, k), int_expbetaD(n, k), int_Lam(n, k), int_dexpbetaD(n, k), int_dexpbetaD_Lam(n, k), int_expbetaD_dLam(n, k), int_expbetaD_dLam_dalpha(n, p), prev_Hessian(p+1, p+1), inv_Hessian(p+1, p+1);
  arma::mat pk_mat(p+1, 1), H_com(1, 1);
  // double H_com_val, bktrkg;
  
  // for (int iter = 0; iter < max_iter; iter++)
  // {
    res.zeros(); Covbeta.zeros(); int_cexpbetaD.zeros(); int_Lam.zeros(); int_cdexpbetaD_1.zeros(); int_cdexpbetaD_2.zeros();
    SY.zeros(); dSY.zeros(); int_D.zeros(); int_expbetaD.zeros(), int_dexpbetaD.zeros(), int_dexpbetaD_Lam.zeros();
    int_expbetaD_dLam_dalpha.zeros();dPhi.zeros(); 
    
    Hessian = 0.0;
    fn = 0.0;
    
    for (int i = 0; i < n; i++)
    {
      for (int ii = 0; ii < p; ii++)
      {
        Covbeta[i] += Covariates(i, ii) * cv_beta(i, ii);
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
        if (Yt(i, j) > 0)
        {
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
            res[i] = res[i] + exp(int_D(i, j) * betaD) * (dNt(i, j)- Yt(i, j) * cv_dLam(i, j)) - int_expbetaD(i, j) * (IV[i] * betaD + Covbeta[i]);
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
            res[i] = res[i] + exp(int_D(i, j) * betaD) * (dNt(i, j) - Yt(i, j) * cv_dLam(i, j)) - Yt(i, j-1) * int_expbetaD(i, j) * (D_status(i, j-1) * betaD + Covbeta[i]);
          }
          
          
          dPhi[i] = dPhi[i] + (dNt(i, j) - cv_dLam(i, j)) * int_D(i, j) * exp(betaD * int_D(i, j)) - 
            Yt(i, j) * int_dexpbetaD(i, j) * (D_status(i, j) * betaD + Covbeta[i]);
        
        }
        
        
        
        
      }
      // printf("%d\n", i);
      
      fn = fn + cv_IV_c[i] * res[i];
      Hessian = Hessian + cv_IV_c[i] * dPhi[i];
    }
    // if (iter > 0)
    // {
    //   H_com = pk*Hessian*pk;
    //   H_com_val = H_com(0, 0);
    //   fn_2 = fn*fn;
    //   prev_fn_2 = prev_fn*prev_fn;
    //   bktrkg = (fn_2/2 - (prev_fn_2)/2)/(pk*fn + H_com_val);
    //   if(bktrkg > eta)
    //   {
    //     betaD = betaD - (1 - contraction) * pk;
    //     if (iter >= max_iter - 1) Convergence = false;
    //     continue;
    //   }
    // }
    // 
    // arma::arma_print(fn);
    // arma::arma_print(Hessian);
    // arma::arma_print(dNt);
    // if(iter == 0)
    //     arma::arma_print(dPhi);
    // pk = -fn/Hessian;
    // 
    // new_betaD = betaD + pk;
    // printf("%f\n", pk);
    // fn_abs = abs(fn);
    // diff = abs(new_betaD - betaD);
    // if(fn_abs < tol || diff < (tol * tol))
    // {
    //   step = iter;
    //   Convergence = true;
    //   break;
    // }
    // betaD = new_betaD;
    // // prev_fn = fn;
    // if (iter >= (max_iter - 1))
    // {
    //   step = iter;
    //   Convergence = false;
    // }
  // }
  
  return Rcpp::List::create(
    Named("fn") = fn
    // Named("Convergence") = Convergence,
    // Named("iter") = step
  );
}
