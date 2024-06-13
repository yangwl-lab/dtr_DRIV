#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;


// [[Rcpp::export]]
SEXP integral2_dLam(
    arma::vec init_parameters,
    arma::vec time,
    arma::vec event,
    arma::vec IV,
    arma::vec IV_c,
    arma::mat Covariates,
    arma::mat D_status,
    arma::vec stime,
    int max_iter,
    double tol,
    double eta,
    double contraction
)
{
  int n = time.size();
  int k = stime.size();
  int p = init_parameters.size() - 1;
  // printf("%d %d %d", n, k, p);
  
  double betaD;
  arma::vec beta(p), res(n), Covbeta(n), fn(p+1), SY(k), dSY(k),A(k), B(k), dLam(k);
  arma::mat int_D(n, k), dNt(n, k), Yt(n, k), K(n, k), int_expbetaD(n, k), int_dexpbetaD(n, k);
  // double H_com_val, bktrkg;
  
  
  res.zeros(); Covbeta.zeros();
  fn.zeros(); SY.zeros(); dSY.zeros(); int_D.zeros(); int_expbetaD.zeros(), int_dexpbetaD.zeros();
  A.zeros(), B.zeros();
  betaD = init_parameters[0];
  beta = init_parameters.subvec(1, p); // index is only adapted to the one dimensional covariates
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
      Covbeta[i] += Covariates(i, ii) * beta[ii];
    }
  }
  
  // arma::arma_print(Covbeta);
  
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
      K(i, j) = exp(betaD * int_D(i, j));
      if (j == 0)
      {
        A[j] = A[j] + (1 - IV[i]);
        B[j] = B[j] + (IV[i]) ;
      } else {
        A[j] = A[j] + K(i, j-1) * (1 - D_status(i, j-1)) * Yt(i, j-1);
        B[j] = B[j] + K(i, j-1) * D_status(i, j-1) * Yt(i, j-1);
      }
    }
  }
  // arma::arma_print(SY);
  // printf("dd\n");
  
  dLam.zeros();
  
  // for (int i = 0; i < n; i++)
  // {
  //   for (int j = 0; j < k; j++)
  //   {
  //     if (j == 0)
  //     {
  //       if (betaD == 0)
  //       {
  //         dLam[j] = dLam[j] + Yt(i, j) * dNt(i, j) / SY[j] - 1 * Covbeta[i] * stime[j] / n;
  //       } else{
  //         dLam[j] = dLam[j] + Yt(i, j) * exp(int_D(i, j) * betaD) * Yt(i, j) * dNt(i, j) / SY[j] -
  //           ((IV[i]*betaD + Covbeta[i]) * IV[i]) * log((exp(betaD*stime[j])*B[j] + A[j])/(B[j]+A[j])) / (betaD*B[j]) -
  //           ((IV[i]*betaD + Covbeta[i]) * (1-IV[i])) * log((A[j] + B[j])/(B[j] + A[j] * exp(-betaD*stime[j]))) / (betaD * A[j]);
  //       }
  //     } else{
  //       if (betaD == 0)
  //       {
  //         dLam[j] = dLam[j]+ Yt(i, j) * dNt(i, j) / SY[j] - Yt(i, j-1) * Covbeta[i] * (stime[j] - stime[j-1]) / SY[j-1];
  //       } else{
  //         // if(j == k-1) printf("%f\n", B[j]);
  //         if(B[j] != 0 && A[j] != 0){
  //           dLam[j] = dLam[j] + Yt(i, j) * exp(int_D(i, j) * betaD) * Yt(i, j) * dNt(i, j) / SY[j] -
  //             (K(i, j-1) * Yt(i, j-1) * (D_status(i, j-1)*betaD + Covbeta[i]) * D_status(i, j-1)) * log((exp(betaD*(stime[j] - stime[j-1]))*B[j] + A[j])/(B[j]+A[j])) / (betaD*B[j]) -
  //             (K(i, j-1) * Yt(i, j-1) * (D_status(i, j-1)*betaD + Covbeta[i]) * (1-D_status(i, j-1))) * log((A[j] + B[j])/(B[j] + A[j] * exp(betaD*(stime[j-1] - stime[j])))) / (betaD * A[j]);
  //         } else{
  //           if(B[j] == 0 && A[j] != 0)
  //           {
  //             dLam[j] = dLam[j] + Yt(i, j) * exp(int_D(i, j) * betaD) * Yt(i, j) * dNt(i, j) / SY[j] -
  //               (K(i, j-1) * Yt(i, j-1) * (D_status(i, j-1)*betaD + Covbeta[i]) * (1-D_status(i, j-1))) * (stime[j] - stime[j-1])/A[j];
  //           } else {
  //             dLam[j] = dLam[j] + Yt(i, j) * exp(int_D(i, j) * betaD) * Yt(i, j) * dNt(i, j) / SY[j] -
  //               (K(i, j-1) * Yt(i, j-1) * (D_status(i, j-1)*betaD + Covbeta[i]) * D_status(i, j-1)) * (stime[j] - stime[j-1]) / (B[j]);
  //           }
  // 
  //         }
  //       }
  //     }
  //   }
  // }
  
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < k; j++)
    {
      if (j == 0)
      {
        dLam[j] = dLam[j] + Yt(i, j) * exp(int_D(i, j) * betaD) * (dNt(i, j) - D_status(i, j)*betaD*stime[j] - Covbeta[i]*stime[j])/ SY[j];
      } else
      {
        dLam[j] = dLam[j] + Yt(i, j) * exp(int_D(i, j) * betaD) * (dNt(i, j) - D_status(i, j)*betaD*(stime[j] - stime[j-1]) - 
          Covbeta[i]*(stime[j] - stime[j-1]))/ SY[j];
      }
    }
  }
  
  
  
  // for (int i = 0; i < n; i++)
  // {
  //   for (int j = 0; j < k; j++)
  //   {
  //     // Integral computation for exp(D\theta) and D * exp(D\theta)
  //     if (j == 0)
  //     {
  //       if (betaD != 0)
  //       {
  //         int_expbetaD(i, j) = (IV[i] > 0.5 ) ? (exp(betaD * int_D(i, j))-1)/betaD:stime[j];
  //         int_dexpbetaD(i, j) = (IV[i] > 0.5) ? (int_D(i, j) * exp(betaD * int_D(i, j))/betaD - int_expbetaD(i, j)/betaD):0;
  //       }
  //       else
  //       { 
  //         int_expbetaD(i, j) = stime[j];
  //         int_dexpbetaD(i, j) = int_D(i, j) * int_D(i, j) / 2;
  //       }
  //       res[i] = res[i] + exp(int_D(i, j) * betaD) * (dNt(i, j)- Yt(i, j) * dLam[j]) - int_expbetaD(i, j) * (IV[i] * betaD + Covbeta[i] + 0.25);
  //     }
  //     else
  //     {
  //       if (betaD != 0)
  //       {
  //         int_expbetaD(i, j) = (D_status(i, j-1) > 0.5) ? (exp(betaD * int_D(i, j)) - exp(betaD * int_D(i, j-1)))/betaD:(stime[j] - stime[j-1]) * exp(betaD * int_D(i, j));
  //         int_dexpbetaD(i, j) = (D_status(i, j-1) > 0.5) ? ((int_D(i, j) * exp(betaD * int_D(i, j)) - int_D(i, j-1) * exp(betaD * int_D(i, j-1)))/betaD - int_expbetaD(i, j)/betaD) : (int_D(i, j)*exp(betaD * int_D(i, j))*(stime[j] - stime[j-1]));
  //       }
  //       else 
  //       {
  //         int_expbetaD(i, j) = stime[j] - stime[j-1];
  //         int_dexpbetaD(i, j) = (stime[j] - stime[j-1]) * int_D(i, j-1) + D_status(i, j-1) * (stime[j] - stime[j-1]) * (stime[j] - stime[j-1])/2;
  //         
  //       }
  //       res[i] = res[i] + exp(int_D(i, j) * betaD) * (dNt(i, j) - Yt(i, j) * dLam[j]) - Yt(i, j-1) * int_expbetaD(i, j) * (D_status(i, j-1) * betaD + Covbeta[i]);
  //     }
  //     
  //     
  //     
  //     
  //     
  //     // if (j == 1 && i ==0) arma::arma_print(res);
  //   }
  //   // printf("%d\n", i);
  //   for (int kk = 0; kk < p+1; kk++)
  //   {
  //     if (kk == 0)
  //     {
  //       fn[kk] = fn[kk] + IV_c[i] * res[i];
  //     }
  //     else
  //     {
  //       fn[kk] = fn[kk] + Covariates(i, kk-1) * res[i];
  //     }
  //   }
  // }
  // arma::arma_print(res);
  
  // if (iter > 0)
  // {
  //     pk_mat.col(0) = pk;
  //     H_com = pk_mat.t()*Hessian*pk_mat;
  //     H_com_val = H_com(0, 0);
  //     fn_2 = arma::pow(fn, 2);
  //     prev_fn_2 = arma::pow(prev_fn, 2);
  //     bktrkg = (arma::sum(fn_2)/2 - arma::sum(prev_fn_2)/2)/(arma::dot(pk, fn) + H_com_val);
  //     if(bktrkg > eta)
  //     {
  //         printf("here\n");
  //         init_parameters = init_parameters - (1 - contraction) * pk;
  //         if (iter >= max_iter - 1) Convergence = false;
  //         continue;
  //     }
  // }
  
  // arma::arma_print(fn);
  // arma::arma_print(Hessian);
  // arma::arma_print(dNt);
  // if(iter == 0)
  //     arma::arma_print(dPhi);
  // printf("This way\n");
  
  
  // asymptotic variance
  
  
  return Rcpp::List::create(
    Named("dlam") = dLam
  );
}