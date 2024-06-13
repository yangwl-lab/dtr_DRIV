#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;


// [[Rcpp::export]]
SEXP integral_est(
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
    
    int step = 0;
    bool Convergence;
    double tmp, tmp_cexpbetaD, tmp_cD, old_cexpbetaD, betaD, Asyvar, hatsigma;
    arma::vec beta(p), res(n), Covbeta(n), int_cexpbetaD(n), int_cdexpbetaD_1(n), int_cdexpbetaD_2(n), fn(p+1), prev_fn(p+1), fn_abs(p+1), SY(k), dSY(k), new_parameters(p+1), diff(p+1), pk(p+1), dLam(k+1);
    arma::mat int_D(n, k), dNt(n, k), Yt(n, k), int_expbetaD(n, k), int_Lam(n, k), int_dexpbetaD(n, k), int_dexpbetaD_Lam(n, k), int_expbetaD_dLam(n, k), int_expbetaD_dLam_dalpha(n, p), dPhi(n, p+1), Hessian(p+1, p+1), prev_Hessian(p+1, p+1), inv_Hessian(p+1, p+1);
    arma::mat pk_mat(p+1, 1), H_com(1, 1);
    // double H_com_val, bktrkg;
    arma::vec fn_2(p+1), prev_fn_2(p+1);
    
    for (int iter = 0; iter < max_iter; iter++)
    {
        res.zeros(); Covbeta.zeros(); int_cexpbetaD.zeros(); int_Lam.zeros(); int_cdexpbetaD_1.zeros(); int_cdexpbetaD_2.zeros();
        fn.zeros(); SY.zeros(); dSY.zeros(); int_D.zeros(); int_expbetaD.zeros(), int_dexpbetaD.zeros(), int_dexpbetaD_Lam.zeros();
        int_expbetaD_dLam_dalpha.zeros();dPhi.zeros(); Hessian.zeros();
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
            }
        }


        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < k; j++)
            {
                for (int kk = 0; kk < p; kk++)
                {
                    int_expbetaD_dLam_dalpha(i, kk) = 0;
                }
                
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
                    }
                    

                    for (int ii = 0; ii < n; ii++)
                    {
                        tmp_cexpbetaD = exp((int_D(i, j) + int_D(ii, j)) * betaD);
                        if (j == 0)
                        {
                            // printf("%d\n", j);
                            tmp_cD = IV[i] + IV[ii];
                            if (betaD != 0)
                            {
                                int_cexpbetaD[ii] = ((tmp_cD) > 0.5) ? (tmp_cexpbetaD - 1) / ((tmp_cD) * betaD):stime[j];
                                int_cdexpbetaD_1[ii] = (IV[i] > 0.5) ? (int_D(i, j) * tmp_cexpbetaD/((tmp_cD)*betaD) - int_cexpbetaD[ii]/((tmp_cD)*betaD)):0;
                                int_cdexpbetaD_2[ii] = (IV[ii] > 0.5) ? (int_D(ii, j) * tmp_cexpbetaD/((tmp_cD)*betaD) - int_cexpbetaD[ii]/((tmp_cD)*betaD)):0;
                            }
                            else
                            {
                                int_cexpbetaD[ii] = stime[j];
                                int_cdexpbetaD_1[ii] = int_D(i, j) * int_D(i, j)/2;
                                int_cdexpbetaD_2[ii] = int_D(ii, j) * int_D(ii, j)/2;
                            }
                        }
                        else
                        {
                            
                            old_cexpbetaD = exp((int_D(i, j-1) + int_D(ii, j-1)) * betaD);
                            tmp_cD = D_status(i, j-1) + D_status(ii, j-1);
                            if (betaD != 0)
                            {
                                int_cexpbetaD[ii] = ((tmp_cD) > 0.5) ? 
                                                    (tmp_cexpbetaD - old_cexpbetaD) / ((tmp_cD) * betaD):
                                                    (stime[j] - stime[j-1]) * exp((int_D(i, j) + int_D(ii, j))*betaD);
                                int_cdexpbetaD_1[ii] = (D_status(i, j-1) > 0.5) ? 
                                                    ((int_D(i, j) * tmp_cexpbetaD - int_D(i, j-1) * old_cexpbetaD)/((tmp_cD)*betaD) - 
                                                    int_cexpbetaD[ii]/((tmp_cD)*betaD)):int_D(i, j-1) * int_cexpbetaD[ii];
                                int_cdexpbetaD_2[ii] = (D_status(ii, j-1) > 0.5) ? 
                                                    ((int_D(ii, j) * tmp_cexpbetaD - int_D(ii, j-1) * old_cexpbetaD)/((tmp_cD)*betaD) - 
                                                    int_cexpbetaD[ii]/((tmp_cD)*betaD)):int_D(ii, j-1) * int_cexpbetaD[ii];
                            }
                            else
                            {
                                // printf("%d\n", j);
                                int_cexpbetaD[ii] = stime[j] - stime[j-1];
                                int_cdexpbetaD_1[ii] = (stime[j] - stime[j-1]) * int_D(i, j-1) + D_status(i, j-1) * (stime[j] - stime[j-1]) * (stime[j] - stime[j-1])/2;
                                int_cdexpbetaD_2[ii] = (stime[j] - stime[j-1]) * int_D(ii, j-1) + D_status(ii, j-1) * (stime[j] - stime[j-1]) * (stime[j] - stime[j-1])/2;
                            }
                        }
                        tmp = dNt(ii, j) * tmp_cexpbetaD - Yt(ii, j) * int_cexpbetaD[ii] * (D_status(ii, j) * betaD + Covbeta[ii]);
                        int_Lam(i, j) = int_Lam(i, j) +  tmp;
                        int_dexpbetaD_Lam(i, j)  = int_dexpbetaD_Lam(i, j) + dNt(ii, j) * int_D(i, j) * tmp_cexpbetaD - 
                                                        Yt(ii, j) * int_cdexpbetaD_1[ii] * (D_status(ii, j) * betaD + Covbeta[ii]);
                        int_expbetaD_dLam(i, j) = int_expbetaD_dLam(i, j) + (dNt(ii, j) * int_D(ii, j) * tmp_cexpbetaD - 
                                                        Yt(ii, j) * int_cdexpbetaD_2[ii] * (D_status(ii, j) * betaD + Covbeta[ii]) - 
                                                        int_cexpbetaD[ii] * Yt(ii, j) * D_status(ii, j)) * SY[j] - tmp * dSY[j];                        
                        for (int kk = 0; kk < p; kk++)
                        {
                            int_expbetaD_dLam_dalpha(i, kk) = int_expbetaD_dLam_dalpha(i, kk) - int_cexpbetaD[ii] * Yt(ii, j) * Covariates(ii, kk)/SY[j];
                        }
                    }
                    
                    int_Lam(i, j) = int_Lam(i, j) / SY[j];
                    int_dexpbetaD_Lam(i, j) = int_dexpbetaD_Lam(i, j) / SY[j];
                    int_expbetaD_dLam(i, j) = int_expbetaD_dLam(i, j) / (SY[j] * SY[j]);

                }
                res[i] = res[i] + exp(int_D(i, j) * betaD) * dNt(i, j) - Yt(i, j) * int_expbetaD(i, j) * (D_status(i, j) * betaD + Covbeta[i]) - Yt(i, j) * int_Lam(i, j);

                for (int kk = 0; kk < p+1; kk++)
                {
                    if (kk == 0)
                    {
                        
                        dPhi(i, kk) = dPhi(i, kk) + dNt(i, j) * int_D(i, j) * exp(betaD * int_D(i, j)) - 
                                        Yt(i, j) * int_dexpbetaD(i, j) * (D_status(i, j) * betaD + Covbeta[i]) - 
                                        Yt(i, j) * int_dexpbetaD_Lam(i, j) - 
                                        Yt(i, j) * int_expbetaD(i, j) * D_status(i, j) - 
                                        Yt(i, j) * int_expbetaD_dLam(i, j);
                    }
                    else
                    {
                        dPhi(i, kk) = dPhi(i, kk) - Yt(i, j) * Covariates(i, kk-1) * int_expbetaD(i, j) - 
                                        Yt(i, j) * int_expbetaD_dLam_dalpha(i, kk-1);
                    }
                    
                }
                
            }
            // printf("%d\n", i);
            for (int kk = 0; kk < p+1; kk++)
            {
                if (kk == 0)
                {
                    fn[kk] = fn[kk] + IV_c[i] * res[i];
                    
                    for (int kkk = 0; kkk < p+1; kkk++)
                    {
                        Hessian(kk, kkk) = Hessian(kk, kkk) + IV_c[i] * dPhi(i, kkk);                        
                    }
                }
                else
                {
                    fn[kk] = fn[kk] + Covariates(i, kk-1) * res[i];
                    for (int kkk = 0; kkk < p+1; kkk++)
                    {
                        // if(iter == 0) arma::arma_print(Covariates(i, kk-1) * dPhi(i, kkk));
                        Hessian(kk, kkk) = Hessian(kk, kkk) + Covariates(i, kk-1) * dPhi(i, kkk);
                    }
                }
            }
        }
        
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
        pk = -arma::solve(Hessian, fn);
        new_parameters = init_parameters + pk;
        fn_abs = arma::abs(fn);
        diff = arma::abs(new_parameters - init_parameters);
        if(arma::sum(fn_abs) < tol || arma::sum(diff) < (tol * tol))
        {
            step = iter;
            Convergence = true;
            break;
        }
        init_parameters = new_parameters;
        // prev_fn = fn;
        if (iter >= (max_iter - 1))
        {
            step = iter;
            Convergence = false;
        }
    }
    
    // asymptotic variance
    Asyvar = 0.0;
    hatsigma = 0.0;
    
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < k; j++)
        {
            hatsigma = hatsigma + ((IV_c[i] * IV_c[i])*exp(2*int_D(i, j)*betaD)*dNt(i, j));
        }
    }
    Asyvar = hatsigma / (Hessian(0, 0)*Hessian(0, 0));
    
    
    dLam.zeros();
    betaD = init_parameters[0];
    beta = init_parameters.subvec(1, p);
    Covbeta.zeros();
    for (int i = 0; i < n; i++)
    {
        for (int ii = 0; ii < p; ii++)
        {
            Covbeta[i] += Covariates(i, ii) * beta[ii];
        }
    }
    
    for (int j = 0; j < k+1; j++)
    {
        for (int i = 0; i < n; i++)
        {
            if (j == 0)
            {
                dLam[j] = dLam[j] + (-IV[i]*betaD-Covbeta[i]) / n;
            } else 
            {
                dLam[j] = dLam[j] + exp(int_D(i, j-1) * betaD) * Yt(i, j-1) * (dNt(i, j-1) -D_status(i, j-1)*betaD - Covbeta[i]) / SY[j-1];
            }
        }
    }

    return Rcpp::List::create(
        Named("x") = init_parameters,
        Named("Convergence") = Convergence,
        Named("iter") = step,
        Named("var") = Asyvar,
        Named("Lam") = dLam
    );
}
