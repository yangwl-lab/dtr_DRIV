library(Rcpp)
sourceCpp("Function/integral/integral2_dlam_deBug.cpp")




integral2_dLam_cpp = function(init_parameters, time, event, IV, 
                            Covariates, Covariates2, D_status, stime, max_iter = 20, tol = 1e-5,
                            contraction = 0.5, eta = 1e-4)
{
  mod = glm(IV ~ Covariates2, family = binomial(link = "logit"))
  IV_c = IV - expit(predict(mod))
  
  out = integral2_dLam(init_parameters = init_parameters, time = time, 
                     event = event, IV = IV, IV_c = IV_c, Covariates = Covariates, 
                     D_status = D_status, stime = stime, max_iter = max_iter,
                     tol = tol, contraction = contraction, eta = eta)
  
  return(out)
}
