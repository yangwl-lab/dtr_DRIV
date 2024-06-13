library(Rcpp)
sourceCpp("Function/integral/integral2.cpp")

integral2_cpp_F  = function(x, time = time, 
                            event = event, IV = IV, IV_c = IV_c, Covariates = Covariates, 
                            D_status = D_status, stime = stime, max_iter = max_iter,
                            tol = tol, contraction = contraction, eta = eta){
  out = integral2(x, time = time, 
                  event = event, IV = IV, IV_c = IV_c, Covariates = Covariates, 
                  D_status = D_status, stime = stime, max_iter = max_iter,
                  tol = tol, contraction = contraction, eta = eta)
  # print(out$fn)
  # print(out$dlam)
  return(out$fn)
}

integral2_cpp = function(init_parameters, time, event, IV, 
                            Covariates, Covariates2, D_status, stime, max_iter = 20, tol = 1e-5,
                            contraction = 0.5, eta = 1e-4)
{
  mod = glm(IV ~ Covariates2, family = binomial(link = "logit"))
  IV_c = IV - expit(predict(mod))
  # print(head(IV_c))
  
  out = nleqslv::nleqslv(init_parameters, integral2_cpp_F, time = time, 
                         event = event, IV = IV, IV_c = IV_c, Covariates = Covariates, 
                         D_status = D_status, stime = stime, max_iter = max_iter,
                         tol = tol, contraction = contraction, eta = eta, method = "Newton",
                         global = "none", jacobian = TRUE, control = list(maxit = max_iter,
                                                                          trace = 1))
  return(out)
}


