library(Rcpp)
sourceCpp("Function/integral/integral_customized_est.cpp")


integral_customized_est_cpp = function(init_parameters, time, event, IV, IV_c, 
                             ConfoundingPart, D_status, stime, max_iter = 50, tol = 1e-5,
                             contraction = 0.5, eta = 1e-4)
{
  out = integral_customized_est(init_parameters = init_parameters, time = time, 
                      event = event, IV = IV, IV_c = IV_c, ConfoundingPart = ConfoundingPart, 
                      D_status = D_status, stime = stime, max_iter = max_iter,
                      tol = tol, contraction = contraction, eta = eta)
  return(out)
}
