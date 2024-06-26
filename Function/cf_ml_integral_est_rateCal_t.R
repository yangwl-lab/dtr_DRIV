library(Rcpp)
library(doParallel)
library(foreach)
sourceCpp("Function/integral/integral_customized_est.cpp")

cf_group = function(nfolds, datasize, seed) {
  cvlist <- list()
  set.seed(seed)
  n <- rep(1:nfolds, ceiling(datasize / nfolds))[1:datasize]
  temp <- sample(n, datasize)
  x <- 1:nfolds
  dataseq <- 1:datasize
  cvlist <- lapply(x, function(x) dataseq[temp == x])
  return(cvlist)
}

calculate_bias <- function(cf_surv_input, Lam, cf_surv_true, tmp_int) {
  bias = sqrt(mean(apply(tmp_int * (cf_surv_input + Lam - cf_surv_true), 1, function(d) return(sum(d)))^2))
  return(bias)
}

cfhaz_ml_integral_est_cpp_rateCal = function(init_parameters, time, event, IV, 
                                             Covariates, Covariates2, D_status, stime, ml_fitting_surv,
                                             ml_fitting_propensity,
                                             ml_fitting_surv_true,
                                             ml_fitting_propensity_true,
                                             max_iter = 20, tol = 1e-5,
                                             contraction = 0.5, eta = 1e-4, nfolds = 10, seed = 5884419, 
                                             rate_propensity, rate_hazard, true_theta, target_biases_hazard, 
                                             target_biases_prop,
                                             learning_rate = 0.1) {
  
  N = length(time)
  b_Covariates = Covariates
  cflist = cf_group(nfolds = nfolds, datasize = N, seed = seed)
  cf_IV_c = rep(0, N)
  cf_surv = matrix(0, nrow = N, ncol = length(stime))
  cf_cumsurv = matrix(0, nrow = N, ncol = length(stime))
  
  for (i in 1:length(cflist)) {
    cat("cv", i, "\n")
    ind = cflist[[i]]
    tmp_df = data.frame(IV = IV[-ind], Covariates2[-ind, , drop = FALSE])
    pred_df = data.frame(IV = IV[ind], Covariates2[ind, , drop = FALSE])
    mod = ml_fitting_propensity(tmp_df, predictx = pred_df)
    # mod = glm(IV ~ ., data = tmp_df, family = binomial(link = "logit"))
    
    tmp_df_surv = list(time = time[-ind], 
                       event = event[-ind],
                       Covariates = Covariates[-ind, , drop = FALSE],
                       Covariates2 = Covariates2[-ind, , drop = FALSE],
                       IV = IV[-ind],
                       D_status = D_status[-ind, ],
                       stime = stime)
    pred_df = data.frame(Covariates[ind, ])
    mod2 = ml_fitting_surv(tmp_df_surv, predictx = pred_df)
    
    # dLam = diff(out$Lam[, 1])
    # print(head(out$Lam))
    # if(any(is.na(dLam))) dLam[which(is.na(dLam))] = 0
    # cf_IV_c[ind] = IV[ind] - expit(predict(mod, newdata = data.frame(Covariates2[ind, ])))
    cf_IV_c[ind] = IV[ind] - mod
    cf_surv[ind, ] = mod2$surv
    cf_cumsurv[ind, ] = mod2$cumsurv
    
    # cf_dLam[ind, ] = matrix(out$Lam[, 1], byrow = T, nrow = length(ind), ncol = length(stime))
    # if(i == 1) print(matrix(out$Lam[, 1], byrow = T, nrow = length(ind), ncol = length(stime))[, 1:10])
  }
  
  mod = ml_fitting_propensity_true(Covariates2)
  cf_IV_c_true = IV - mod
  mod2_true = ml_fitting_surv_true(Covariates, stime = stime)
  cf_surv_true = mod2_true$surv
  cf_cumsurv_true = mod2_true$cumsurv
  out = list()
  k = 1
  
  int_D = t(apply(D_status, 1, function(d) return(cumsum(d * c(0, diff(stime))))))
  Y_t = sapply(stime, function(d) ifelse(time >= d, 1, 0))
  tmp_int = Y_t
  
  max_iterations = 20
  rate_bounds = c(-1, 5)
  registerDoParallel(cores = detectCores() - 1) # Use all but one core
  try({
  for (iteration in 1:max_iterations) {
    cat("Iteration", iteration, "\n")
    
    out <- foreach(j = 1:ncol(rate_propensity), .combine = "c", .multicombine = TRUE) %:%
      foreach(jj = 1:nrow(rate_hazard), .combine = "list", .multicombine = TRUE) %dopar% {
        cat("comb", j, jj, "\n")
        
        cf_surv_input = cf_surv * rate_hazard[j, jj] + cf_surv_true * (1 - rate_hazard[j, jj])
        cf_cumsurv_input = cf_cumsurv * rate_hazard[j, jj] + cf_cumsurv_true * (1 - rate_hazard[j, jj])
        cf_IV_c_input = cf_IV_c * rate_propensity[j, jj] + cf_IV_c_true * (1 - rate_propensity[j,jj])
        
        est_result <- integral_customized_est(0.1, time = time, event = event, IV = IV, IV_c = cf_IV_c_input, 
                                              D_status = D_status, stime = stime,
                                              ConfoundingPart = cf_surv_input, max_iter = max_iter,
                                              tol = tol, contraction = contraction, eta = eta)
        
        mse_propensity = sqrt(mean((cf_IV_c_input - cf_IV_c_true)^2))
        Lam = matrix(est_result$Lam[, 1], nrow = N, ncol = length(stime), byrow = TRUE)
        
        current_bias_surv = calculate_bias(cf_surv_input, Lam, cf_surv_true, tmp_int)
        current_bias_propensity = mse_propensity
        
        # Adjust the rate based on the current bias for survival and propensity models
        new_rate_hazard = rate_hazard[j, jj] - learning_rate * (current_bias_surv - target_biases_hazard[jj])
        new_rate_propensity = rate_propensity[j, jj] - (current_bias_propensity - target_biases_prop[j])
        
        # Ensure rates are within bounds
        new_rate_hazard = max(min(new_rate_hazard, rate_bounds[2]), rate_bounds[1])
        new_rate_propensity = max(min(new_rate_propensity, rate_bounds[2]), rate_bounds[1])
        
        cat("Updated rate_hazard", rate_hazard[j, jj], current_bias_surv, "Updated rate_propensity", rate_propensity[j, jj], mse_propensity,"\n")
        
        est_result$mse_propensity = mse_propensity
        est_result$mse_surv = current_bias_surv
        est_result$rate_hazard = new_rate_hazard
        est_result$rate_propensity = new_rate_propensity
        est_result$diff_propensity = (current_bias_propensity - target_biases_prop[j])
        est_result$diff_hazard = (current_bias_surv - target_biases_hazard[jj])
        est_result
      }
    results <- out
    
    # Extract and update rates
    rate_hazard <- matrix(sapply(results, function(x) x$rate_hazard), nrow = nrow(rate_hazard), ncol = ncol(rate_hazard), byrow = TRUE)
    # browser()
    rate_propensity <- matrix(sapply(results, function(x) x$rate_propensity), nrow = nrow(rate_propensity), ncol = ncol(rate_propensity), byrow = TRUE)
    
    
    # Check for convergence
    if (all(sapply(out, function(x) x$diff_hazard) < 1e-3) && all(sapply(out, function(x) x$diff_propensity) < 1e-3)) {
      cat("Converged in", iteration, "iterations\n")
      break
    }
  }
  })
  
  stopImplicitCluster()
  return(out)
}
