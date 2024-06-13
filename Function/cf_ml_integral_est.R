library(Rcpp)
# sourceCpp("Function/integral/integral_est.cpp")
sourceCpp("Function/integral/integral_customized_est.cpp")


cf_group = function(nfolds, datasize, seed) {
  cvlist <- list()
  set.seed(seed)
  n <- rep(1:nfolds,ceiling(datasize/nfolds))[1:datasize]    
  temp <- sample(n,datasize)  
  x <- 1:nfolds
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x])
  return(cvlist)
}

cfhaz_ml_integral_est_cpp = function(init_parameters, time, event, IV, 
                                  Covariates, Covariates2, D_status, stime, ml_fitting_surv,
                                  ml_fitting_propensity,
                                  max_iter = 20, tol = 1e-5,
                                  contraction = 0.5, eta = 1e-4, nfolds = 10, seed = 5884419)
{
  N = length(time)
  # b_Covariates = sapply(BasisFun, function(d) return(d(Covariates)))
  b_Covariates = Covariates
  cflist = cf_group(nfolds = nfolds, datasize = N, seed = seed)
  cf_IV_c = rep(0, N)
  cf_surv = matrix(0, nrow = N, ncol = length(stime))
  # ind_control = apply(D_status, 1, function(d) return(all(d == 0)))
  for (i in 1:length(cflist)) {
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
    out = ml_fitting_surv(tmp_df_surv, predictx = pred_df)
    
    # dLam = diff(out$Lam[, 1])
    # print(head(out$Lam))
    # if(any(is.na(dLam))) dLam[which(is.na(dLam))] = 0
    # cf_IV_c[ind] = IV[ind] - expit(predict(mod, newdata = data.frame(Covariates2[ind, ])))
    cf_IV_c[ind] = IV[ind] - mod
    cf_surv[ind, ] = out
    
    # cf_dLam[ind, ] = matrix(out$Lam[, 1], byrow = T, nrow = length(ind), ncol = length(stime))
    # if(i == 1) print(matrix(out$Lam[, 1], byrow = T, nrow = length(ind), ncol = length(stime))[, 1:10])
  }
  
  
  out = integral_customized_est(0.1, time = time, event = event, IV = IV, IV_c = cf_IV_c, 
                                D_status = D_status, stime = stime,
                                ConfoundingPart = cf_surv, max_iter = max_iter,
                                tol = tol, contraction = contraction, eta = eta)
  
  return(out)
}
