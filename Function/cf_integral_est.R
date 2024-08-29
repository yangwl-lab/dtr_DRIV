library(Rcpp)
# sourceCpp("Function/integral/integral_est.cpp")
sourceCpp("Function/integral/cf_integral_est.cpp")


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

cfhaz_integral_est_cpp = function(init_parameters, time, event, IV, 
                                  Covariates, Covariates2, D_status, stime,
                                  max_iter = 20, tol = 1e-5,
                                  contraction = 0.5, eta = 1e-4, nfolds = 10, seed = 5884419)
{
  N = length(time)
  # b_Covariates = sapply(BasisFun, function(d) return(d(Covariates)))
  b_Covariates = Covariates
  cflist = cf_group(nfolds = nfolds, datasize = N, seed = seed)
  cf_IV_c = rep(0, N)
  cf_beta = matrix(0, nrow = N, ncol = dim(b_Covariates)[2], byrow = T)
  cf_dLam = matrix(0, nrow = N, ncol = length(stime), byrow = T)
  for (i in 1:length(cflist)) {
    ind = cflist[[i]]
    tmp_df = data.frame(IV = IV[-ind], Covariates2[-ind, , drop = FALSE])
    mod = glm(IV ~ ., data = tmp_df, family = binomial(link = "logit"))
    IV_c = IV[-ind] - expit(predict(mod))
    
    out = integral2_est(init_parameters = init_parameters, time = time[-ind], 
                        event = event[-ind], IV = IV[-ind], IV_c = IV_c, Covariates = b_Covariates[-ind, , drop = FALSE], 
                        D_status = D_status[-ind, ], stime = stime, max_iter = max_iter,
                        tol = tol, contraction = contraction, eta = eta)
    # dLam = diff(out$Lam[, 1])
    # print(head(out$Lam))
    # if(any(is.na(dLam))) dLam[which(is.na(dLam))] = 0
    cf_IV_c[ind] = IV[ind] - expit(predict(mod, newdata = data.frame(Covariates2[ind, , drop = FALSE])))
    # cf_IV_c[ind] = IV[ind]
    cf_beta[ind, ] = matrix(out$x[-1], byrow = T, nrow = length(ind), ncol = length(out$x[-1]))
    
    cf_dLam[ind, ] = matrix(out$dLam[, 1], byrow = T, nrow = length(ind), ncol = length(stime))
    # if(i == 1) print(matrix(out$Lam[, 1], byrow = T, nrow = length(ind), ncol = length(stime))[, 1:10])
  }
  
  
  out = cf_integral_est(0.01, time = time, event = event, IV = IV, cf_IV_c = cf_IV_c,
                         Covariates = b_Covariates, D_status = D_status, stime = stime,
                         cf_beta = cf_beta, cf_dLam = cf_dLam, max_iter = max_iter,
                         tol = tol, contraction = contraction, eta = eta)
  
  return(out)
}
