InitCovariates = function(N, p){
  L = matrix(runif(N*(p-1)), ncol = p-1, nrow = N)
  U = exp(L[, p-1])*rexp(N)
  return(cbind(L, U))
}

InitAssignment = function(N, p, Covariates, gamma){
  Z_p = expit(Covariates[, 1:(p-1)] %*% gamma[1:p-1])
  return(Z = rbinom(N, 1, Z_p))
}

SurvTime = function(N, p, Covariates, W, theta, alpha, Z, r_T_D){
  W_copy = W
  W[W <= 0] = Inf
  W_copy[W_copy <= 0] = 0
  # T_D = nleqslv(rep(0, n), F_time, w = W, z = Z, Cov = Cov, uni = runif(n))$x
  T_D = r_T_D/(0.25 + theta*Z + Covariates %*% alpha)
  T_D_ind = T_D >= W
  if(any(T_D_ind)){
    T_D[T_D_ind] = ((T_D*(0.25 + theta*Z + Covariates %*% alpha)- 
                       (theta * Z * W_copy - theta* (1-Z)*W_copy))/
                      (0.25 + theta*(1-Z) + Covariates %*% alpha))[T_D_ind]
  }
  
  return(as.vector(T_D))
}

SwitchingTime = function(N, p, Covariates, Z, beta, r_T_D, theta){
  T0 = r_T_D/(0.25 + Covariates %*% beta)
  T1 = r_T_D/(0.25 + theta + Covariates %*% beta)
  W = ifelse(as.logical(Z), rexp(N)/(exp(-T1)*(0.1 + Z + Covariates %*% beta)),
             rexp(N)/(exp(-T0)*(0.1 - (1 - Z) + Covariates %*% beta)))
  return(W)
}

CensoringTime = function(N, p, Covariates) {
  C = rexp(N)/(0.1 + Covariates[, 1:(p-1)] %*% rep(0.05, p-1))
}


SimuArg = New_SimuArg(nrep = 1000, N = 1600, p = 3, Scenario = "endogenous", max_t = 3,
                      theta = 0.1, InitCovariates = InitCovariates,
                      InitAssignment = InitAssignment,
                      SurvTime = SurvTime,
                      SwitchingTime = SwitchingTime,
                      CensoringTime = CensoringTime,
                      Control = list(json_save = T, Annotation = "S2"),
                      gamma = rep(c(1, -1), length.out = 2),
                      alpha = rep(0.25, 3),
                      beta = rep(0.25, 3))

DataGenerating(SimuArg)


dat = jsonlite::read_json("DataGenerated/800S2/1.json", simplifyVector = T)
ModelPar = New_ModelPar(N = 800, p = 9, max_t = 6, dat = dat, method = "ITT")
DataFitting(ModelPar)
result = TRTSWE(dat, max_t = 6, methods = c("ITT"), 
       Covariates2 = dat$Covariates)


results = SimuRun(SimuArg, methods = c("ITT", "remove.control"),
                  Control = list(grid = 100))
apply(results$SimuResults$ITT$Coef, 1, mean)
apply(results$SimuResults$remove.control$Coef, 1, mean)
apply(results$SimuResults$DRIV.s$Coef, 1, mean)
apply(results$SimuResults$DRIV.s$Coef, 1, sd)
mean(sqrt(results$SimuResults$DRIV.s$Var))

interval = cbind(results$SimuResults$DRIV.s$Coef[1, ] - 1.96*sqrt(results$SimuResults$DRIV.s$Var),
  results$SimuResults$DRIV.s$Coef[1, ] + 1.96*sqrt(results$SimuResults$DRIV.s$Var))

sum(interval[, 1] <= 0.1 & interval[, 2] >=0.1)/1000



