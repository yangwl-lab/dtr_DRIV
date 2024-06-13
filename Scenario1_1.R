source("ClassBuilder.R")
source("Generic.R")
source("Functions.R")
source("Valid.R")
source("Function/integral2_est.R")
source("Function/integral2.R")
source("Function/integral.R")
source("Function/cf_integral.R")
source("Function/cf_integral_est.R")
source("Function/integral_customized.R")


InitCovariates = function(N, p){
  L = matrix(runif(N*(p)), ncol = p, nrow = N)
  U = rexp(N)
  return(cbind(L, U))
}

InitCovariates2 = function(N, p){
  L = matrix(runif(N*(p)), ncol = p, nrow = N)
  U = exp(L[, p])
  return(cbind(L, U))
}


InitCovariates3 = function(N, p){
  L = matrix(runif(N*(p)), ncol = p, nrow = N)
  U = L[, p]*rnorm(N, sd = 0.5)
  return(cbind(L, U))
}




InitAssignment = function(N, p, Covariates, gamma){
  Z_p = expit(Covariates[, 1:(p)] %*% gamma[1:p])
  return(Z = rbinom(N, 1, Z_p))
}

InitAssignment2 = function(N, p, Covariates, gamma){
  X = qnorm(Covariates[, 1:(p)])
  x1 = exp(X[, 1]/2)
  x2 = X[, 2]/(1 + exp(X[, 1]))+1
  tmp_cov = cbind(x1, x2)
  Z_p = expit(tmp_cov %*% gamma[1:p])
  return(Z = rbinom(N, 1, Z_p))
}

SurvTime = function(N, p, Covariates, W, theta, alpha, Z){
  W_copy = W
  W[W <= 0] = Inf
  W_copy[W_copy <= 0] = 0
  # T_D = nleqslv(rep(0, n), F_time, w = W, z = Z, Cov = Cov, uni = runif(n))$x
  T = rexp(N)
  T_0 = T/(0.25 + Covariates %*% alpha)
  T_D = T/(0.25 + theta*Z + Covariates %*% alpha)
  T_D_ind = T_D >= W
  if(any(T_D_ind)){
    T_D[T_D_ind] = ((T_D*(0.25 + theta*Z + Covariates %*% alpha)- 
                       (theta * Z * W_copy - theta* (1-Z)*W_copy))/
                      (0.25 + theta*(1-Z) + Covariates %*% alpha))[T_D_ind]
  }
  
  return(list(T_D = as.vector(T_D),
              T_0 = T_0))
}

SwitchingTime = function(N, p, Covariates, Z, beta){
  W = ifelse(as.logical(Z), rexp(N)/(0.5*(0.1 + Z + Covariates %*% beta)),
             rexp(N)/(0.5*(0.1 - (1 - Z) + Covariates %*% beta)))
  return(W)
}

CensoringTime = function(N, p, Covariates) {
  C = rexp(N)/(0.1 + Covariates[, 1:(p)] %*% rep(0.05, p))
}


SimuArg = New_SimuArg(nrep = 1000, N = 1600, p = 2, p_U = 1, Scenario = "exogenous", max_t = 3,
                      theta = 0.1, InitCovariates = InitCovariates,
                      InitAssignment = InitAssignment,
                      SurvTime = SurvTime,
                      SwitchingTime = SwitchingTime,
                      CensoringTime = CensoringTime,
                      Control = list(json_save = T),
                      gamma = rep(c(1, -1), length.out = 2),
                      alpha = rep(0.25, 3),
                      beta = rep(0.25, 3))

DataGenerating(SimuArg)

SimuArg2 = New_SimuArg(nrep = 1000, N = 1600, p = 2, p_U = 1, Scenario = "exogenous", max_t = 3,
                      theta = 0.1, InitCovariates = InitCovariates2,
                      InitAssignment = InitAssignment,
                      SurvTime = SurvTime,
                      SwitchingTime = SwitchingTime,
                      CensoringTime = CensoringTime,
                      Control = list(json_save = T,
                                     Annotation = "Correlated"),
                      gamma = rep(c(0.1, 0.1), length.out = 2),
                      alpha = c(0.25, 0.25, 0.75),
                      beta = rep(0.25, 3))
DataGenerating(SimuArg2)

SimuArg3 = New_SimuArg(nrep = 1000, N = 1600, p = 2, p_U = 1, Scenario = "exogenous", max_t = 3,
                       theta = 0.1, InitCovariates = InitCovariates3,
                       InitAssignment = InitAssignment,
                       SurvTime = SurvTime,
                       SwitchingTime = SwitchingTime,
                       CensoringTime = CensoringTime,
                       Control = list(json_save = T,
                                      Annotation = "Correlated2"),
                       gamma = rep(c(0.1, 0.1), length.out = 2),
                       alpha = c(0.25, 0.25, 0.25),
                       beta = rep(0.25, 3))
DataGenerating(SimuArg3)

dat = jsonlite::read_json("DataGenerated/1600Correlated2/3.json", simplifyVector = T)
# 
# ModelPar = New_ModelPar(N = 1600, p = 3, max_t = 6, dat = dat, method = "DRIV.s2",
#                         Control = list(init_parameters = c(0.1, 0.25, 0.25),
#                                        max_iter = 1))
# DataFitting(ModelPar)


results = SimuRun(SimuArg, methods = c("ITT", "remove.control", "counterfactual", "DRIV.s"),
                  Control = list(grid = 100))
results2 = SimuRun(SimuArg2, methods = c("ITT", "remove.control", "counterfactual", "DRIV.s"),
                  Control = list(grid = 100))
results3 = SimuRun(SimuArg3, methods = c("ITT", "remove.control", "counterfactual", "DRIV.s"),
                  Control = list(grid = 100,
                                 init_parameters = c(0.1, 0.25, 0.25)))


results4 = list()
results4$Coef = vector(length = 1000)
results4$Var = vector(length = 1000)
results4$Convergence = vector(length = 1000)
for(j in 1:1000) {
  dat = jsonlite::read_json(paste0("DataGenerated/1600/", j, ".json"), simplifyVector = T)
  fit_set = sample(1:1600, 1600*0.5)
  dat1 = list(Covariates = dat$Covariates[fit_set, ],
              Z = dat$Z[fit_set],
              W = dat$W[fit_set],
              T_D_c = dat$T_D_c[fit_set, ],
              T_D = dat$T_D[fit_set],
              T_0 = dat$T_0[fit_set, ],
              C = dat$C[fit_set, ],
              event = dat$event[fit_set, ])
  mod = TRTSWE(dat1, max_t = 3, methods = "DRIV.s")
  mod_IV = glm(dat1$Z ~dat1$Covariates, family = binomial(link = "logit"))
  dat2 = list(Covariates = dat$Covariates[-fit_set, ],
              Z = dat$Z[-fit_set],
              W = dat$W[-fit_set],
              T_D_c = dat$T_D_c[-fit_set, ],
              T_D = dat$T_D[-fit_set],
              T_0 = dat$T_0[-fit_set, ],
              C = dat$C[-fit_set, ],
              event = dat$event[-fit_set, ])
  stime = sort(dat2$T_D_c)
  stime = unique(stime)
  PredictedPropensityScore = expit(predict.glm(mod_IV, newdata = as.data.frame(dat2$Covariates)))
  PredictedCovariates = matrix(dat2$Covariates%*% mod$DRIV.s$Estim$Coef[2:4, ], 
                               nrow = 1600, ncol = length(stime), byrow = FALSE)
  out = TRTSWE(dat2, max_t = 3, methods = "customized", 
                  PredictedPropensityScore = PredictedPropensityScore,
                  PredictedCovariates = PredictedCovariates, 
                  Control = list(init_parameters = 0.1))
  results4$Coef[j] = out$customized$Estim$Coef
  results4$Var[j] = out$customized$Estim$Var
  results4$Convergence[j] = out$customized$Estim$Convergence
  cat("[[ rep ", j, "\t", out$customized$Estim$Coef, "\t", 
      out$customized$Estim$Var, "\t", out$customized$Estim$Convergence, "]]\n")
}
mean(results4$Coef)
sd(results4$Coef)
mean(sqrt(results4$Var))
mean((0.1 < results4$Coef + 1.96*sqrt(results4$Var) &
        0.1 > results4$Coef - 1.96*sqrt(results4$Var)))



print(results, Comp_parameters = c(0, 0.25,0.25))
apply(results$SimuResults$DRIV.s.control$Coef, 1, mean)- c(0.1, 0.25, 0.25)
apply(results$SimuResults$remove.control$Coef, 1, mean)- c(0.1, 0.25, 0.25)
apply(results$SimuResults$DRIV.s$Coef, 1, sd)
mean(sqrt(results$SimuResults$DRIV.s$Var))
mean(results$SimuResults$DRIV.s$Coef[1, ] + 1.96*sqrt(results$SimuResults$DRIV.s$Var)>= 0.1 &
       results$SimuResults$DRIV.s$Coef[1, ] - 1.96*sqrt(results$SimuResults$DRIV.s$Var)<= 0.1)


result = TRTSWE(dat, max_t = 6, methods = c("DRIV.s"), 
                Control = list(init_parameters = c(0.1, 0.25, 0.25),
                               max_iter =20, 
                               tol = 1e-5))
system.time(out <- TRTSWE(dat, max_t = 3, methods = c("DRIV.cf.hz.est", "DRIV.s", "DRIV.cf.hz"), 
                          Control = list(init_parameters = c(0.1, 0.25, 0.25, 0.25),
                                         max_iter =20, 
                                         tol = 1e-5,
                                         grid = 1000)))
apply(results$DRIV$Coef[, results$DRIV$Convergence], 1, mean)
apply(results$ITT$Coef, 1, mean)
apply(results$remove$Coef, 1, mean)

apply(results$DRIV$Coef[, results$DRIV$Convergence], 1, sd)
apply(results$ITT$Coef, 1, sd)
apply(results$remove$Coef, 1, sd)

mean(sqrt(results$DRIV$Var[results$DRIV$Convergence]))
apply(sqrt(results$ITT$Var), 1, mean)
apply(sqrt(results$remove$Var), 1, mean)

k1 = cbind(results$DRIV$Coef[1, results$DRIV$Convergence] - 1.96*sqrt(results$DRIV$Var[results$DRIV$Convergence]),
           results$DRIV$Coef[1, results$DRIV$Convergence] + 1.96*sqrt(results$DRIV$Var[results$DRIV$Convergence]))
mean(k1[, 1] <= 0.1 & k1[, 2] >= 0.1)

k2 = cbind(results$ITT$Coef[1, ] - 1.96*sqrt(results$ITT$Var[1, ]),
           results$ITT$Coef[1, ] + 1.96*sqrt(results$ITT$Var[1, ]))
mean(k2[, 1] <= 0.1 & k2[, 2] >= 0.1)


k3 = cbind(results$remove$Coef[1, ] - 1.96*sqrt(results$remove$Var[1, ]),
           results$remove$Coef[1, ] + 1.96*sqrt(results$remove$Var[1, ]))
mean(k3[, 1] <= 0.1 & k3[, 2] >= 0.1)


