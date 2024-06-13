source("ClassBuilder.R")
source("Generic.R")
source("Functions.R")
source("Valid.R")
source("Function/integral2_est.R")
source("Function/integral2.R")
source("Function/integral.R")
source("Function/cf_integral.R")
source("Function/cf_integral_est.R")


InitCovariates = function(N, p){
  return(matrix(runif(N*p), ncol = p, nrow = N))
}

InitAssignment = function(N, p, Covariates, gamma){
  Z_p = expit(Covariates[, 1:(p-1)] %*% gamma[1:p-1])
  return(Z = rbinom(N, 1, Z_p))
}

SurvTime = function(N, p, Covariates, W, theta, alpha, Z){
  W_copy = W
  W[W <= 0] = Inf
  W_copy[W_copy <= 0] = 0
  # T_D = nleqslv(rep(0, n), F_time, w = W, z = Z, Cov = Cov, uni = runif(n))$x
  T_D = rexp(N)/(0.25 + theta*Z + Covariates %*% alpha)
  T_D_ind = T_D >= W
  if(any(T_D_ind)){
    T_D[T_D_ind] = ((T_D*(0.25 + theta*Z + Covariates %*% alpha)- 
                       (theta * Z * W_copy - theta* (1-Z)*W_copy))/
                      (0.25 + theta*(1-Z) + Covariates %*% alpha))[T_D_ind]
  }
  
  return(as.vector(T_D))
}

SwitchingTime = function(N, p, Covariates, Z, beta){
  W = ifelse(as.logical(Z), rexp(N)/(0.5*(0.1 + Z + Covariates %*% beta)),
             rexp(N)/(0.5*(0.1 - (1 - Z) + Covariates %*% beta)))
  return(W)
}

CensoringTime = function(N, p, Covariates) {
  C = rexp(N)/(0.1 + Covariates[, 1:(p-1)] %*% rep(0.05, p-1))
}


SimuArg = New_SimuArg(nrep = 100, N = 1600, p = 3, Scenario = "exogenous", max_t = 3,
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


dat = jsonlite::read_json("DataGenerated/1600/1.json", simplifyVector = T)

ModelPar = New_ModelPar(N = 1600, p = 3, max_t = 6, dat = dat, method = "DRIV.s2",
                        Control = list(init_parameters = c(0.1, 0.25, 0.25),
                                       max_iter = 1))
DataFitting(ModelPar)


results = SimuRun(SimuArg, methods = c("ITT", "DRIV.cf.hz.est"),
        Control = list(grid = 100))
apply(results$SimuResults$DRIV.s$Coef, 1, mean)
apply(results$SimuResults$DRIV.s$Coef, 1, sd)
mean(sqrt(results$SimuResults$DRIV.s$Var))
mean(results$SimuResults$DRIV.s$Coef[1, ] + 1.96*sqrt(results$SimuResults$DRIV.s$Var)>= 0.1 &
       results$SimuResults$DRIV.s$Coef[1, ] - 1.96*sqrt(results$SimuResults$DRIV.s$Var)<= 0.1)

apply(results$SimuResults$DRIV.cf.hz.est$Coef, 1, mean)[1]
apply(results$SimuResults$DRIV.cf.hz.est$Coef, 1, sd)
mean(sqrt(results$SimuResults$DRIV.cf.hz.est$Var))
mean(results$SimuResults$DRIV.cf.hz.est$Coef[1, ] + 1.96*sqrt(results$SimuResults$DRIV.cf.hz.est$Var)>= 0.1 &
       results$SimuResults$DRIV.cf.hz.est$Coef[1, ] - 1.96*sqrt(results$SimuResults$DRIV.cf.hz.est$Var)<= 0.1)

result = TRTSWE(dat, max_t = 6, methods = c("DRIV"), 
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


