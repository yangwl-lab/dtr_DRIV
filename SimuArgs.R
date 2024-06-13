source("Settings.R")
library(LongCART)
library(randomForestSRC)


SimuArg = New_SimuArg(nrep = 100, N = 1600, p = 2, p_U = 1, Scenario = "exogenous", max_t = 3,
                      theta = 0.1, 
                      unmeasured_Confounding = unmeasured_Confounding,
                      InitCovariates = InitCovariates,
                      InitAssignment = InitAssignment,
                      SurvTime = SurvTime,
                      SwitchingTime = SwitchingTime,
                      CensoringTime = CensoringTime,
                      Control = list(json_save = T),
                      gamma = rep(c(1, -1), length.out = 2),
                      alpha = rep(0.25, 3),
                      beta = rep(0.25, 3),
                      diffcoef = 0.5)
# DataGenerating(SimuArg)


SimuArg_dep = New_SimuArg(nrep = 1000, N = 1600, p = 2, p_U = 1, Scenario = "exogenous", max_t = 3,
                      theta = 0.1, 
                      unmeasured_Confounding = unmeasured_Confounding,
                      InitCovariates = InitCovariates,
                      InitAssignment = InitAssignment,
                      SurvTime = SurvTime,
                      SwitchingTime = SwitchingTime,
                      CensoringTime = CensoringTime,
                      Control = list(json_save = T,
                                     Annotation = "dep"),
                      gamma = rep(c(1, -1), length.out = 2),
                      alpha = rep(0.25, 3),
                      beta = rep(0.25, 3),
                      diffcoef = 0.5)
# DataGenerating(SimuArg_dep)


SimuArg_KangSchafer = New_SimuArg(nrep = 100, N = 1600, p = 1, p_U = 1, Scenario = "exogenous", max_t = 0.8,
                          theta = 0.5, 
                          unmeasured_Confounding = unmeasured_Confounding,
                          InitCovariates = InitCovariates,
                          InitAssignment = InitAssignment_Nonlinear,
                          SurvTime = SurvTime_Nonlinear,
                          SwitchingTime = SwitchingTime,
                          CensoringTime = CensoringTime,
                          Control = list(json_save = T,
                                         Annotation = "KangSchafer"),
                          beta = rep(0.5, 2),
                          diffcoef = 3)
# DataGenerating(SimuArg_KangSchafer)

ml_fitting_SurvCART = function(data, predictx, stime) {
  names_val = colnames(data)
  names_val = names_val[-which(names_val == "time" | names_val == "event")]
  data$patid = 1:nrow(data)
  data$event = ifelse(data$event, 1, 0)
  mod = SurvCART(data, patid = "patid", timevar = "time", censorvar = "event",
                 gvars = names_val, tgvars = rep(1, length(names_val)))
  
}


ml_fitting_rfsrc = function(data, predictx) {
  ind_control = apply(data$D_status, 1, function(d) return(all(d == 0)))
  p = ncol(data$Covariates)
  tmp_data = data.frame(time = data$time, event = data$event,
                        data$Covariates)
  colnames(tmp_data) = c("time", "event", paste0("X", 1:p))
  colnames(predictx) = paste0("X", 1:p)
  mod = rfsrc(Surv(time, event) ~ ., data = tmp_data)
  pre_mod = predict(mod, newdata = predictx)
  pre_mod_stime = pre_mod$time.interest
  pre_mod_chf = pre_mod$chf
  pre_mod_chf = cbind(pre_mod_chf[, 1], t(apply(pre_mod_chf, 1, diff)))
  tmp_chf = matrix(0, nrow = nrow(predictx), ncol = length(data$stime))
  k = 1
  for(j in 1:length(pre_mod_stime)) {
    tmp_chf[, data$stime >= pre_mod_stime[j]] = pre_mod_chf[, j]
  }
  return(tmp_chf)
}

ml_fitting_DRIV = function(data, predictx) {
  tmp_df = data.frame(IV = data$IV, data$Covariates2)
  p = ncol(data$Covariates)
  mod = glm(IV ~ ., data = tmp_df, family = binomial(link = "logit"))
  IV_c = data$IV - expit(predict(mod))
  init_parameters = runif(1 + p)
  out = integral2_est(init_parameters = init_parameters, time = data$time, 
                      event = data$event, IV = data$IV, IV_c = IV_c, Covariates = data$Covariates, 
                      D_status = data$D_status, stime = data$stime, max_iter = 50,
                      tol = 1e-5, contraction = 0.5, eta = 1e-4)
  fitting = matrix(as.matrix(predictx) %*% out$x[-1], nrow = nrow(predictx), 
                   ncol = length(data$stime), byrow = FALSE)
  return(fitting)
}




ml_fitting_propensity = function(data, predictx) {
  colnames(data) = c("IV", "X1")
  colnames(predictx) = c("IV", "X1")
  mod = rpart(IV~., data = data, method = "class")
  return(predict(mod, newdata = predictx)[, 2])
}


ml_fitting_propensity_logit = function(data, predictx){
  colnames(data) = c("IV", "X1")
  colnames(predictx) = c("IV", "X1")
  mod = glm(IV ~ ., data = data, family = binomial(link = "logit"))
  return(expit(predict(mod, newdata = predictx)))
}

trial = SimuRun(SimuArg_KangSchafer, methods = c("DRIV.cf.hz.ml.est", "DRIV.s"),
        ml_fitting_surv = ml_fitting_rfsrc,
        ml_fitting_propensity = ml_fitting_propensity_logit)
