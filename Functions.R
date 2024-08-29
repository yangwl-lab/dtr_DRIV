SimuRun = function(SimuArg, methods, 
                   Control = list(grid = 100, 
                                  max_iter = 20, 
                                  tol = 1e-5, 
                                  contraction = 0.5, 
                                  eta = 1e-4), ...){
  Validate_method(methods)
  args = rlang::dots_list(N = SimuArg$initials$N,
                          p = SimuArg$initials$p,
                          max_t = SimuArg$initials$max_t,
                          Control = Control, !!!list(...), .homonyms = "first")
  # attach out to SimuArg, then construct summary.SimuArg
  out = vector("list", length(methods)) 
  names(out) = methods
  templist = list(Coef = matrix(0, nrow = args$p+1, ncol = SimuArg$initials$nrep),
                  Var = matrix(0, nrow = args$p+1, ncol = SimuArg$initials$nrep))
  for (kk in methods) {
    if(kk %in% c("DRIV", "DRIV.s", "DRIV.s2", "DRIV.cf.hz.est", "DRIV.s.control", "DRIV.cf.hz.ml.est")){
      out[[kk]] = list(Coef = matrix(0, nrow = args$p+1, ncol = SimuArg$initials$nrep),
                       Var = rep(0, SimuArg$initials$nrep),
                       Convergence = rep(FALSE, SimuArg$initials$nrep))
      next
    }
    out[[kk]] = templist
  }
  for (i in 1:SimuArg$initials$nrep) {
    data = jsonlite::read_json(paste0(SimuArg$Control$save_path, 
                                      "DataGenerated/",
                                      SimuArg$initials$N,
                                      SimuArg$Control$Annotation, 
                                      "/", i, ".json"), simplifyVector = T)
    data$Covariates = data$Covariates[, 1:args$p, drop = FALSE]
    colnames(data$Covariates) = paste0("X", 1:args$p)
    args$dat = data
    cat("[[rep", i)
    for (kk in methods) {
      args$method = kk
      ModelPar = do.call(New_ModelPar, args)
      mod = DataFitting(ModelPar)
      out[[kk]]$Coef[, i] = mod$Coef
      cat("\t", out[[kk]]$Coef[1, i])
      if (kk %in% c("DRIV", "DRIV.s", "DRIV.s2", "DRIV.cf.hz.est", "DRIV.s.control", "DRIV.cf.hz.ml.est")){
        if(is.null(mod$Var)) next
        out[[kk]]$Var[i] = mod$Var
        out[[kk]]$Convergence[i] = mod$Convergence
        cat("\t", out[[kk]]$Var[i])
        cat("\t", out[[kk]]$Convergence[i])
        next
      }
      out[[kk]]$Var[, i] = mod$Var
    }
    cat("]]\n")
  }
  SimuArg$methods = methods
  SimuArg$SimuResults = out
  structure(SimuArg, class = "SimuResults")
}


SimuRun_rateCal = function(SimuArg, methods = "DRIV.cf.hz.ml.est.rateCal",
                   Control = list(grid = 100, 
                                  max_iter = 20, 
                                  tol = 1e-5, 
                                  contraction = 0.5, 
                                  eta = 1e-4), sequence,...){
  Validate_method(methods)
  args = rlang::dots_list(N = SimuArg$initials$N,
                          p = SimuArg$initials$p,
                          max_t = SimuArg$initials$max_t,
                          Control = Control, !!!list(...), .homonyms = "first")
  rate_propensity = args$rate_propensity
  rate_hazard = args$rate_hazard
  # attach out to SimuArg, then construct summary.SimuArg
  out = vector("list", length(methods)) 
  names(out) = methods
  templist = list(Coef = matrix(0, nrow = args$p+1, ncol = SimuArg$initials$nrep),
                  Var = matrix(0, nrow = args$p+1, ncol = SimuArg$initials$nrep))
  for (kk in methods) {
    if(kk %in% c("DRIV", "DRIV.s", "DRIV.s2", "DRIV.cf.hz.est", "DRIV.s.control", "DRIV.cf.hz.ml.est")){
      out[[kk]] = list(Coef = matrix(0, nrow = args$p+1, ncol = SimuArg$initials$nrep),
                       Var = rep(0, SimuArg$initials$nrep),
                       Convergence = rep(FALSE, SimuArg$initials$nrep))
      next
    }
    
    if(kk %in% c("DRIV.cf.hz.ml.est.rateCal")){
      out[[kk]] = list(Coef = matrix(0, nrow = nrow(rate_propensity)*ncol(rate_hazard), ncol = SimuArg$initials$nrep),
                       Var = matrix(0, nrow = nrow(rate_propensity)*ncol(rate_hazard), ncol = SimuArg$initials$nrep),
                       Convergence = matrix(FALSE, nrow = nrow(rate_propensity)*ncol(rate_hazard), ncol = SimuArg$initials$nrep),
                       mse_propensity = matrix(0, nrow = nrow(rate_propensity)*ncol(rate_hazard), ncol = SimuArg$initials$nrep),
                       mse_surv = matrix(0, nrow = nrow(rate_propensity)*ncol(rate_hazard), ncol = SimuArg$initials$nrep),
                       rate_hazard = matrix(0, nrow = nrow(rate_propensity)*ncol(rate_hazard), ncol = SimuArg$initials$nrep),
                       rate_propensity = matrix(0, nrow = nrow(rate_propensity)*ncol(rate_hazard), ncol = SimuArg$initials$nrep))
      next
    }
    out[[kk]] = templist
  }
  for (i in sequence) {
    data = jsonlite::read_json(paste0(SimuArg$Control$save_path, 
                                      "DataGenerated/",
                                      SimuArg$initials$N,
                                      SimuArg$Control$Annotation, 
                                      "/", i, ".json"), simplifyVector = T)
    data$Covariates = data$Covariates[, 1:args$p, drop = FALSE]
    colnames(data$Covariates) = paste0("X", 1:args$p)
    args$dat = data
    cat("[[rep", i)
    for (kk in methods) {
      args$method = kk
      ModelPar = do.call(New_ModelPar, args)
      mod = DataFitting(ModelPar)
      out[[kk]]$Coef[, i] = mod$Coef
      cat("\t", out[[kk]]$Coef[, i], "\n")
      if (kk %in% c("DRIV", "DRIV.s", "DRIV.s2", "DRIV.cf.hz.est", "DRIV.s.control", "DRIV.cf.hz.ml.est")){
        if(is.null(mod$Var)) next
        out[[kk]]$Var[i] = mod$Var
        out[[kk]]$Convergence[i] = mod$Convergence
        cat("\t", out[[kk]]$Var[i])
        cat("\t", out[[kk]]$Convergence[i])
        next
      }
      if(kk %in% c("DRIV.cf.hz.ml.est.rateCal")) {
        out[[kk]]$Var[, i] = mod$Var
        out[[kk]]$Convergence[, i] = mod$Convergence
        out[[kk]]$mse_propensity[, i] = mod$mse_propensity
        out[[kk]]$mse_surv[, i] = mod$mse_surv
        # browser()
        out[[kk]]$rate_hazard[, i] = mod$rate_hazard
        out[[kk]]$rate_propensity[, i] = mod$rate_propensity
        cat("\t", out[[kk]]$Var[, i], "\n")
        cat("\t", out[[kk]]$Convergence[, i], "\n")
        next
      }
      out[[kk]]$Var[, i] = mod$Var
    }
    cat("]]\n")
  }
  SimuArg$methods = methods
  SimuArg$SimuResults = out
  structure(SimuArg, class = "SimuResults")
}


TRTSWE = function(dat, max_t, methods, 
                 Control = list(), ...)
{
  N = nrow(dat$Covariates)
  p = ncol(dat$Covariates)
  Validate_method(methods)
  args = rlang::dots_list(N = N, p = p, max_t = max_t,
                          dat = dat, Control = Control,
                          !!!list(...), .homonyms = "first")
  ModelPar = vector("list", length(methods)) 
  names(ModelPar) = methods
  
  for (kk in methods) {
    args$method = kk
    ModelPar[[kk]] = do.call(New_ModelPar, args)
    mod = DataFitting(ModelPar[[kk]])
    ModelPar[[kk]]$Estim = mod
  }
  structure(ModelPar, class = "TRTSWE")
}

