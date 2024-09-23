Validate_scenario = function(x){
  values = unclass(x)
  stopifnot(all(values %in% c("exogenous", "endogenous")))
  stopifnot(length(x) <= 1)
}


Validate_method = function(x){
  values = unclass(x)
  stopifnot(all(values %in% c("ITT", "remove", "recensor", "TimeVar", "DRIV", 
                              "DRIV.s", "DRIV.s.nleqslv", "DRIV.cf.hz", "DRIV.cf.hz.est", "DRIV.cf.chz","cIV",
                              "DRIV.s.control", "remove.control", "counterfactual",
                              "customized", "DRIV.cf.hz.ml.est", "DRIV.cf.hz.ml.est.rateCal",
                              "DRIV.cf.hz.ml.est.ave",
                              "DRIV.cf.hz.est.ave")))
  # stopifnot(length(x) <= 1)
}
