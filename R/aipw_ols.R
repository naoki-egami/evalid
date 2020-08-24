# Augumented IPW Estimator with OLS
aipw_ols <- function(formula_outcome, treatment,
                     formula_weights,
                     model, weights_type = "calibration",
                     exp_data,
                     pop_data, weight_max, ...){

  if(weights_type == "ipw"){
    weights <- weights_ipw(formula_weights = formula_weights,
                           exp_data = exp_data,
                           pop_data = pop_data, weight_max = weight_max, ...)
  }else if(weights_type == "calibration"){
    weights <- weights_cal(formula_weights = formula_weights,
                           exp_data = exp_data,
                           pop_data = pop_data, weight_max = weight_max, ...)
  }

  ## do ols projection, get models back
  proj <- proj_ols(formula_outcome = formula_outcome,
                   exp_data = exp_data, pop_data = pop_data,
                   treatment = treatment, return_models = TRUE, ...)

  ## get residuals
  exp_data$resid <- exp_data$Y -
    ifelse(exp_data$treatment == 1,
           ## if true project model Y1
           predict(proj$model[[1]], newdata = exp_data),
           ## if false, project model Y0
           predict(proj$model[[2]], newdata = exp_data))
  ## add wls on residuals to projection estimate
  res <- proj$pate_proj + wls(exp_data,
                              ## run model on residuals + treatment_var
                              update(resid ~ 1, paste("~ . +", treatment)),
                              treatment = treatment, weights = weights, ...)$est

  return(res)
}
