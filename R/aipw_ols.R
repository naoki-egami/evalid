# Augumented IPW Estimator with OLS
aipw_ols <- function(formula_outcome,
                     formula_weights,
                     treatment,
                     exp_data, pop_data,
                     weights_type = "calibration",
                     weights_max,
                     pop_weights = NULL,
                     boot = TRUE, sims = 1000, boot_ind = NULL,
                     numCores = 1, seed = 1234){

  # We only support boot = TRUE
  if(weights_type == "logit"){
    ipw_weights <- weights_logit(formula_weights = formula_weights,
                                 exp_data = exp_data, pop_data = pop_data,
                                 weights_max = weights_max, pop_weights = pop_weights)
  }else if(weights_type == "calibration"){
    ipw_weights <- weights_cal(formula_weights = formula_weights,
                               exp_data = exp_data, pop_data = pop_data,
                               calfun = "raking", weights_max = weights_max,
                               pop_weights = pop_weights)
  }

  ## Prepare
  formula_proj <- update(formula_outcome, paste("~ . -", treatment))

  out <- gen_bootstrap(est = aipw_ols_base, numCores = numCores, sims = sims,
                       formula_outcome = formula_proj,
                       formula_weights = NULL,
                       treatment = treatment,
                       exp_data = exp_data,
                       pop_data = pop_data,
                       weights = ipw_weights,
                       pop_weights = pop_weights,
                       boot_ind = boot_ind,
                       seed = seed)

  return(out)
}


# Augumented IPW Estimator with OLS
aipw_ols_base <- function(x,
                          formula_outcome,
                          formula_weights,
                          treatment,
                          exp_data,
                          pop_data,
                          weights,
                          pop_weights,
                          boot_ind,
                          seed){

  # Boostrap the data
  boot_use <- boot_data(x = x, seed = seed, data = exp_data, boot_ind = boot_ind)
  exp_data_b <- exp_data[boot_use, , drop = FALSE]
  weights_b <- weights[boot_use]

  ## do ols projection, get models back
  proj <- proj_ols(formula_outcome = formula_outcome,
                   treatment = treatment,
                   exp_data = exp_data_b,
                   pop_data = pop_data,
                   pop_weights = pop_weights,
                   boot = FALSE,
                   sims = 2, boot_ind = NULL,
                   numCores = 1, seed = 1234)

  Y <- exp_data_b[, all.vars(formula_outcome)[1]]
  treat <- exp_data_b[, treatment]

  ## get residuals
  exp_data_b$resid <- Y -
    ifelse(treat == 1,
           ## if true project model Y1
           predict(proj$lm_mod_1, newdata = exp_data_b),
           ## if false, project model Y0
           predict(proj$lm_mod_0, newdata = exp_data_b))
  ## add wls on residuals to projection estimate
  formula_res <- update(resid ~ 1, paste("~ . +", treatment))
  res_fit <- lm_robust(formula_res, data = exp_data_b, weights = weights_b)
  out <- proj$est + coef(res_fit)[treatment]

  return(out)
}




# # Augumented IPW Estimator with OLS
# aipw_ols_base <- function(formula_outcome, treatment,
#                           formula_weights,
#                           model, weights_type = "calibration",
#                           exp_data,
#                           pop_data, weight_max,
#                           pop_weights = NULL, ...){
#
#
#
#   if(weights_type == "logit"){
#     ipw_weights <- weights_logit(formula_weights = formula_weights,
#                                  exp_data = exp_data, pop_data = pop_data,
#                                  weight_max = weight_max, pop_weights = pop_weights)
#   }else if(weights_type == "calibration"){
#     ipw_weights <- weights_cal(formula_weights = formula_weights,
#                                exp_data = exp_data, pop_data = pop_data,
#                                calfun = "raking", weight_max = weight_max,
#                                pop_weights = pop_weights)
#   }
#
#   ## do ols projection, get models back
#   formula_proj <- update(formula_outcome, paste("~ . -", treatment))
#   proj <- proj_ols(formula_outcome = formula_proj,
#                    treatment = treatment,
#                    exp_data = exp_data,
#                    pop_data = pop_data,
#                    pop_weights = pop_weights,
#                    boot = FALSE,
#                    sims = 2, boot_ind = NULL,
#                    numCores = 1, seed = 1234)
#
#   Y <- exp_data[, all.vars(formula_outcome)[1]]
#   treat <- exp_data[, treatment]
#
#   ## get residuals
#   exp_data$resid <- Y -
#     ifelse(treat == 1,
#            ## if true project model Y1
#            predict(proj$lm_mod_1, newdata = exp_data),
#            ## if false, project model Y0
#            predict(proj$lm_mod_0, newdata = exp_data))
#   ## add wls on residuals to projection estimate
#   formula_res <- update(resid ~ 1, paste("~ . +", treatment))
#   res_fit <- lm_robust(formula_res, data = exp_data, weights = ipw_weights)
#   out <- proj$est + coef(res_fit)[treatment]
#
#   return(out)
# }

