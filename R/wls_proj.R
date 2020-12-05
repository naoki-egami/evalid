wls_proj <- function(formula_outcome,
                     formula_weights,
                     treatment,
                     exp_data, pop_data,
                     weights_type = "calibration",
                     weights_max,
                     pop_weights = NULL,
                     boot = TRUE, sims = 1000, boot_ind = NULL,
                     numCores = 1, seed = 1234,
                     ...){

  ## Naoki: I have not included pop_weights into the estimation

  # We only support boot = TRUE
  if(weights_type == "logit"){
    ipw_weights <- weights_logit(formula_weights = formula_weights,
                                 exp_data = exp_data, pop_data = pop_data,
                                 weights_max = weights_max, pop_weights = pop_weights)
  }else if(weights_type == "calibration"){
    ipw_weights <- weights_cal(formula_weights = formula_weights,
                               exp_data = exp_data, pop_data = pop_data,
                               calfun = "raking", weights_max = weights_max,
                               pop_weights = pop_weights,
                               ...)
  }

  out_m <- gen_bootstrap(est = wls_proj_base, numCores = numCores, sims = sims,
                         formula_outcome = formula_outcome, formula_weights = NULL,
                         treatment = treatment, exp_data = exp_data, pop_data = pop_data,
                         weights = ipw_weights, pop_weights = pop_weights,
                         boot_ind = boot_ind, seed = seed)

  out <- list(out_m = out_m, ipw_weights = ipw_weights)
  return(out)
}



wls_proj_base <- function(x,
                          formula_outcome, formula_weights = NULL,
                          treatment, exp_data, pop_data,
                          weights = NULL,
                          pop_weights = NULL, boot_ind = NULL, seed = 1234) {

  # Boostrap the data
  boot_use <- boot_data(x = x, seed = seed, data = exp_data, boot_ind = boot_ind)
  exp_data_b <- exp_data[boot_use, , drop = FALSE]
  weights_b <- weights[boot_use]

  w_1 <- as.numeric(weights_b[exp_data_b[, treatment] == 1])
  exp_1 <- exp_data_b[exp_data_b[, treatment] == 1, , drop = FALSE]
  w_0 <- as.numeric(weights_b[exp_data_b[, treatment] == 0])
  exp_0 <- exp_data_b[exp_data_b[, treatment] == 0, , drop = FALSE]

  ## Model Y1 with weights
  lm_mod_1 <- lm_robust(formula_outcome, data = exp_1, weights = w_1)

  ## Model Y0 with weights
  lm_mod_0 <- lm_robust(formula_outcome, data = exp_0, weights = w_0)

  ## Project Y1
  lm_proj_1 <- predict(lm_mod_1, newdata = pop_data)

  ## Project Y0
  lm_proj_0 <- predict(lm_mod_0, newdata = pop_data)

  out <- weighted.mean(lm_proj_1  - lm_proj_0, w = pop_weights)

  return(out)
}
