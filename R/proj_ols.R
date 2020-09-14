## OLS Projection
proj_ols <- function(formula_outcome, exp_data, pop_data,
                     treatment,
                     pop_weights = NULL,
                     boot = FALSE, boot_ind = NULL,
                     numCores = 1, seed = 1234) {

  if(boot == FALSE){
    # I will add model-based inference later.
  }else{

    out <- gen_bootstrap(est = proj_ols_base, numCores = numCores, boot = boot,
                         formula_outcome = formula_outcome, formula_weights = NULL,
                         treatment = treatment, exp_data = exp_data, pop_data = pop_data,
                         weights = NULL, pop_weights = pop_weights,
                         boot_ind = boot_ind, seed = seed)
  }

  return(out)
}

proj_ols_base <- function(x,
                          formula_outcome,
                          formula_weights = NULL,
                          treatment,
                          exp_data,
                          pop_data,
                          weights = NULL,
                          pop_weights,
                          boot_ind,
                          seed) {


  # Boostrap the data
  boot_use <- boot_data(x = x, seed = seed, data = exp_data, boot_ind = boot_ind)
  exp_data_b <- exp_data[boot_use, , drop = FALSE]

  ## Run models separately on Y1 and Y0
  lm_mod_1 <- lm(formula_outcome,
                 data = exp_data_b[exp_data_b[, treatment] == 1, , drop = FALSE])
  lm_mod_0 <- lm(formula_outcome,
                 data = exp_data_b[exp_data_b[, treatment] == 0, , drop = FALSE])

  ## Project on to population
  lm_proj_1 <- predict(lm_mod_1, newdata = pop_data)
  lm_proj_0 <- predict(lm_mod_0, newdata = pop_data)

  ## Mean difference in projected Y1 and Y0
  out <- weighted.mean(lm_proj_1  - lm_proj_0, w = pop_weights)

  ## Return models for aipw if requested
  return(out)
}
