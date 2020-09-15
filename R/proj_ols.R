## OLS Projection
proj_ols <- function(formula_outcome, exp_data, pop_data,
                     treatment,
                     pop_weights = NULL,
                     boot = TRUE,
                     sims = 1000, boot_ind = NULL,
                     numCores = 1, seed = 1234) {

  if(boot == FALSE){
    ## Run models separately on Y1 and Y0
    lm_mod_1 <- lm(formula_outcome,
                   data = exp_data[exp_data[, treatment] == 1, , drop = FALSE])
    lm_mod_0 <- lm(formula_outcome,
                   data = exp_data[exp_data[, treatment] == 0, , drop = FALSE])

    coef_t <- mvrnorm(n = sims, mu = coef(lm_mod_1), Sigma = vcov(lm_mod_1))
    coef_c <- mvrnorm(n = sims, mu = coef(lm_mod_0), Sigma = vcov(lm_mod_0))

    fit_X <- paste0("~", as.character(formula(lm_mod_1))[3])

    X_pop <- model.matrix(formula(fit_X), data = pop_data)

    predict_t <- X_pop %*% t(coef_t)
    predict_c <- X_pop %*% t(coef_c)

    pate_sim <- apply(predict_t - predict_c, 2, weighted.mean, w = pop_weights)
    est <- mean(pate_sim, na.rm = TRUE)
    se <- sd(pate_sim, na.rm = TRUE)
    ci_lower <- quantile(pate_sim, prob = 0.025, na.rm = TRUE)
    ci_upper <- quantile(pate_sim, prob = 0.975, na.rm = TRUE)

    # results
    out <- list(est = est, se = se, ci_lower = ci_lower, ci_upper = ci_upper,
                lm_mod_1 = lm_mod_1, lm_mod_0 = lm_mod_0)

    return(out)
  }else{

    out <- gen_bootstrap(est = proj_ols_base, numCores = numCores, sims = sims,
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
