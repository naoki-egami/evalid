aipw_bart <- function(formula_outcome,
                      formula_weights,
                      treatment,
                      exp_data, pop_data,
                      weights_type = "calibration",
                      weights_max,
                      pop_weights = NULL,
                      boot = TRUE, sims = 1000, boot_ind = NULL,
                      numCores = 1, seed = 1234) {

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

  ## Determine outcome variable name and covariates
  ## (note, might instead track outcome var name)
  outcome_var <- all.vars(formula_outcome)[1]
  covariates <- all.vars(formula_outcome)[-1]

  ## Fit bart model
  bart_fit <- bartc(response = exp_data[, outcome_var],
                    treatment = exp_data[, treatment],
                    confounders = as.matrix(exp_data[, covariates]),
                    method.trt = "none", keepTrees = TRUE)

  # ###############################
  # Projection to the Population
  # ###############################
  std_pop_weights <- pop_weights / sum(pop_weights)
  ## project Y1
  bart_proj_1 <- t(predict(bart_fit,
                           newdata = pop_data[, covariates],
                           type = "y.1")) * std_pop_weights

  ## project Y0
  bart_proj_0 <- t(predict(bart_fit,
                           newdata = pop_data[, covariates],
                           type = "y.0")) * std_pop_weights

  # #######################
  # Predict on to sample
  # #######################
  pred_y1_samp <- t(predict(bart_fit, newdata = exp_data[, covariates, drop = FALSE], type = "y.1"))
  pred_y0_samp <- t(predict(bart_fit, newdata = exp_data[, covariates, drop = FALSE], type = "y.0"))

  ## Collect residuals
  pred <- pred_y1_samp*as.numeric(exp_data[, treatment] == 1) +
    pred_y0_samp*as.numeric(exp_data[, treatment] == 0)
  res <- as.vector(exp_data[, outcome_var]) - pred
  res_w1 <- apply(res * as.numeric(exp_data[, treatment] == 1)* ipw_weights, 2, sum)/sum(as.numeric(exp_data[, treatment] == 1)*ipw_weights)
  res_w0 <- apply(res * as.numeric(exp_data[, treatment] == 0)* ipw_weights, 2, sum)/sum(as.numeric(exp_data[, treatment] == 0)*ipw_weights)

  estimate <- apply((bart_proj_1 - bart_proj_0), 2, mean) + (res_w1 - res_w0)

  est <- mean(estimate, na.rm = TRUE)
  se <- sd(estimate, na.rm = TRUE)
  ci_lower <- quantile(estimate, 0.025, na.rm = TRUE)
  ci_upper <- quantile(estimate, 0.975, na.rm = TRUE)

  out <- list(est = est, se = se, ci_lower = ci_lower, ci_upper = ci_upper)
  return(out)
}
