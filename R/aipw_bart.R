aipw_bart <- function(formula_outcome, treatment,
                      formula_weights,
                      model, weights_type = "calibration",
                      exp_data, pop_data, bart_pates, 
                      pop_weights = NULL, ...) {

  outcome_var <- all.vars(formula_outcome)[attr(terms(formula_outcome), "response")]
  covariates <- labels(terms(formula_outcome))
  
  ## assumes that bart_pates projection already incorporates pop_weights
  if(is.null(pop_weights)) pop_weights <- rep(1, nrow(pop_data))

  if(weights_type == "ipw"){
    weights <- weights_ipw(formula_weights = formula_weights,
                           exp_data = exp_data,
                           pop_data = pop_data, weight_max = weight_max, 
                           pop_weights = pop_weights, ...)
  }else if(weights_type == "calibration"){
    weights <- weights_cal(formula_weights = formula_weights,
                           exp_data = exp_data,
                           pop_data = pop_data, weight_max = weight_max, 
                           pop_weights = pop_weights, ...)
  }

  ## Predict on to sample
  pred_y1_samp = predict(model, newdata = exp_data[, covariates], type = "y.1")
  pred_y0_samp = predict(model, newdata = exp_data[, covariates], type = "y.0")

  ## Collect residuals
  residuals <- lapply(1:nrow(pred_y1_samp), function(i) {
    resid0 <- as.vector(exp_data[, outcome_var])
    resid0 <- resid0 - ifelse(exp_data[treatment_var] == 1,
                              ## if true, Y1 prediction
                              pred_y1_samp[i, ],
                              ## if false, Y0 prediction
                              pred_y0_samp[i, ])
  })

  ## do wLS on residuals (for each draw of posterior)
  w_residuals <- unlist(lapply(residuals, function(resid) {
    est <- wls(formula_outcome = resid ~ treatment,
               exp_data = data.frame(resid = as.vector(resid),
                                     treatment = exp_data[, treatment]),

               treatment = "treatment",
               weights = weights)$est
    return(est)
  }))

  # calculate augmented values
  pate_aipw <- bart_pates + w_residuals
  pate_ci_aipw = quantile(pate_aipw, c(alpha/2, 1 - (alpha/2)))

  res <- data.frame(type = "bart_aipw",
                    est = mean(pate_aipw),
                    se = sd(pate_aipw),
                    ci_lower = pate_ci_aipw[1],
                    ci_upper = pate_ci_aipw[2],
                    n = nrow(exp_data))
  return(res)
}
