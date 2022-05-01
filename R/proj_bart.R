## Outcome-based Estimator with BART
proj_bart <- function(formula_outcome,
                      treatment,
                      exp_data, pop_data,
                      pop_weights = NULL,
                      sims = 1000) {

  ## Determine outcome variable name and covariates
  ## (note, might instead track outcome var name)
  outcome_var <- all.vars(formula_outcome)[1]
  covariates <- all.vars(formula_outcome)[-1]

  ## Fit bart model
  bart_fit <- bartc(response = exp_data[, outcome_var],
                    treatment = exp_data[, treatment],
                    confounders = as.matrix(exp_data[, covariates]),
                    method.trt = "none", keepTrees = TRUE, n.samples = ceiling(sims/10))

  # ## breaks the projection down in to chunks to avoid overly large matrix
  # pop_data$chunks <- cut(1:nrow(pop_data), max(2, nrow(pop_data) %/% 1000))
  #

  ## project Y1
  bart_proj_1 <- t(t(predict(bart_fit,
                             newdata = pop_data[, covariates],
                             type = "y.1")) * pop_weights)

  ## project Y0
  bart_proj_0 <- t(t(predict(bart_fit,
                             newdata = pop_data[, covariates],
                             type = "y.0")) * pop_weights)

  ## distribution of PATE estimates
  bart_pates <- apply(bart_proj_1 - bart_proj_0, 1, mean, na.rm = TRUE)

  est <- mean(bart_pates, na.rm = TRUE)
  se <- sd(bart_pates, na.rm = TRUE)
  ci_lower <- quantile(bart_pates, 0.025, na.rm = TRUE)
  ci_upper <- quantile(bart_pates, 0.975, na.rm = TRUE)

  out <- list(est = est, se = se, ci_lower = ci_lower, ci_upper = ci_upper)
  return(out)
}
