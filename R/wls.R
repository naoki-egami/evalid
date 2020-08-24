#' Weighted Least Squares
#' @param formula_outcome Formula for outcome model
#' @param treatment Name of the treatment variable
#' @param exp_data Experimental Data
#' @param weights Weights should be a vector of numeric values with the same length with `exp_data`
#' @export

wls <- function(formula_outcome, treatment, exp_data, weights = NULL, ...){

  ## runs lm robust with weights
  ## can pass in lm_robust params such as cluster
  res <- lm_robust(formula_outcome, data = exp_data, weights = weights, ...)

  ## selects treatment results
  est = res$coefficients[treatment]
  se = res$std.error[treatment]
  lm_w_ci = confint(res)[treatment, ]

  return(data.frame(type = "wls", est = est, se = se,
                    ci_lower = lm_w_ci[1], ci_upper = lm_w_ci[2],
                    n = res$N))
}
