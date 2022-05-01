#' Weighted Least Squares

wls <- function(formula_outcome,
                treatment, exp_data,
                weights = NULL,
                pop_weights = NULL,
                boot = TRUE, sims = 1000, boot_ind = NULL,
                numCores = 1, seed = 1234){

  ## Naoki: I have not included pop_weights into the estimation

  if(boot == FALSE){
    ## runs lm robust with weights
    ## can pass in lm_robust params such as cluster
    res <- lm_robust(formula_outcome, data = exp_data, weights = weights)

    ## Standard Errors from lm_robust
    est = res$coefficients[treatment]
    se = res$std.error[treatment]
    lm_w_ci = confint(res)[treatment, ]

    # results
    out <- list(est = est, se = se, ci_lower = lm_w_ci[1], ci_upper = lm_w_ci[2])
  }else{

    out <- gen_bootstrap(est = wls_base, numCores = numCores, sims = sims,
                         formula_outcome = formula_outcome, formula_weights = NULL,
                         treatment = treatment, exp_data = exp_data, pop_data = NULL,
                         weights = weights, pop_weights = pop_weights,
                         boot_ind = boot_ind, seed = seed)
  }

  return(out)
}

#' Basis of Weighted Least Squares
#' @param formula_outcome Formula for outcome model
#' @param treatment Name of the treatment variable
#' @param exp_data Experimental Data
#' @param weights Weights should be a vector of numeric values with the same length with `exp_data`

wls_base <- function(x,
                     formula_outcome, formula_weights = NULL,
                     treatment, exp_data, pop_data = NULL,
                     weights = NULL,
                     pop_weights = NULL, boot_ind = NULL, seed = 1234){

  ## Naoki: I have not included pop_weights into the estimation

  # Boostrap the data
  boot_use <- boot_data(x = x, seed = seed, data = exp_data,
                        boot_ind = boot_ind)
  exp_data_b <- exp_data[boot_use, , drop = FALSE]
  weights_b <- weights[boot_use]

  ## runs lm robust with weights
  ## can pass in lm_robust params such as cluster
  res <- lm_robust(formula_outcome, data = exp_data_b, weights = weights_b)

  ## selects treatment results
  est = res$coefficients[treatment]
  out <- est
  return(out)
}
