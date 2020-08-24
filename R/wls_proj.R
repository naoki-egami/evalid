wls_projection <- function(formula_outcome, treatment,
                           exp_data, pop_data, weights, ...) {

  ## Model Y1 with weights
  lm_mod_1 <- lm_robust(formula_outcome, weights = weights[exp_data[, treatment_var] == 1],
                        data = exp_data[exp_data[, treatment_var] == 1, , drop = FALSE],
                        ...)

  ## Model Y0 with weights
  lm_mod_0 <- lm_robust(formula_outcome, weights = weights[exp_data[, treatment_var] == 0],
                        data = exp_data[exp_data[, treatment_var] == 0, , drop = FALSE],
                        ...)

  ## Project Y1
  lm_proj_1 <- predict(lm_mod_1, newdata = pop_data)

  ## Project Y0
  lm_proj_0 <- predict(lm_mod_0, newdata = pop_data)

  pate_proj <- mean(lm_proj_1  - lm_proj_0)

  return(pate_proj)
}
