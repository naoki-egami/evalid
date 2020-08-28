## OLS Projection
proj_ols <- function(formula_outcome, exp_data, pop_data,
                     treatment = "treatment", 
                     pop_weights = NULL, 
                     return_models = FALSE, ...) {
  
  if(is.null(pop_weights)) pop_weights <- rep(1, nrow(pop_data))

  ## Run models separately on Y1 and Y0
  lm_mod_1 <- lm(formula_outcome,
                 data = exp_data[exp_data[, treatment_var] == 1, , drop = FALSE])
  lm_mod_0 <- lm(formula_outcome,
                 data = exp_data[exp_data[, treatment_var] == 0, , drop = FALSE])

  ## Project on to population
  lm_proj_1 <- predict(lm_mod_1, newdata = data.frame(pop_data))
  lm_proj_0 <- predict(lm_mod_0, newdata = data.frame(pop_data))

  ## Mean difference in projected Y1 and Y0
  pate_proj <- weighted.mean(lm_proj_1  - lm_proj_0, w = pop_weights)

  ## Return models for aipw if requested
  if(!return_models) {
    return(pate_proj)
  } else {
    return(list(pate_proj = pate_proj,
                model = list("lm_mod_1" = lm_mod_1,
                             "lm_mod_0" = lm_mod_0)))
  }
}
