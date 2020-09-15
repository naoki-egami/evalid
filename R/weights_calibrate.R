### Helper function for creating targets from auxiliary information and formula for

### calibration
create_targets <- function (target_design, target_formula) {
  target_mf <- model.frame(target_formula, model.frame(target_design))
  target_mm <- model.matrix(target_formula, target_mf)
  wts <- weights(target_design)
  colSums(target_mm * wts) / sum(wts)
}

## Construct calibration weights
weights_cal <- function(formula_weights, exp_data, pop_data,
                        calfun = "raking", weights_max = Inf,
                        pop_weights = NULL) {

  if(is.null(pop_weights)) pop_weights <- rep(1, nrow(pop_data))

  ## gather the covariate terms
  covariates <- labels(terms(formula_weights))

  ## Make calibration formula
  formula = as.formula(paste0("~ ", paste(covariates, collapse = " + ")))

  ## Make survey designs for calibration
  suppressWarnings(sample_surv <- svydesign(ids = ~ 1, data = exp_data))
  pop_surv <- svydesign(ids = ~ 1, data = pop_data, weights = pop_weights)

  ### Population targets (normalized to counts in sample size)
  targets <- create_targets(pop_surv, formula) * nrow(sample_surv)

  ### Weighted survey designs
  PATE_cal <- calibrate(design = sample_surv,
                        formula = formula,
                        population = targets,
                        calfun = calfun)

  ## Trim weights
  if(is.finite(weights_max)) {
    PATE_cal <- trimWeights(PATE_cal, upper = weights_max)
  }

  return(weights(PATE_cal))
}
