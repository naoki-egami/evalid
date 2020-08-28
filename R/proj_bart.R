proj_bart <- function(formula_outcome, treatment,
                     exp_data, pop_data, alpha = 0.05, 
                     pop_weights = NULL,
                     return_models = FALSE, ...) {

  ## Determine outcome variable name and covariates
  ## (note, might instead track outcome var name)
  outcome_var <- all.vars(formula_outcome)[attr(terms(formula_outcome), "response")]
  covariates <- labels(terms(formula_outcome))
  
  if(is.null(pop_weights)) pop_weights <- rep(1, nrow(pop_data))

  ## Fit bart model
  bart_fit <- bartc(response = exp_data[, outcome_var],
                    treatment = exp_data[, treatment_var],
                    confounders = as.matrix(exp_data[, covariates]),
                    method.trt = "none", keepTrees = TRUE,
                    ...)

  # ## breaks the projection down in to chunks to avoid overly large matrix
  # pop_data$chunks <- cut(1:nrow(pop_data), max(2, nrow(pop_data) %/% 1000))
  #

  ## project Y1
  bart_proj_1 <- rowSums(predict(bart_fit,
                                 newdata = pop_data[, covariates],
                                 type = "y.1")) %*% pop_weights / sum(pop_weights)

  ## project Y0
  bart_proj_0 <- rowSums(predict(bart_fit,
                                 newdata = pop_data[, covariates],
                                 type = "y.0")) %*% pop_weights / sum(pop_weights)

  ## distribution of PATE estimates
  bart_pates <- bart_proj_1 - bart_proj_0

  pate <- mean(bart_pates)
  pate_se <- sd(bart_pates)
  pate_ci <- quantile(bart_pates, c(alpha/2, 1 - (alpha/2)))

  if(return_models == FALSE){
    res <- data.frame(type = "proj_bart", est = pate, se = pate_se,
                      ci_lower = pate_ci[1], ci_upper = pate_ci[2],
                      n = nrow(exp_data))
  }else{
    res <- data.frame(type = "proj_bart", est = pate, se = pate_se,
                      ci_lower = pate_ci[1], ci_upper = pate_ci[2],
                      n = nrow(exp_data),
                      model = bart_fit, bart_pates = bart_pates)
  }
  return(res)
}
