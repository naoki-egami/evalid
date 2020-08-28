## Calculate logit weights
weights_ipw <- function(formula_weights, exp_data, pop_data, weight_max = Inf,
                        pop_weights = NULL, ...) {

  ## gather the covariate terms
  covariates <- labels(terms(formula_weights))

  ## make stacked dataset, define sampling indicator "S"
  full_data <- rbind(data.frame(exp_data[, covariates], S = 1),
                     data.frame(pop_data[, covariates], S = 0))
  
  if(is.null(pop_weights)) {
    full_data$pop_weights_ipw = rep(1, nrow(full_data))
  } else {
    full_data$pop_weights_ipw = c(rep(1, nrow(sample)), pop_weights)
  }

  ## make formula of S on collapsed covariates
  formula_ipw <- as.formula(paste0("S ~ ", paste0(covariates, collapse = " + ")))

  ## get weights
  ipw_fit <- glm(formula_ipw, data = full_data, family = binomial(link = "logit"), 
                 weights = pop_weights_ipw)
  p_expt  <- predict(ipw_fit, type = "response")
  weights <- 1 / p_expt * (1 - p_expt) / mean(full_data$S == 0)
  weights_sample <- weights[full_data$S == 1]
  weights_sample <- weights_sample/mean(weights_sample)

  ## trim weights
  # if(is.finite(weight_max)) {
  #   weights_sample <- weights(trimWeights(svydesign(~1, data = sample,
  #                                                   weights = weights_sample),
  #                                         upper = weight_max))
  # }

  weights_sample[weights_sample > weight_max] <- weight_max


  return(weights_sample)
}
