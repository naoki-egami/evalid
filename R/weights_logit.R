## Calculate logit weights
weights_logit <- function(formula_weights, exp_data, pop_data, weights_max = Inf,
                          pop_weights = NULL, ...) {

  ## gather the covariate terms
  covariates <- labels(terms(formula_weights))

  ## make stacked dataset, define sampling indicator "S"
  full_data <- rbind(data.frame(exp_data[, covariates], S = 1),
                     data.frame(pop_data[, covariates], S = 0))

  if(is.null(pop_weights)) {
    full_data$pop_weights_logit = rep(1, nrow(full_data))
  } else {
    full_data$pop_weights_logit = c(rep(1, nrow(exp_data)), pop_weights)
  }

  ## make formula of S on collapsed covariates
  formula_logit <- as.formula(paste0("S ~ ", paste0(covariates, collapse = " + ")))

  ## get weights
  logit_fit <- glm(formula_logit, data = full_data, family = binomial(link = "logit"),
                   weights = pop_weights_logit)
  p_expt  <- predict(logit_fit, type = "response")
  weights <- 1 / p_expt * (1 - p_expt) / mean(full_data$S == 0)
  weights_sample <- weights[full_data$S == 1]
  weights_sample <- weights_sample/mean(weights_sample)

  ## trim weights
  # if(is.finite(weight_max)) {
  #   weights_sample <- weights(trimWeights(svydesign(~1, data = sample,
  #                                                   weights = weights_sample),
  #                                         upper = weight_max))
  # }

  weights_sample[weights_sample > weights_max] <- weights_max


  return(weights_sample)
}
