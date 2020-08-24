#' Estimating the Target Population Average Treatment Effect (TPATE)
#' @param formula_outcome Formula for outcome model
#' @param formula_weights Formula for sampling weights
#' @param treatment Name of the treatment variable
#' @param exp_data Experimental Data
#' @param pop_data Population Data
#' @param estimator_type Whether we use a paired-choice conjoint design
#' @import estimatr
#' @import bartCause
#' @import survey
#' @import tidyverse
#' @return \code{model_pAMCE} returns an object of \code{pAMCE} class.
#'  \itemize{
#'    \item \code{AMCE}: Estimates of the pAMCE for all factors.
#'    \item \code{boot_AMCE}: Estimates of the pAMCE for all factors in each bootstrap sample.
#'    \item \code{boot_coef}: Estimates of coefficients for the linear probability model in each bootstrap sample.
#'    \item \code{approach}: "model_based"
#'    \item \code{input}: Input into the function.
#'    \item \code{...}: Values for internal use.
#'  }
#' @description \code{tpate} implements the effect-generalization.
#' @references Egami and Hartman. (2020+). Elements of External Validity: Framework, Design, and Analysis
#' @export


tpate <- function(formula_outcome,
                  formula_weights,
                  treatment,
                  exp_data,
                  pop_data,
                  type = "weighting",
                  outcome_type = "ols",
                  weights_type = "calibration",
                  ...) {

  ## Error handling that needs to be done:
  ## - are covariates present in both population and exp data (unless using wLS)

  ## Note:
  ## Need to take in outcome_var for weighted methods, right now assumes "Y"
  ## Need to take in alpha level, right now assumes 0.05 (for the bart posteriors)

  ## Naoki: I will work on Error Handling Later.

  ## ipw-logit
  if(estimator_type == "sate") {
    res <- wls(exp_data,
               ## Run formula on outcome + treatment, no weights
               formula_outcome = update(formula(Y ~ 1), paste("~ . +", treatment)),
               treatment = treatment, weights = NULL, ...)$est
  }

  ## ipw-logit
  if(estimator_type == "ipw-logit") {
    ## get logit weights
    weights <- ipw_weights(exp_data, pop_data, formula_weights)

    res <- wls(exp_data,
               ## Run formula on outcome + treatment, with weights
               formula_outcome = update(formula(Y ~ 1), paste("~ . +", treatment)),
               treatment = treatment, weights = weights, ...)$est
  }

  ## ipw-cal
  if(estimator_type == "ipw-cal") {
    ## get calibration weights
    weights <- cal_weights(exp_data, pop_data, formula_weights)
    res <- wls(exp_data,
               ## Run formula on outcome + treatment, with weights
               formula_outcome = update(formula(Y ~ 1), paste("~ . +", treatment)),
               treatment = treatment, weights = weights, ...)$est
  }

  ## wls-controls-logit
  if(estimator_type == "wls-controls-logit") {
    ## get logit weights
    weights <- ipw_weights(exp_data, pop_data, formula_weights)
    res <- wls(exp_data,
               ## Run formula on formula_outcome + treatment, with weights
               formula_outcome = update(formula_outcome, paste("~ . +", treatment)),
               treatment = treatment, weights = weights, ...)$est
  }

  ## wls-controls-cal
  if(estimator_type == "wls-controls-cal") {
    ## get calibration weights
    weights <- cal_weights(exp_data, pop_data, formula_weights)
    res <- wls(exp_data,
               ## Run formula on formula_outcome + treatment, with weights
               formula_outcome = update(formula_outcome, paste("~ . +", treatment)),
               treatment = treatment, weights = weights, ...)$est
  }

  ## outcome-ols
  if(estimator_type == "outcome-ols") {
    res <- ols_proj(exp_data, pop_data,
                    formula_outcome = formula_outcome)
  }

  ## outcome-bart
  ## would like to be able to pass in the following for bart
  ## nsample = 1000, nburn = 1000, nchains = 4
  if(estimator_type == "outcome-bart") {
    res <- bart_projection(exp_data, pop_data,
                           formula_outcome = formula_outcome,
                           treatment = treatment, ...)
  }

  ## "dr-bart-logit"
  if(estimator_type == "dr-bart-logit") {
    ## get logit weights
    weights <- ipw_weights(exp_data, pop_data, formula_weights = formula_weights)
    ## do bart projection and bart aipw
    res <- bart_projection(exp_data, pop_data,
                           formula_outcome = formula_outcome,
                           treatment = treatment,
                           do_aipw = TRUE, weights = weights, ...)
  }

  ## "dr-bart-cal"
  if(estimator_type == "dr-bart-cal") {
    ## get calibration weights
    weights <- cal_weights(exp_data, pop_data, formula_weights = formula_weights)
    ## do bart projection and bart aipw
    res <- bart_projection(exp_data, pop_data,
                           formula_outcome = formula_outcome,
                           treatment = treatment,
                           do_aipw = TRUE, weights = weights, ...)
  }

  ## dr-ols-logit
  if(estimator_type == "dr-ols-logit") {
    ## get logit weights
    weights <- ipw_weights(exp_data, pop_data, formula_weights = formula_weights)
    ## do ols projection, get models back
    proj <- ols_proj(exp_data, pop_data,
                     formula_outcome = formula_outcome,
                     return_models = TRUE, ...)

    ## get residuals
    exp_data$resids <- exp_data$Y -
      ifelse(exp_data$treatment == 1,
             ## if true project model Y1
             predict(proj$mod_1, newdata = exp_data),
             ## if false, project model Y0
             predict(proj$mod_0, newdata = exp_data))
    ## add wls on residuals to projection estimate
    res <- proj$pate_proj +
      wls(exp_data,
          ## run model on residuals + treatment
          update(resids ~ 1, paste("~ . +", treatment)),
          treatment = treatment, weights = weights, ...)$est
  }

  ## dr-ols-cal
  if(estimator_type == "dr-ols-cal") {
    ## get calibration weights
    weights <- cal_weights(exp_data, pop_data, formula_weights = formula_weights)
    ## do ols projection, get models back
    proj <- ols_proj(exp_data, pop_data, formula_outcome = formula_outcome, return_models = TRUE, ...)
    ## get residuals
    exp_data$res <- exp_data$Y -
      ifelse(exp_data$treatment == 1,
             ## if true project model Y1
             predict(proj$mod_1, newdata = exp_data),
             ## if false, project model Y0
             predict(proj$mod_0, newdata = exp_data))
    ## add wls on residuals to projection estimate
    res <- proj$pate_proj +
      wls(exp_data,
          ## run model on residuals + treatment
          update(res ~ 1, paste("~ . +", treatment)),
          treatment = treatment, weights = weights, ...)$est
  }

  ## dr-wls-logit
  if(estimator_type == "dr-wls-logit") {
    ## get logit weights
    weights <- ipw_weights(exp_data, pop_data, formula_weights = formula_weights)
    ## get doubly robust wLS projection
    res <- wls_projection(exp_data, pop_data,
                          formula_outcome = formula_outcome,
                          treatment = treatment,
                          weights = weights, ...)
  }

  ## dr-wls-cal
  if(estimator_type == "dr-wls-cal") {
    ## get calibratoin weights
    weights <- cal_weights(exp_data, pop_data, formula_weights = formula_weights)
    ## get doubly robust wLS projection
    res <- wls_projection(exp_data, pop_data,
                          formula_outcome = formula_outcome,
                          treatment = treatment,
                          weights = weights, ...)
  }

  return(res)
}
