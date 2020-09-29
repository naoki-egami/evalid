#' Estimating the Target Population Average Treatment Effect (TPATE)
#' @param formula_outcome Formula for outcome model (should include treatment)
#' @param formula_weights Formula for sampling weights
#' @param treatment Name of the treatment variable
#' @param exp_data Experimental Data
#' @param pop_data Population Data
#' @param est_type Estimator Type; `ipw`, `outcome-ols`, `outcome-bart`, `dr-ols`, `dr-bart`, or `wls-proj`
#' @param weights_type Weights Type; `logit`, `calibration`
#' @param weights_max Default is `Inf`
#' @param boot Whether you do bootstrap (`TRUE` or `FALSE`)
#' @param id_cluster Identifies for cluster bootstrap
#' @param sims The number of simulations
#' @param numCores The number of cores we use
#' @param seed seed (default = `1234`)
#' @param compute_sate whether we compute SATE
#' @import estimatr
#' @import bartCause
#' @import survey
#' @import tidyverse
#' @importFrom pbmcapply pbmclapply
#' @importFrom pbapply pblapply
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach "%dopar%" "%do%" foreach
#' @importFrom MASS mvrnorm
#' @return \code{tpate} returns the following values.
#'  \itemize{
#'    \item \code{sate}: Estimates of the SATE
#'    \item \code{tpate}: Estimates of the T-PATE
#'  }
#' @description \code{tpate} implements the effect-generalization.
#' @references Egami and Hartman. (2020+). Elements of External Validity: Framework, Design, and Analysis
#' @export


tpate <- function(formula_outcome,
                  formula_weights,
                  treatment,
                  exp_data,
                  pop_data,
                  est_type = "wls",
                  weights_type = "calibration", weights_max = Inf,
                  boot = TRUE, id_cluster = NULL,
                  sims = 1000, numCores = NULL, seed = 1234,
                  compute_sate = TRUE,
                  ...) {

  ## Error handling that needs to be done:
  ## - are covariates present in both population and exp data (unless using wLS)

  # Check 1: Don't allow variable called, "weights"
  # Check 2: Require researchers to transform variables beforehand

  ## Note:
  ## Need to take in outcome_var for weighted methods, right now assumes "Y"
  ## Need to take in alpha level, right now assumes 0.05 (for the bart posteriors)

  ## Naoki: I will work on Error Handling Later.
  ## Naoki: I think bartCause is risky to use. We will change it later.

  # ###################
  # Setup
  # ###################
  if(is.null(numCores) == FALSE){
    if(numCores >= detectCores()) numCores <- detectCores() - 1
  }
  if(is.null(numCores)) numCores <- detectCores() - 1


  ## If pop_weights are null, assign as all 1
  pop_weights <- NULL
  if(is.null(pop_weights)) pop_weights <- rep(1, nrow(pop_data))

  # Keep names of outcome variable
  outcome <- all.vars(formula_outcome)[1]

  # prepare for Bootstrap
  # For now, I assume that we only do the complete randomization.
  boot_ind <- id_cluster

  if(compute_sate == TRUE) {
    ## Estimate SATE with difference in means
    formula_sate <- as.formula(paste0(outcome, "~", paste0(treatment, collapse = " + ")))
    ## Run formula on outcome + treatment, no weights
    sate_fit <- wls(formula_outcome = formula_sate,
                    treatment = treatment,
                    exp_data = exp_data,
                    weights = NULL, pop_weights = pop_weights,
                    boot = boot,
                    sims = sims, boot_ind = boot_ind,
                    numCores = numCores, seed = seed)
  } else{
    sate_fit <- NA
  }


  # 1. Weighting-based Estimators
  if(est_type == "ipw" | est_type == "wls"){
    if(weights_type == "logit") {
      ## get logit weights
      ipw_weights <- weights_logit(formula_weights = formula_weights,
                                   exp_data = exp_data, pop_data = pop_data,
                                   weights_max = weights_max, pop_weights = pop_weights)
    }else if(weights_type == "calibration"){
      ## get calibration weights
      ipw_weights <- weights_cal(formula_weights = formula_weights,
                                 exp_data = exp_data, pop_data = pop_data,
                                 calfun = "raking", weights_max = weights_max,
                                 pop_weights = pop_weights)
    }

    if(est_type == "ipw"){
      # 1.1 IPW (est_type = "ipw")
      formula_weighting <- as.formula(paste0(outcome, "~", paste0(treatment, collapse = " + ")))
    }else if(est_type == "wls"){
      # 1.2 wls (est_type = "wls")
      # We are not currently accepting interactions between X and treatments
      formula_weighting <- update(formula_outcome, paste("~ . +", treatment))
    }
    tpate_fit <- wls(formula_outcome = formula_weighting,
                     treatment = treatment,
                     exp_data = exp_data,
                     weights = ipw_weights, pop_weights = pop_weights,
                     boot = boot,
                     sims = sims, boot_ind = boot_ind,
                     numCores = numCores, seed = seed)
  }

  # 2. Outcome-based Estimators
  if(est_type == "outcome-ols" | est_type == "outcome-bart"){
    formula_proj <- update(formula_outcome, paste("~ . -", treatment))
    if(is.null(pop_weights)) pop_weights <- rep(1, nrow(pop_data))


    if(est_type == "outcome-ols") {
      # 2.1: OLS (outcome-ols)
      tpate_fit <- proj_ols(formula_outcome = formula_proj,
                            treatment = treatment,
                            exp_data = exp_data,
                            pop_data = pop_data,
                            pop_weights = pop_weights,
                            boot = boot,
                            sims = sims, boot_ind = boot_ind,
                            numCores = numCores, seed = seed)
    }else if(est_type == "outcome-bart"){
      # 2.2: BART (outcome-bart)
      ## would like to be able to pass in the following for bart
      ## nsample = 1000, nburn = 1000, nchains = 4
      tpate_fit <- proj_bart(formula_outcome = formula_proj,
                             treatment = treatment,
                             exp_data = exp_data,
                             pop_data = pop_data,
                             pop_weights = pop_weights,
                             sims = sims)
    }
  }

  # 3. Doubly-Robust Estimators
  if(is.null(pop_weights)) pop_weights <- rep(1, nrow(pop_data))
  if(est_type == "dr-ols") {
    formula_proj <- update(formula_outcome, paste("~ . -", treatment))
    tpate_fit <- aipw_ols(formula_outcome = formula_proj,
                          formula_weights = formula_weights,
                          treatment = treatment,
                          exp_data = exp_data,
                          pop_data = pop_data,
                          weights_type = weights_type,
                          weights_max = weights_max,
                          pop_weights = pop_weights,
                          boot = TRUE, sims = sims, boot_ind = boot_ind,
                          numCores = numCores, seed = seed)
  }else if(est_type == "dr-bart") {
    formula_proj <- update(formula_outcome, paste("~ . -", treatment))
    tpate_fit <- aipw_bart(formula_outcome = formula_proj,
                           formula_weights = formula_weights,
                           treatment = treatment,
                           exp_data = exp_data,
                           pop_data = pop_data,
                           weights_type = weights_type,
                           weights_max = weights_max,
                           pop_weights = pop_weights,
                           boot = TRUE, sims = sims, boot_ind = boot_ind,
                           numCores = numCores, seed = seed)
  }else if(est_type == "wls-proj") {
    formula_proj <- update(formula_outcome, paste("~ . -", treatment))
    tpate_fit <- wls_proj(formula_outcome = formula_proj,
                          formula_weights = formula_weights,
                          treatment = treatment,
                          exp_data = exp_data,
                          pop_data = pop_data,
                          weights_type = weights_type,
                          weights_max = weights_max,
                          pop_weights = pop_weights,
                          boot = TRUE, sims = sims, boot_ind = boot_ind,
                          numCores = numCores, seed = seed)
  }

  out <- list(sate = sate_fit, tpate = tpate_fit)

  return(out)
}
