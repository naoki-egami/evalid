#' Estimating the Target Population Average Treatment Effect (TPATE)
#' @param outcome Mame of the outcome variable
#' @param treatment Name of the treatment variable
#' @param covariates Names of covariates we adjust for.
#' @param covariates_exp (Optional) Names of covariates measured only in the experimental data. These covariates are used only in `wls`
#' @param exp_data Experimental data. This should be `data.frame`
#' @param pop_data Population data. This should be `data.frame`
#' @param est_type Estimator Type. Should be one of the six types (`ipw`, `wls`, `outcome-ols`, `outcome-bart`, `dr-ols`, `dr-bart`). Weighting-based estimators are `ipw` and `wls`. Outcome-based estimators are `outcome-ols` and `outcome-bart.` Doubly-robust estimators are `dr-ols` and `dr-bart`.
#' @param weights_type Weights Type. Should be one of the two types (`calibration`, `logit`). We recommend `calibration` for its stability in many applications.
#' @param weights_max Maximum weights. Default is `Inf`.
#' @param boot Logical. Whether we use bootstrap (`TRUE` or `FALSE`). For doubly robust estimators (`dr_ols` and `dr-bart`), we always use bootstrap (`boot = TRUE`).
#' @param id_cluster Identifies for cluster bootstrap. The length should be the same as the number of rows of `exp_data`.
#' @param sims The number of bootstrap replications (when `boot = TRUE)`. The number of quasi-bayesian posterior samples (when `boot = FALSE`).
#' @param numCores The number of cores we use for parallel computing.
#' @param seed seed (default = `1234`)
#' @param compute_sate Logical. Whether we compute the SATE within the package.
#' @import estimatr
#' @import bartCause
#' @import survey
#' @import tidyverse
#' @importFrom pbmcapply pbmclapply
#' @importFrom pbapply pblapply
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach "%dopar%" "%do%" foreach
#' @importFrom MASS mvrnorm
#' @return \code{tpate} returns the following values.
#'  \itemize{
#'    \item \code{sate}: Estimates of the SATE
#'    \item \code{tpate}: Estimates of the T-PATE
#'    \item \code{ipw_weights}: Estimated weights (when weighting-based or doubly robust estimators are used)
#'  }
#' @description \code{tpate} implements the effect-generalization.
#' @references Egami and Hartman. (2022+). Elements of External Validity: Framework, Design, and Analysis.
#' @export


tpate <- function(outcome,
                  treatment,
                  covariates,
                  covariates_exp = NULL,
                  exp_data,
                  pop_data,
                  est_type = "wls",
                  weights_type = "calibration",
                  weights_max = Inf,
                  boot = TRUE, id_cluster = NULL,
                  sims = 500, numCores = NULL, seed = 1234,
                  compute_sate = FALSE) {

  # Housekeeping

  ## data.frame
  if(("data.frame" %in% class(exp_data)) == FALSE){
    stop("'exp_data' should be 'data.frame'")
  }else{
    class(exp_data) <- "data.frame"
  }
  if(("data.frame" %in% class(pop_data)) == FALSE){
    stop("'pop_data' should be 'data.frame'")
  }else{
    class(pop_data) <- "data.frame"
  }

  ## est_type
  if((est_type %in% c("ipw","wls", "outcome-ols", "outcome-bart", "dr-ols", "dr-bart")) == FALSE){
    stop("`est_type` should be one of `ipw`, `wls`, `outcome-ols`, `outcome-bart`, `dr-ols`, or `dr-bart`.")
  }

  # ###################
  # Setup
  # ###################
  formula_outcome <- as.formula(paste0(outcome, "~", paste(c(treatment, covariates), collapse = "+")))
  formula_weights <- as.formula(paste0("~", paste(c(covariates), collapse = "+")))

  if(is.null(numCores) == FALSE){
    if(numCores >= detectCores()) numCores <- detectCores() - 1
  }
  if(is.null(numCores)) numCores <- detectCores() - 1


  ## If pop_weights are null, assign as all 1
  pop_weights <- NULL
  if(is.null(pop_weights)) pop_weights <- rep(1, nrow(pop_data))

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
      formula_weighting <- as.formula(paste0(outcome, "~", paste(c(treatment, covariates, covariates_exp), collapse = "+")))
    }
    tpate_fit <- wls(formula_outcome = formula_weighting,
                     treatment = treatment,
                     exp_data = exp_data,
                     weights = ipw_weights,
                     pop_weights = pop_weights,
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
    ipw_weights <- NULL
  }

  # 3. Doubly-Robust Estimators
  if(is.null(pop_weights)) pop_weights <- rep(1, nrow(pop_data))
  if(est_type == "dr-ols" | est_type == "dr-bart"){
    if(est_type == "dr-ols") {
      formula_proj <- update(formula_outcome, paste("~ . -", treatment))
      tpate_fit_b <- aipw_ols(formula_outcome = formula_proj,
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
      tpate_fit_b <- aipw_bart(formula_outcome = formula_proj,
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
    tpate_fit <- tpate_fit_b$out_m
    ipw_weights <- tpate_fit_b$ipw_weights
  }

  out <- list(sate = sate_fit,
              tpate = tpate_fit,
              ipw_weights = ipw_weights)

  return(out)
}
