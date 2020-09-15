#' Generic Bootstrap function
#' @param est Function of Estimators
#' @param formula_outcome Formula for outcome model
#' @param formula_weights Formula for sampling weights
#' @param treatment Name of the treatment variable
#' @param exp_data Experimental Data
#' @param pop_data Population Data
#' @export


gen_bootstrap <- function(est, numCores, sims = 1000,
                          formula_outcome,
                          formula_weights,
                          treatment,
                          exp_data,
                          pop_data,
                          weights = NULL,
                          pop_weights = NULL,
                          boot_ind = NULL,
                          seed = 1234){

  message(paste("\nBootstrap (", sims, "):\n", sep = ""))


  if(Sys.info()[['sysname']] == 'Windows') {

    # ------------
    # Windows
    # ------------

    if (numCores == 1){
      fit_boot <- pblapply(1:sims, function(x)
        est(x,
            formula_outcome = formula_outcome,
            formula_weights = formula_weights,
            treatment = treatment,
            exp_data = exp_data,
            pop_data = pop_data,
            weights = weights,
            pop_weights = pop_weights,
            boot_ind = boot_ind,
            seed = seed))
    }else {

      cl <- makeCluster(numCores)
      registerDoParallel(cl)
      # pb <- txtProgressBar(max = boot, style = 3)
      # progress <- function(n) setTxtProgressBar(pb, n)
      # opts <- list(progress = progress)
      #
      #     cl <- makeCluster(numCores)
      #     registerDoParallel(cl)
      #     on.exit(stopCluster(cl))

      fit_boot <- foreach(i = 1:sims,
                          .packages = c("estimatr")) %dopar% {
                            est(x = i,
                                formula_outcome = formula_outcome,
                                formula_weights = formula_weights,
                                treatment = treatment,
                                exp_data = exp_data,
                                pop_data = pop_data,
                                weights = weights,
                                pop_weights = pop_weights,
                                boot_ind = boot_ind,
                                seed = seed)
                          }
      #close(pb)
      stopCluster(cl)
    }
  }else{

    # ------------
    # Mac
    # ------------
    fit_boot <- pbmclapply(seq(1:sims), function(x) est(x,
                                                        formula_outcome = formula_outcome,
                                                        formula_weights = formula_weights,
                                                        treatment = treatment,
                                                        exp_data = exp_data,
                                                        pop_data = pop_data,
                                                        weights = weights,
                                                        pop_weights = pop_weights,
                                                        boot_ind = boot_ind,
                                                        seed = seed),
                           mc.cores = numCores)
  }

  # combine results
  estimate <- unlist(fit_boot)

  est <- mean(estimate, na.rm = TRUE)
  se <- sd(estimate, na.rm = TRUE)
  ci_lower <- quantile(estimate, 0.025, na.rm = TRUE)
  ci_upper <- quantile(estimate, 0.975, na.rm = TRUE)

  out <- list(est = est, se = se, ci_lower = ci_lower, ci_upper = ci_upper)
  return(out)
}


#' Generate Bootstrap data
#' @param x
#' @param seed
#' @param data
#' @param boot_ind
#' @export

boot_data <- function(x, seed, data, boot_ind){

  # boot_ind should be the same length as data

  # set seed
  seed.b <- 1000*x + seed
  set.seed(seed.b)

  # Stratified Bootstrap (bootstrap within strata)
  uniq_strata <- unique(boot_ind)
  pos_strata <- lapply(uniq_strata, function(x) which(x == boot_ind))
  boot_which <- lapply(pos_strata, function(x) sample(x, size = length(x), replace = TRUE))
  # data_b <- data[unlist(boot_which), , drop = FALSE]
  boot_use <- unlist(boot_which)
  # Block Bootstrap (bootstrap blocks)
  # boot_id0 <- sample(unique(boot_ind), size = length(unique(boot_ind)), replace = TRUE)
  # # create bootstap sample with sapply
  # boot_which <- sapply(boot_id0, function(x) which(boot_ind == x))
  # new_boot_id <- rep(seq(1:length(boot_id0)), times = unlist(lapply(boot_which, length)))
  # data_b <- data[unlist(boot_which), , drop = FALSE]
  return(boot_use)
}

