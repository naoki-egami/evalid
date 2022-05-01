# Generic Bootstrap function
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
        tryCatch(est(x,
            formula_outcome = formula_outcome,
            formula_weights = formula_weights,
            treatment = treatment,
            exp_data = exp_data,
            pop_data = pop_data,
            weights = weights,
            pop_weights = pop_weights,
            boot_ind = boot_ind,
            seed = seed), error = function(e) NA))
    } else {

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
                            tryCatch(est(x = i,
                                formula_outcome = formula_outcome,
                                formula_weights = formula_weights,
                                treatment = treatment,
                                exp_data = exp_data,
                                pop_data = pop_data,
                                weights = weights,
                                pop_weights = pop_weights,
                                boot_ind = boot_ind,
                                seed = seed), error = function(e) NA)
                          }
      #close(pb)
      stopCluster(cl)
    }
  }else{

    # ------------
    # Mac
    # ------------
    fit_boot <- pbmclapply(seq(1:sims), function(x) {
      tryCatch(est(x,
                   formula_outcome = formula_outcome,
                   formula_weights = formula_weights,
                   treatment = treatment,
                   exp_data = exp_data,
                   pop_data = pop_data,
                   weights = weights,
                   pop_weights = pop_weights,
                   boot_ind = boot_ind,
                   seed = seed), error = function(e) NA)
      }, mc.cores = numCores)
  }

  # combine results
  estimate <- unlist(fit_boot)

  est <- mean(estimate, na.rm = TRUE)
  se <- sd(estimate, na.rm = TRUE)
  ci_lower <- quantile(estimate, 0.025, na.rm = TRUE)
  ci_upper <- quantile(estimate, 0.975, na.rm = TRUE)

  # count the number of NAs
  num_NA <- sum(is.na(estimate))

  # out <- list(est = est, se = se, ci_lower = ci_lower, ci_upper = ci_upper, num_NA = num_NA)
  out <- list(est = est, se = se, ci_lower = ci_lower, ci_upper = ci_upper)
  return(out)
}

boot_data <- function(x, seed, data, boot_ind){

  # boot_ind should be the same length as data

  # set seed
  seed.b <- 1000*x + seed
  set.seed(seed.b)

  if(is.null(boot_ind) == FALSE){
    # Block Bootstrap (within each treatment)
    boot_id0 <- sample(unique(boot_ind), size = length(unique(boot_ind)), replace = TRUE)
    # # create bootstap sample with sapply
    boot_which <- sapply(boot_id0, function(x) which(boot_ind == x))
    boot_use <- unlist(boot_which)
  }else if(is.null(boot_ind) == TRUE){
    boot_use <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
  }

  return(boot_use)
}

