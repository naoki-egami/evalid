# simulate
library(MASS)
library(boot)

dgp_exp_pop <- function(total_size,
                        K, K_W, coef_tau, coef_select,
                        seed_use = 1234, sim_index  = 1, sd_Y = 0.1){

  set.seed(seed_use)
  coef_0 <- rnorm(n = K+1)
  # coef_tau <- rnorm(n = K_W)*tau_beta
  # coef_select <- rnorm(n = K_W)*select_beta

  set.seed(sim_index)
  # Super-Population
  # Sigma_X0  <- rnorm(n = K)
  # Sigma_X <-  (Sigma_X0 %*% t(Sigma_X0))/100
  X  <- mvrnorm(n = total_size, mu = rep(0, K), Sigma = diag(K))
  X_W  <- X[, 1:K_W]

  # Selection Model
  select_prob  <- inv.logit(X_W %*% coef_select)
  select_var <- rbinom(n = total_size, size = 1, prob = select_prob)

  # Potential Outcomes (for both experimental and population data)
  Y_0 <- cbind(1, X) %*% coef_0 + rnorm(n = total_size, sd = sd_Y)
  tau <- X_W %*% coef_tau + rnorm(n = total_size, sd = sd_Y)
  Y_1 <- Y_0 + tau

  # Population Data
  pop_data <-  cbind(Y_1, Y_0, X)[select_var == 0, ]
  colnames(pop_data) <- c("Y_1", "Y_0",
                          paste0("X_",seq(from  = 1, to = ncol(X))))

  # Experimental Data
  exp_size  <- sum(select_var)
  treatment <- rbinom(n = exp_size, size = 1, prob =  0.5)

  # observe outcomes
  Y <- treatment*Y_1[select_var == 1] + (1-treatment)*Y_0[select_var == 1]

  # experimental data
  exp_data <- cbind(Y, Y_1[select_var == 1], Y_0[select_var == 1],
                    treatment, X[select_var == 1,])
  colnames(exp_data) <- c("Y", "Y_1", "Y_0",
                          "treatment", paste0("X_",seq(from  = 1, to = ncol(X))))

  return(list(pop_data =  pop_data,  exp_data =  exp_data))

}
