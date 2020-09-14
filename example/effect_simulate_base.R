
# Simulation Main Code
rm(list = ls())
library(estimatr)
source("effect_dgp.R")

K <- 5
sim <- 1000
total_size <- 1000
# Simulate Data
dgp <-  dgp_exp_pop(total_size = total_size,
                    K = K,
                    K_W = 5,
                    coef_select = rep(1, 5),
                    coef_tau = rep(1, 5),
                    seed_use = 1000, sim_index = sim,
                    sd_Y = 0.05)

exp_data <- as.data.frame(dgp$exp_data)
pop_data <- as.data.frame(dgp$pop_data)

# True Pate
pate_true <- mean(pop_data$Y_1 - pop_data$Y_0)


# Check function
formula_outcome <- Y ~ treatment + X_1 + X_2 + X_3 + X_4 + X_5

formula_outcome <- Y ~ X_1*X_2 + X_3 + X_4 + X_5
update(formula_outcome, paste("~ . +", treatment))

update(formula_outcome, paste("~ . -", treatment))

formula_weights <- S ~ X_1 + X_2 + X_3 + X_4 + X_5
treatment <- "treatment"
exp_data <- exp_data
pop_data <- pop_data
type = "weighting"

# SATE
