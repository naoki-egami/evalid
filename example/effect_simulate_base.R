
# Simulation Main Code
rm(list = ls())d
source("example/effect_dgp.R")

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

formula_outcome <- Y ~ treatment + X_1 + X_2 + X_3 + X_4 + X_5
formula_weights <- S ~ X_1 + X_2 + X_3 + X_4 + X_5


# 1. Weighting-based Estimator
ipw_logit <- tpate(formula_outcome = formula_outcome,
                   formula_weights = formula_weights,
                   exp_data = exp_data, pop_data = pop_data,
                   treatment = "treatment",
                   est_type = "ipw", weights_type = "logit",
                   sims = 1000)
ipw_cal <- tpate(formula_outcome = formula_outcome,
                   formula_weights = formula_weights,
                   exp_data = exp_data, pop_data = pop_data,
                   treatment = "treatment",
                   est_type = "ipw", weights_type = "calibration",
                   sims = 1000)

wls_logit <- tpate(formula_outcome = formula_outcome,
                   formula_weights = formula_weights,
                   exp_data = exp_data, pop_data = pop_data,
                   treatment = "treatment",
                   est_type = "wls", weights_type = "logit",
                   sims = 1000)
wls_cal <- tpate(formula_outcome = formula_outcome,
                 formula_weights = formula_weights,
                 exp_data = exp_data, pop_data = pop_data,
                 treatment = "treatment",
                 est_type = "wls", weights_type = "calibration",
                 sims = 1000)

# 2. Outcome-based Estimator
out_ols <- tpate(formula_outcome = formula_outcome,
                 formula_weights = formula_weights,
                 exp_data = exp_data, pop_data = pop_data,
                 treatment = "treatment",
                 est_type = "outcome-ols",
                 sims = 1000)
out_bart <- tpate(formula_outcome = formula_outcome,
                  formula_weights = formula_weights,
                  exp_data = exp_data, pop_data = pop_data,
                  treatment = "treatment",
                  est_type = "outcome-bart",
                  sims = 1000)

# 2. Outcome-based Estimator
out_ols <- tpate(formula_outcome = formula_outcome,
                 formula_weights = formula_weights,
                 exp_data = exp_data, pop_data = pop_data,
                 treatment = "treatment",
                 est_type = "outcome-ols",
                 sims = 1000)
out_bart <- tpate(formula_outcome = formula_outcome,
                  formula_weights = formula_weights,
                  exp_data = exp_data, pop_data = pop_data,
                  treatment = "treatment",
                  est_type = "outcome-bart",
                  sims = 1000)


# 3. DR Estimator
dr_logit_ols <- tpate(formula_outcome = formula_outcome,
                 formula_weights = formula_weights,
                 exp_data = exp_data, pop_data = pop_data,
                 treatment = "treatment", weights_type = "logit",
                 est_type = "dr-ols",
                 sims = 1000)
dr_logit_bart <- tpate(formula_outcome = formula_outcome,
                  formula_weights = formula_weights,
                  exp_data = exp_data, pop_data = pop_data,
                  treatment = "treatment", weights_type = "logit",
                  est_type = "dr-bart",
                  sims = 1000)

dr_cal_ols <- tpate(formula_outcome = formula_outcome,
                      formula_weights = formula_weights,
                      exp_data = exp_data, pop_data = pop_data,
                      treatment = "treatment", weights_type = "calibration",
                      est_type = "dr-ols",
                      sims = 1000)
dr_cal_bart <- tpate(formula_outcome = formula_outcome,
                       formula_weights = formula_weights,
                       exp_data = exp_data, pop_data = pop_data,
                       treatment = "treatment", weights_type = "calibration",
                       est_type = "dr-bart",
                       sims = 1000)

est_tab <- as.matrix(rbind(ipw_logit$tpate,
                           ipw_cal$tpate,
                           wls_logit$tpate,
                           wls_cal$tpate,
                           out_ols$tpate,
                           out_bart$tpate,
                           dr_logit_ols$tpate,
                           dr_logit_bart$tpate,
                           dr_cal_ols$tpate,
                           dr_cal_bart$tpate))

plot(seq(1:10), est_tab[,1], pch = 19, ylim = c(-2, 1), xaxt = "n")
segments(seq(1:10), unlist(est_tab[, 3]),
         seq(1:10), unlist(est_tab[, 4]))
abline(h = pate_true, lwd = 2, col = "red")
Axis(side = 1, at = seq(1:10),
     labels = c("ipw_logit",
                "ipw_cal",
                "wls_logit",
                "wls_cal",
                "out_ols",
                "out_bart",
                "dr_logit_ols",
                "dr_logit_bart",
                "dr_cal_ols",
                "dr_cal_bart"))
