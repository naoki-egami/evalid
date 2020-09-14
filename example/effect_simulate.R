
# Simulation Main Code
rm(list = ls())
library(estimatr)
source("effect_dgp.R")
source("estimator/pate_estimator.R")

total_size <- 6000
sim_size <-  100
pate_true <-  c()

sate_estimate <- matrix(NA, nrow = sim_size, ncol = 4)
pate_ipw <- matrix(NA, nrow = sim_size, ncol = 4)
pate_wLS <- matrix(NA, nrow = sim_size, ncol = 4)
pate_wLS_int <- matrix(NA, nrow = sim_size, ncol = 4)
pate_ipw_m <- matrix(NA, nrow = sim_size, ncol = 4)
pate_wLS_m <- matrix(NA, nrow = sim_size, ncol = 4)
pate_pro_r <- matrix(NA, nrow = sim_size, ncol = 4)
pate_pro_r_m <- matrix(NA, nrow = sim_size, ncol = 4)
pate_pro_lm <- matrix(NA, nrow = sim_size, ncol = 4)
pate_pro_lm_m <- matrix(NA, nrow = sim_size, ncol = 4)
pate_aipw_os <- matrix(NA, nrow = sim_size, ncol = 4)
pate_aipw_o <- matrix(NA, nrow = sim_size, ncol = 4)
pate_aipw_s <- matrix(NA, nrow = sim_size, ncol = 4)
pate_aipw_no <- matrix(NA, nrow = sim_size, ncol = 4)
pate_wls_pro_os <- matrix(NA, nrow = sim_size, ncol = 4)
pate_wls_pro_o <- matrix(NA, nrow = sim_size, ncol = 4)
pate_wls_pro_s <- matrix(NA, nrow = sim_size, ncol = 4)
pate_wls_pro_no <- matrix(NA, nrow = sim_size, ncol = 4)

K <- 5

for(sim in 1:sim_size){
  
  if(sim %% 10  == 0){
    cat(paste0(sim,"..."))  
  }
  
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
  pate_true[sim] <- mean(pop_data$Y_1 - pop_data$Y_0)
  
  # SATE estimate
  fit_sate <- lm_robust(Y ~ treatment, data = exp_data)
  sate_estimate[sim, 1:4] <- summary(fit_sate)$coefficient["treatment", c(1,2,5,6)]
  
  # ######################
  # Correct Models
  # ######################
  
  # 1.IPW 
  ipw_e <- wLS(formula_Y = Y ~ treatment, 
               formula_S = S ~ X_1 + X_2 + X_3 + X_4 + X_5, 
               treat.name = "treatment", exp_data = exp_data, 
               pop_data  = pop_data)
  pate_ipw[sim, 1:4] <- ipw_e$pate[c(1,2,5,6)]
  
  # setup formula 
  formula_wls <- as.formula(paste0("Y ~ treatment + ", 
                                   paste(paste0("X_",seq(from = 1, to = K)), 
                                         collapse = "+")))
  formula_p <- as.formula(paste0("Y ~ ", 
                                 paste(paste0("X_",seq(from = 1, to = K)), 
                                       collapse = "+")))
  # 2.wLS
  wLS_e <- wLS(formula_wls, 
               formula_S = S ~ X_1 + X_2 + X_3 + X_4 + X_5, 
               treat.name = "treatment", 
               exp_data = exp_data, 
               pop_data  = pop_data)
  pate_wLS[sim, 1:4] <- wLS_e$pate[c(1,2,5,6)]
  
  # 3. projection 
  # 3.1 ranger (random  forrests): seems not-consistent when the true model is linear!
  pro_r_pate <- project(formula_p, type = "ranger",
                        exp_data = exp_data,
                        pop_data  = pop_data)
  pate_pro_r[sim, 1] <- pro_r_pate$pate
  
  # 3.2 projection with lm
  pro_lm_e <- project(formula_p, type = "lm",
                      exp_data = exp_data,
                      pop_data  = pop_data)
  pate_pro_lm[sim, 1:4] <- pro_lm_e$pate
  
  # 4. wLS_projection (doubly robust)
  # 4.1 both correct
  wls_pro_os <- wLS_project(formula_p,
                            formula_S = S ~ X_1 + X_2 + X_3 + X_4 + X_5,
                            type = "lm",
                            exp_data = exp_data,
                            pop_data  = pop_data)
  pate_wls_pro_os[sim, 1:4] <- wls_pro_os$pate
  
  # 4.2 only weighting 
  wls_pro_s <- wLS_project(Y ~ X_1 + X_2,
                           formula_S = S ~ X_1 + X_2 + X_3 + X_4 + X_5,
                           type = "lm",
                           exp_data = exp_data,
                           pop_data  = pop_data)
  pate_wls_pro_s[sim, 1:4] <- wls_pro_s$pate
  
  # 4.3 only outcome
  wls_pro_o <- wLS_project(formula_p,
                           formula_S = S ~ X_1 + X_2,
                           type = "lm",
                           exp_data = exp_data,
                           pop_data  = pop_data)
  pate_wls_pro_o[sim, 1:4] <- wls_pro_o$pate
  
  # 5. augumented IPW (doubly robust)
  # 5.1 both are correct
  aipw_os <- aipw(formula_p,
                formula_S = S ~ X_1 + X_2 + X_3 + X_4 + X_5,
                type = "lm",
                exp_data = exp_data,
                pop_data  = pop_data)
  pate_aipw_os[sim, 1] <- aipw_os$pate
  
  # 5.2 only weighting correct
  aipw_s <- aipw(Y ~ X_1 + X_2,
               formula_S = S ~ X_1 + X_2 + X_3 + X_4 + X_5,
               type = "lm",
               exp_data = exp_data,
               pop_data  = pop_data)
  pate_aipw_s[sim, 1] <- aipw_s$pate
  
  # 5.3 only outcome correct
  aipw_o <- aipw(formula_p,
               formula_S = S ~ X_1 + X_2,
               type = "lm",
               exp_data = exp_data,
               pop_data  = pop_data)
  pate_aipw_o[sim, 1] <- aipw_o$pate
  
  
  # # ######################
  # # Misspecified Models
  # # ######################
  
  # 1. misspecified IPW
  ipw_e_m <- wLS(Y ~ treatment,
                 formula_S = S ~ X_1 + X_2,
                 treat.name = "treatment", exp_data = exp_data,
                 pop_data  = pop_data)
  pate_ipw_m[sim, 1:4] <- ipw_e_m$pate[c(1,2,5,6)]
  
  # 2. misspecified wLS
  wLS_e_m <- wLS(formula_wls,
                 formula_S = S ~ X_1 + X_2,
                 treat.name = "treatment",
                 exp_data = exp_data,
                 pop_data  = pop_data)
  pate_wLS_m[sim, 1:4] <- wLS_e_m$pate[c(1,2,5,6)]
  
  # 3. misspecified projection
  # 3.1 projection with ranger
  pro_r_e_m <- project(Y ~ X_1 + X_2, type = "ranger",
                       exp_data = exp_data,
                       pop_data  = pop_data)
  pate_pro_r_m[sim, 1] <- pro_r_e_m$pate
  
  # 3.2 projection with lm
  pro_lm_e_m <- project(Y ~ X_1 + X_2, type = "lm",
                        exp_data = exp_data,
                        pop_data  = pop_data)
  pate_pro_lm_m[sim, 1:4] <- pro_lm_e_m$pate
  
  # 4. misspecified wLS_projection (both models are wrong)
  wls_pro_no <- wLS_project(Y ~ X_1 + X_2,
                            formula_S = S ~ X_1 + X_2,
                            type = "lm",
                            exp_data = exp_data,
                            pop_data  = pop_data)
  pate_wls_pro_no[sim, 1:4] <- wls_pro_no$pate
  
  
  # 5. misspecified augumented IPW (both models are wrong)
  aipw_no <- aipw(Y ~ X_1 + X_2,
                formula_S = S ~ X_1 + X_2,
                type = "lm",
                exp_data = exp_data,
                pop_data  = pop_data)
  pate_aipw_no[sim, 1] <- aipw_no$pate
}
