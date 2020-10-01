
exp_size <- 5000 # the target size of the experimental data
prob_exp <- runif(n = nrow(exp_data), min = 0.01, max = 0.98) # this should be replaced by the estimated sampling probability for each of the experimental data.

mean_prob <- mean(prob_exp) # compute its mean
exp_id <- sample(seq(from = 1, to = nrow(exp_data)),
                 size = exp_size*(1/mean_prob), replace = TRUE) # bootstrap the experimental data
sample_var <- rbinom(n = length(exp_id), size = 1, prob = prob_exp[exp_id])
use_exp <- exp_id[sample_var == 1]
boot_exp <- exp_data[use_exp, , drop = FALSE]
