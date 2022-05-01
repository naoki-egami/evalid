evalid: Methods for Improving External Validity of Randomized Experiments
=========================================================================

**Description:**

R package `evalid` provides two key functions. (1) `tpate`: three
classes of estimators for the target population average treatment effect
(T-PATE) to conduct effect-generalization. (2) `pct`: a partial
conjunction test to conduct a sign-generalization test. Detailes of the
methods are described in Egami and Hartman (2022+).

**Authors:**

-   [Naoki Egami](https://naokiegami.com/)
-   [Erin Hartman](https://erinhartman.com/)

**References:**

-   Egami and Hartman. (2022+). [Elements of External Validity:
    Framework, Design, and
    Analysis.](https://naokiegami.com/paper/external_full.pdf)
    Conditionally accepted at *American Political Science Review*

Installation Instructions
-------------------------

You can install the most recent development version using the `devtools`
package. First you have to install `devtools` using the following code.
Note that you only have to do this once:

``` r
if(!require(devtools)) install.packages("devtools")
```

Then, load `devtools` and use the function `install_github()` to install
`factorEx`:

``` r
library(devtools)
install_github("naoki-egami/evalid", dependencies = TRUE)
```

Examples
--------

-   **Effect-Generalization**
    -   Weighting-based Estimator: `est_type` is `ipw` or `wls`
    -   Outcome-based Estimator: `est_type` is `outcome-ols` or
        `outcome-bart`
    -   Doubly Robust Estimator: `est_type` is `dr-ols` or `dr-bart`
-   **Sign-Generalization**
    -   How to use `pct()` function

Effect-Generalization via `tpate()`
-----------------------------------

The goal is to estimate the target population average treatment effect
(TPATE)

### Prepare data

When using marginal distributions, `target_dist` should be a list and
each element should have a factor name. Within each list, a `numeric`
vector should have the same level names as those in `data`.

``` r
## Load the package and data
library(evalid)
data("BroockmanKalla")

exp_data <- BroockmanKalla$exp_data # experimental data
pop_data <- BroockmanKalla$pop_data # target population data

# covariates measured in both experimental and population data
covariates <- c("vf_female", "vf_black", "vf_white", 
                "religious_t0_factor_More_than_once_a_week", "religious_t0_factor_Never", 
                "religious_t0_factor_Once_a_week", "religious_t0_factor_Once_or_twice_a_month", 
                "ideology_t0_factor_Conservative", "ideology_t0_factor_Liberal",
                "ideology_t0_factor_Moderate", "ideology_t0_factor_Very_conservative", 
                "pid_t0_factor_Independent", "pid_t0_factor_Lean_Democrat", 
                "pid_t0_factor_Lean_Republican", "pid_t0_factor_Not_very_strong_Democrat",
                "pid_t0_factor_Not_very_strong_Republican", "pid_t0_factor_Strong_Democrat", 
                "vf_age_bucket_a_18to34", "vf_age_bucket_b_35to49", "vf_age_bucket_c_50to64")
```

### Weighting-based Estimators

We can estimate the pAMCE with `design_pAMCE` with
`target_type = "marginal"`. Use `factor_name` to specify for which
factors we estimate the pAMCE.

``` r
## IPW 
ipw <- tpate(outcome = "trans.tolerance.dv.t1",
             treatment = "treat_ind", 
             covariates = covariates, 
             exp_data = exp_data, 
             pop_data = pop_data,
             id_cluster = exp_data$cluster, # cluster-standard errors 
             est_type = "ipw", # type of estimator 
             weights_max = 10) 

unlist(ipw$tpate)
```

    ##            est             se  ci_lower.2.5% ci_upper.97.5% 
    ##     0.42286638     0.22421762     0.01416949     0.88733270

``` r
## Weighted Least Squares (can include covariates measured only in the experiment)
covariates_exp <- c("miami_trans_law_t0", "miami_trans_law2_t0",
                    "therm_trans_t0", "gender_norms_sexchange_t0", "gender_norms_moral_t0", "gender_norms_abnormal_t0",
                    "ssm_t0", "therm_obama_t0", "therm_gay_t0", "vf_democrat",
                    "exposure_gay_t0", "exposure_trans_t0", "sdo_scale",
                    "gender_norm_daugher_t0", "gender_norm_looks_t0",
                    "gender_norm_rights_t0", "therm_afams_t0", "vf_hispanic",
                    "survey_language_es","cluster_level_t0_scale_mean")
wls <- tpate(outcome = "trans.tolerance.dv.t1",
             treatment = "treat_ind", 
             covariates = covariates, # covariates measured in both the experimental and population data
             covariates_exp = covariates_exp, # covariates measured only in the experimental data
             exp_data = exp_data, 
             pop_data = pop_data,
             id_cluster = exp_data$cluster, # cluster-standard errors 
             est_type = "wls", # type of estimator 
             weights_max = 10) 
unlist(wls$tpate)
```

    ##            est             se  ci_lower.2.5% ci_upper.97.5% 
    ##     0.25118647     0.08473368     0.10046551     0.41433086

### Outcome-based Estimators

We can estimate the pAMCE with `design_pAMCE` with
`target_type = "marginal"`. Use `factor_name` to specify for which
factors we estimate the pAMCE.

``` r
## OLS Outcome-based Estimator 
outcome_ols <- tpate(outcome = "trans.tolerance.dv.t1",
                     treatment = "treat_ind", 
                     covariates = covariates, # covariates measured in both the experimental and population data
                     exp_data = exp_data, 
                     pop_data = pop_data,
                     id_cluster = exp_data$cluster, # cluster-standard errors 
                     est_type = "outcome-ols") # type of estimator 
unlist(outcome_ols$tpate)
```

    ##            est             se  ci_lower.2.5% ci_upper.97.5% 
    ##     0.23738683     0.17278517    -0.09360796     0.56412235

``` r
## BART Outcome-based Estimator 
outcome_bart <- tpate(outcome = "trans.tolerance.dv.t1",
                     treatment = "treat_ind", 
                     covariates = covariates, # covariates measured in both the experimental and population data
                     exp_data = exp_data, 
                     pop_data = pop_data,
                     id_cluster = exp_data$cluster, # cluster-standard errors 
                     est_type = "outcome-bart") # type of estimator
```

    ## fitting response model via method 'bart'

``` r
unlist(outcome_bart$tpate)
```

    ##            est             se  ci_lower.2.5% ci_upper.97.5% 
    ##     0.12580007     0.08896522    -0.04532356     0.29275451

### Doubly Robust Estimators

We can estimate the pAMCE with `design_pAMCE` with
`target_type = "marginal"`. Use `factor_name` to specify for which
factors we estimate the pAMCE.

``` r
## Doubley Robust Estimator (outcome model is OLS)
dr_ols <- tpate(outcome = "trans.tolerance.dv.t1",
             treatment = "treat_ind", 
             covariates = covariates, 
             exp_data = exp_data, 
             pop_data = pop_data,
             id_cluster = exp_data$cluster, # cluster-standard errors 
             est_type = "dr-ols", # type of estimator 
             weights_max = 10) 
unlist(dr_ols$tpate)
```

    ##            est             se  ci_lower.2.5% ci_upper.97.5% 
    ##     0.35433831     0.18752754    -0.01366399     0.72158388

``` r
## Doubley Robust Estimator (outcome model is BART)
dr_bart <- tpate(outcome = "trans.tolerance.dv.t1",
             treatment = "treat_ind", 
             covariates = covariates, 
             exp_data = exp_data, 
             pop_data = pop_data,
             id_cluster = exp_data$cluster, # cluster-standard errors 
             est_type = "dr-bart", # type of estimator 
             weights_max = 10) 
```

    ## fitting response model via method 'bart'

``` r
unlist(dr_bart$tpate)
```

    ##            est             se  ci_lower.2.5% ci_upper.97.5% 
    ##     0.34305283     0.16520622     0.03071198     0.68082746

Sign-Generalization via `pct()`
-------------------------------

### When a substantive theory predicts a positive causal effect (i.e., null hypothesis: a causal effect is less than or equal to zero)

``` r
library(estimatr)
data("Bisgaard")

## Opposition Supporters (The Bisgaard theory predicts positive effects for this subgroup)
opp_data <- Bisgaard$opp_data  ## cleaned data for 12 different variations 

## Compute one-sided p-values for each variation
pv_opp <- c()
for(i in 1:12){
  data_variation <- opp_data[[i]]
  fit_variation <- lm_robust(outcome ~ treatment, data = data_variation)  # compute the causal effect for each variation
  z <- summary(fit_variation)$coef["treatment", "Estimate"]/summary(fit_variation)$coef["treatment", "Std. Error"]
  pv_opp[i] <- 1 - pnorm(z)  # one-sided p-value for positive effects
}

## Partial Conjunction Test
pct(pv_opp)
```

    ##    threshold      p_value h_num
    ## 1          1 5.595524e-14    11
    ## 2          2 7.083223e-14     1
    ## 3          3 5.805023e-11    12
    ## 4          4 1.820434e-10     8
    ## 5          5 4.434087e-07     4
    ## 6          6 5.739760e-05    10
    ## 7          7 5.739760e-05     7
    ## 8          8 9.520086e-05     9
    ## 9          9 1.897268e-04     2
    ## 10        10 4.581073e-04     6
    ## 11        11 4.581073e-04     3
    ## 12        12 2.207919e-01     5

### When a substantive theory predicts a negative causal effect (i.e., null hypothesis: a causal effect is larger than or equal to zero)

``` r
## Government Supporters (The Bisgaard theory predicts negative effects for this subgroup)
gov_data <- Bisgaard$gov_data  ## cleaned data for 12 different variations 

## Compute one-sided p-values for each variation
pv_gov <- c()
for(i in 1:12){
  data_variation <- gov_data[[i]]
  fit_variation <- lm_robust(outcome ~ treatment, data = data_variation)  # compute the causal effect for each variation
  z <- summary(fit_variation)$coef["treatment", "Estimate"]/summary(fit_variation)$coef["treatment", "Std. Error"]
  pv_gov[i] <- pnorm(z)  # one-sided p-value for negative effects
}

## Partial Conjunction Test
pct(pv_gov)
```

    ##    threshold      p_value h_num
    ## 1          1 2.555573e-27     9
    ## 2          2 1.111781e-22    11
    ## 3          3 9.092034e-17     1
    ## 4          4 1.561087e-14     8
    ## 5          5 2.790215e-12     4
    ## 6          6 4.467047e-10     7
    ## 7          7 5.255725e-08    12
    ## 8          8 3.178489e-03     2
    ## 9          9 1.515678e-01     5
    ## 10        10 2.475779e-01     6
    ## 11        11 2.513056e-01    10
    ## 12        12 2.513056e-01     3
