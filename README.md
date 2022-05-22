evalid: Methods for Improving External Validity of Randomized Experiments
=========================================================================

**Description:**

R package `evalid` provides two key functions. (1) `tpate()`: three
classes of estimators for the target population average treatment effect
(T-PATE) to conduct effect-generalization. (2) `pct()`: a partial
conjunction test to conduct a sign-generalization test. Details of the
methods are described in Egami and Hartman (2022+).

**Authors:**

-   [Naoki Egami](https://naokiegami.com/)
-   [Erin Hartman](https://erinhartman.com/)

**Reference:**

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
`evalid`:

``` r
library(devtools)
install_github("naoki-egami/evalid", dependencies = TRUE)
```

Overview
--------

Please read Section “The Proposed Approach toward External Validity:
Outline” in Egami and Hartman (2022+) for the overview of effect- and
sign-generalization.

**Effect-Generalization** 
- Weighting-based Estimator: IPW estimator
(`ipw`) or Weighted Least Squares (`wls`) 
- Outcome-based Estimator:
OLS-based estimator (`outcome-ols`) or BART-based estimator
(`outcome-bart`) 
- Doubly Robust Estimator: Doubly robust estimator with
OLS-based outcome model (`dr-ols`) or Doubly robust estimator with
BART-based outcome model (`dr-bart`)

**Sign-Generalization** 
- How to use `pct()` function

Effect-Generalization via `tpate()`
-----------------------------------

The goal is to estimate the target population average treatment effect
(T-PATE) in the target population data. Please read Section
“Effect-Generalization” in Egami and Hartman (2022+) for details and
practical recommendations for effect-generalization. We use Broockman
and Kalla (2018) as an example. See Section “Empirical Applications” in
Egami and Hartman (2022+) as well.

### Prepare data

We need two data sets. (1) `exp_data`: the experimental data that
includes outcome, treatment, and covariates we adjust for. (2)
`pop_data`: the target population data that includes covariates we
adjust for.

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

Weighting-based estimators (`ipw` and `wls`) estimate the T-PATE via
estimating weights that balance covariates in the experimental data and
population data. In addition to the ignorability of sampling and
treatment effect heterogeneity, this class of estimators assumes weights
are correctly modeled. Importantly, weighted least squares estimators
(`wls`) can incorporate covariates measured only in the experiment to
increase efficiency (using argument `covariates_exp`).

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

Outcome-based estimators (`outcome-ols` and `outcome-bart`) estimate the
T-PATE via modeling outcome regression. In addition to the ignorability
of sampling and treatment effect heterogeneity, this class of estimators
assumes outcome models are correctly specified.

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
    ##     0.12958598     0.08913481    -0.03087554     0.30615813

### Doubly Robust Estimators

Doubly robust estimators (`dr-ols` and `dr-bart`) estimate the T-PATE
via modeling both weights and outcome regression. Doubly robust
estimators are consistent to the T-PATE if one of the two (weights or an
outcome model) is correctly specified. Note that this estimator also
makes the ignorability of sampling and treatment effect heterogeneity
assumption for causal identification.

``` r
## Doubly Robust Estimator (outcome model is OLS)
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
## Doubly Robust Estimator (outcome model is BART)
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
    ##     0.34082774     0.17697526    -0.01074015     0.67573733

Sign-Generalization via `pct()`
-------------------------------

The goal is to test the sign of the target population average treatment
effect (T-PATE) using design of purposive variations within an
experiment. Please read Section “Sign-Generalization” in Egami and
Hartman (2022+) for details and practical recommendations for
sign-generalization. We use Bisgaard (2019) as an example. See Section
“Empirical Applications” in Egami and Hartman (2022+) as well.

### Prepare data

We need to compute the SATE for each purposive variation. To facilitate
this process, we created a cleaned data set for each purposive variation
(in total, 12 for opposition supporters and government supporters,
respectively).

``` r
data("Bisgaard")

## Opposition Supporters (The Bisgaard theory predicts the positive T-PATE for this subgroup)
opp_data <- Bisgaard$opp_data  ## cleaned data for 12 different variations 

head(opp_data[[1]][, c("outcome", "treatment")])
```

    ##   outcome treatment
    ## 1    0.50         1
    ## 2    0.50         1
    ## 3    1.00         1
    ## 4    0.50         1
    ## 5    0.50         1
    ## 6    0.75         1

``` r
## Government Supporters (The Bisgaard theory predicts the negative T-PATE for this subgroup)
gov_data <- Bisgaard$gov_data  ## cleaned data for 12 different variations 

head(gov_data[[1]][, c("outcome", "treatment")])
```

    ##   outcome treatment
    ## 1    1.00         1
    ## 2    0.50         1
    ## 3    0.50         0
    ## 4    0.25         1
    ## 5    0.50         0
    ## 6    0.50         1

### Testing whether the T-PATE is positive

First, we consider cases when a substantive theory predicts the T-PATE
is positive. Therefore, a null hypothesis is that the T-PATE is less
than or equal to zero. Therefore, when we compute one-sided p-values, we
use `1 - prnom(z)` where `z` is the ratio of an estimated causal effect
and its standard error.

``` r
library(estimatr)
## Compute one-sided p-values for each variation
pv_opp <- c()
for(i in 1:12){
  data_variation <- opp_data[[i]]
  fit_variation <- lm_robust(outcome ~ treatment, data = data_variation)  # compute the SATE for each variation
  z <- summary(fit_variation)$coef["treatment", "Estimate"]/summary(fit_variation)$coef["treatment", "Std. Error"]
  pv_opp[i] <- 1 - pnorm(z)  # one-sided p-value for positive effects
}

## Partial Conjunction Test
round(pct(pv_opp), 2)
```

    ##    threshold p_value h_num
    ## 1          1    0.00    11
    ## 2          2    0.00     1
    ## 3          3    0.00    12
    ## 4          4    0.00     8
    ## 5          5    0.00     4
    ## 6          6    0.00    10
    ## 7          7    0.00     7
    ## 8          8    0.00     9
    ## 9          9    0.00     2
    ## 10        10    0.00     6
    ## 11        11    0.00     3
    ## 12        12    0.22     5

### Testing whether the T-PATE is negative

Next, we consider cases when a substantive theory predicts the T-PATE is
negative. Therefore, a null hypothesis is that the T-PATE is larger than
or equal to zero. Therefore, when we compute one-sided p-values, we use
`prnom(z)` where `z` is the ratio of an estimated causal effect and its
standard error.

``` r
## Compute one-sided p-values for each variation
pv_gov <- c()
for(i in 1:12){
  data_variation <- gov_data[[i]]
  fit_variation <- lm_robust(outcome ~ treatment, data = data_variation)  # compute the SATE for each variation
  z <- summary(fit_variation)$coef["treatment", "Estimate"]/summary(fit_variation)$coef["treatment", "Std. Error"]
  pv_gov[i] <- pnorm(z)  # one-sided p-value for negative effects
}

## Partial Conjunction Test
round(pct(pv_gov), 2)
```

    ##    threshold p_value h_num
    ## 1          1    0.00     9
    ## 2          2    0.00    11
    ## 3          3    0.00     1
    ## 4          4    0.00     8
    ## 5          5    0.00     4
    ## 6          6    0.00     7
    ## 7          7    0.00    12
    ## 8          8    0.00     2
    ## 9          9    0.15     5
    ## 10        10    0.25     6
    ## 11        11    0.25    10
    ## 12        12    0.25     3
