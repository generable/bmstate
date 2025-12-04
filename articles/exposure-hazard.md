# Exposure-hazard multistate modeling with bmstate

This vignette is work in progress.

## Data simulation

### Setup

``` r
library(bmstate)
#> Attached bmstate 0.2.11. Type ?bmstate to get started.
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(tidyr)
library(tibble)
library(ggplot2)
```

``` r
# Setup
params <- list(
  N_subject = 600,
  iter_warmup = 800,
  iter_sampling = 400,
  chains = 1
)

# Names of event states
event_state_names <- function(model) {
  tm <- model$system$tm()
  tm$states_df() |>
    filter(!source) |>
    pull(state)
}
```

### True data-generating model

``` r
# True beta
create_true_covariate_effects <- function(mod) {
  NTT <- mod$system$num_states() - 1
  C <- length(mod$covs())
  bh_true <- matrix(0, NTT, C)
  bh_true[1, 2] <- 1
  bh_true[2, 2] <- -1
  bh_true[3, 1] <- 0.3
  sn <- event_state_names(mod)
  rownames(bh_true) <- paste0("Effect on ", sn)
  colnames(bh_true) <- mod$covs()
  df <- data.frame(bh_true) |>
    rownames_to_column("event") |>
    pivot_longer(cols = -event, names_to = "covariate", values_to = "beta")
  beta_pk <- list(CL = c(0.3, -0.3), V2 = c(0.3))
  list(df = df, matrix = bh_true, pk = beta_pk)
}

# True baseline hazard parameters
create_true_baseline_hazard <- function(mod) {
  tm <- mod$system$tm()

  # Spline weights
  w_true <- matrix(0, 7, 9)
  for (j in 1:7) {
    ww <- rep(0, 9)
    if (tm$trans_df()$trans_type[j] == 2) {
      ww <- rep(-1, 9)
      ww[3:5] <- 1
    }
    if (tm$trans_df()$trans_type[j] == 3) {
      ww <- rep(-1, 9)
      ww[4:8] <- 1
    }
    w_true[j, ] <- ww
  }

  # Intercept
  w0_true <- 0.5 * 1e-3
  w0_true_vec <- rep(w0_true, 7)
  w0_true_vec[3] <- 0.1 * w0_true
  w0_true_vec[5] <- 5 * w0_true
  w0_true_vec[7] <- 20 * w0_true

  # Return
  list(w0 = w0_true_vec, w = w_true)
}

# True data-generating model
create_true_model <- function() {
  # Create models
  sn <- c("Healthy", "Bleed", "Stroke", "Dead")
  tm <- transmat_diamond(state_names = sn)
  t3yr <- 3 * 365.25
  haz_covs <- c("age")
  pk_covs <- list(
    CL = c("CrCL", "age"), V2 = "weight"
  )
  create_msm(
    tm,
    hazard_covs = haz_covs, pk_covs = pk_covs, num_knots = 8, t_max = t3yr
  )
}

# Create oracle fit draws
create_oracle_fit <- function(similar_fit, beta_true, h0_true) {
  checkmate::assert_true(similar_fit$is_point_estimate())
  weights <- similar_fit$get_draws("weights")
  weights[1, , ] <- h0_true$w
  log_w0 <- similar_fit$get_draws("log_w0")
  log_w0[1, ] <- log(h0_true$w0)

  beta_oth <- similar_fit$get_draws("beta_oth")
  beta_oth[1, , ] <- t(beta_true$matrix[, 1])
  beta_oth <- posterior::rvar(beta_oth)

  beta_auc <- similar_fit$get_draws("beta_auc")
  beta_auc[1, , ] <- t(beta_true$matrix[, 2])

  lp <- similar_fit$get_draws("lp__")
  lp[] <- NA
  log_z_pk <- similar_fit$get_draws("log_z_pk")
  log_z_pk[, , ] <- 0
  log_sig_pk <- similar_fit$get_draws("log_sig_pk")
  log_sig_pk[, ] <- 0
  log_mu_pk <- similar_fit$get_draws("log_mu_pk")
  log_mu_pk[, ] <- 0
  sigma_pk <- similar_fit$get_draws("sigma_pk")
  sigma_pk[] <- 0.3
  beta_CL <- similar_fit$get_draws("beta_CL")
  beta_CL[1, ] <- beta_true$pk$CL
  beta_V2 <- similar_fit$get_draws("beta_V2")
  beta_V2[1, ] <- beta_true$pk$V2

  # Create fit
  draws <- list(
    beta_oth = beta_oth,
    beta_auc = beta_auc,
    beta_CL = beta_CL,
    beta_V2 = beta_V2,
    weights = weights,
    log_z_pk = log_z_pk,
    log_sig_pk = log_sig_pk,
    log_mu_pk = log_mu_pk,
    sigma_pk = sigma_pk,
    log_w0 = log_w0,
    lp__ = lp
  )
  sd <- similar_fit$get_data()
  mod <- similar_fit$model
  MultistateModelFit$new(similar_fit$data, sd, mod, draws,
    info = "Oracle fit"
  )
}
```

### Data simulation

``` r
mod_true <- create_true_model()
beta_true <- create_true_covariate_effects(mod_true)
h0_true <- create_true_baseline_hazard(mod_true)
simdat <- mod_true$simulate_data(
  params$N_subject,
  beta_haz = beta_true$matrix,
  beta_pk = beta_true$pk,
  w0 = h0_true$w0,
  w = h0_true$w
)
#> Recompiling Stan model
#> Using stan file at /home/runner/work/_temp/Library/bmstate/stan/msm.stan
#> Warning in check_normalized_covariate(x_norm, "ss_auc"): Normalized ss_auc has
#> maximum absolute value 11.11106, are you sure normalization of covariates is
#> correct?
#> Generating 600 paths
covs_dh <- unique(c(mod_true$data_covs(), "dose_amt"))
simdat_dh <- simdat$paths$subset_covariates(covs_dh, renamed_old = "dose", renamed_new = "dose_amt")
simdat_death <- as_single_event(simdat_dh, "Dead", null_state = "Healthy")

sa <- simdat$paths$subject_df$ss_auc
mod_true$set_auc_normalizers(loc = mean(sa), scale = stats::sd(sa))
#> setting auc normalizers to loc = 2071.2843, scale = 1770.30317

simdat <- mod_true$simulate_data(
  params$N_subject,
  beta_haz = beta_true$matrix,
  beta_pk = beta_true$pk,
  w0 = h0_true$w0,
  w = h0_true$w
)
#> Generating 600 paths
```

``` r
print(simdat)
#> A JointData object: 
#> PathData object with 600 paths
#>  * States = {Healthy, Bleed, Stroke, Dead}
#>  * Covariates = {age, CrCL, weight, dose, t_pre, t_post, conc_pre, conc_post, ss_auc, ka, CL, V2, pk_lloq}
#> 
#> A DosingData object with 600 subjects
```

## Modeling

### Defining models

``` r
NK <- 4
tm <- mod_true$system$tm()
pk_covs <- list(
  ka = mod_true$pk_model$ka_covs(),
  CL = mod_true$pk_model$CL_covs(),
  V2 = mod_true$pk_model$V2_covs()
)

# Exposure-hazard multistate model
mod_ms_eh <- create_msm(
  tm,
  hazard_covs = setdiff(mod_true$covs(), "ss_auc"),
  pk_covs = pk_covs,
  num_knots = NK, t_max = mod_true$get_tmax()
)

# Dose-hazard multistate model
mod_ms_dh <- create_msm(
  tm,
  hazard_covs = covs_dh, num_knots = NK, t_max = mod_true$get_tmax()
)

# Survival model
tm0 <- transmat_survival(tm$states[c(1, 4)])
mod_death <- create_msm(
  tm0,
  hazard_covs = covs_dh, num_knots = NK, t_max = mod_true$get_tmax()
)
```

``` r
# Inference model knots
t3yr <- mod_true$get_tmax()
tt1 <- simdat$paths$transition_times()
tt2 <- simdat_dh$transition_times()
tt3 <- simdat_death$transition_times()
mod_ms_eh$set_knots(t3yr, tt1, NK)
mod_ms_dh$set_knots(t3yr, tt2, NK)
mod_death$set_knots(t3yr, tt3, NK)
```

``` r
# Oracle fit
fit_prelim <- fit_stan(mod_true, simdat, method = "optimize", init = 0, iter = 100)
#> Shortest time interval (0.101164891162739) is smaller than delta_grid (1.09575). Consider increasing n_grid or decreasing t_max of the model.
#> Using stan file at /home/runner/work/_temp/Library/bmstate/stan/msm.stan
#> setting auc normalizers to loc = 1946.42405, scale = 1030.76555
#> setting max conc = 7812.11765
#> Initial log joint probability = -913364 
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#> Error evaluating model log probability: Non-finite gradient. 
#> Error evaluating model log probability: Non-finite gradient. 
#>       99        -10054      0.533178       236.108           1           1      119    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      100      -10047.7      0.189833       130.448       0.955       0.955      120    
#> Optimization terminated normally:  
#>   Maximum number of iterations hit, may not be at an optima 
#> Finished in  0.7 seconds.
fit_true <- create_oracle_fit(fit_prelim, beta_true, h0_true)
fit_true$covariate_effects()
#>   covariate      beta target_state_idx target_state
#> 1       age  0.0 ± NA                2        Bleed
#> 2       age  0.0 ± NA                3       Stroke
#> 3       age  0.3 ± NA                4         Dead
#> 4    ss_auc  1.0 ± NA                2        Bleed
#> 5    ss_auc -1.0 ± NA                3       Stroke
#> 6    ss_auc  0.0 ± NA                4         Dead
fit_true$plot_pk()
#> PK simulation
```

![](exposure-hazard_files/figure-html/unnamed-chunk-8-1.png)

``` r
fit_true$plot_h0()
```

![](exposure-hazard_files/figure-html/unnamed-chunk-8-2.png)

### Fitting various models

``` r
fit_ms_eh <- fit_stan(mod_ms_eh, simdat, method = "optimize", init = 0)
#> Shortest time interval (0.101164891162739) is smaller than delta_grid (1.09575). Consider increasing n_grid or decreasing t_max of the model.
#> Using stan file at /home/runner/work/_temp/Library/bmstate/stan/msm.stan
#> setting auc normalizers to loc = 1946.42405, scale = 1030.76555
#> setting max conc = 7812.11765
#> Initial log joint probability = -913214 
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#> Error evaluating model log probability: Non-finite gradient. 
#> Error evaluating model log probability: Non-finite gradient. 
#>       99        -11598      0.111203       1284.24      0.4398      0.4398      127    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      199        -10650     0.0147766       3142.56       0.231       0.231      245    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      299      -10119.2      0.102406       1843.96           1           1      354    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      399      -9700.65     0.0939468       1590.68           1           1      479    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      499      -9210.67     0.0359452       1998.75      0.8871      0.8871      587    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      599      -8661.82      0.128189       3063.91           1           1      695    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      699      -8217.94      0.217147       607.564           1           1      809    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      799      -7993.53     0.0402209       867.503           1           1      924    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      899      -7817.11    0.00780137       2229.05      0.6751      0.6751     1032    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      999      -7670.43     0.0109006       653.835           1           1     1147    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1099      -7557.53     0.0442656        1474.1       0.339           1     1263    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1199      -7488.16     0.0392131        1038.8           1           1     1372    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1299      -7402.29    0.00842692       790.551           1           1     1483    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1399      -7301.99    0.00629007       654.448           1           1     1594    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1499      -7212.89     0.0385246       980.136           1           1     1703    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1599      -7153.81     0.0102807       561.833           1           1     1810    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1699      -7105.35    0.00122472       307.033           1           1     1927    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1799      -7075.92     0.0808023       2341.81           1           1     2039    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1899      -7028.52    0.00288425       536.327           1           1     2148    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1999      -6998.48    0.00666188       317.419           1           1     2261    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     2000      -6998.35    0.00343086       303.829           1           1     2262    
#> Optimization terminated normally:  
#>   Maximum number of iterations hit, may not be at an optima 
#> Finished in  11.7 seconds.
fit_ms_dh <- fit_stan(mod_ms_dh, simdat_dh, method = "optimize", init = 0)
#> Shortest time interval (0.0999999999999943) is smaller than delta_grid (1.09575). Consider increasing n_grid or decreasing t_max of the model.
#> Using stan file at /home/runner/work/_temp/Library/bmstate/stan/msm.stan
#> Initial log joint probability = -823515 
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>       99      -6166.77      0.324794       95.4285           1           1      127    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      199      -6121.13     0.0145094       18.6248           1           1      238    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      299      -6093.87     0.0128702       26.3508      0.3991      0.3991      346    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      399      -6081.71    0.00929491       10.9813           1           1      456    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      499      -6080.61    0.00239326       8.66947           1           1      566    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      599      -6080.33    0.00493426       5.41942           1           1      672    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      699      -6080.04     0.0457803       8.99354           1           1      782    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      799      -6079.82     0.0244148       12.3921      0.2345           1      891    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      899      -6079.65    0.00134296       1.55655           1           1     1002    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      999      -6079.59    0.00514431       2.65507           1           1     1119    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1099      -6079.57    0.00460345       1.78312           1           1     1234    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1199      -6079.54     0.0105793       2.77461           1           1     1342    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1299      -6079.36      0.118515       9.08166           1           1     1451    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1399      -6079.13   0.000727898        3.0097           1           1     1558    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1499      -6079.04   0.000239679       1.89155      0.6807      0.6807     1666    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1555      -6079.03   0.000894302      0.452686           1           1     1728    
#> Optimization terminated normally:  
#>   Convergence detected: relative gradient magnitude is below tolerance 
#> Finished in  5.2 seconds.
fit_death <- fit_stan(mod_death, simdat_death, method = "optimize", init = 0)
#> Using stan file at /home/runner/work/_temp/Library/bmstate/stan/msm.stan
#> Initial log joint probability = -290039 
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>       99       -3436.6      0.013958       2.28045           1           1      119    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      124      -3436.51    0.00200116     0.0855125      0.8082      0.8082      149    
#> Optimization terminated normally:  
#>   Convergence detected: relative gradient magnitude is below tolerance 
#> Finished in  0.1 seconds.
fit_ms_eh$plot_h0()
```

![](exposure-hazard_files/figure-html/unnamed-chunk-9-1.png)

``` r
fit_ms_eh$covariate_effects()
#>   covariate       beta target_state_idx target_state
#> 1       age -0.16 ± NA                2        Bleed
#> 2       age -1.03 ± NA                3       Stroke
#> 3       age -0.17 ± NA                4         Dead
#> 4    ss_auc  3.33 ± NA                2        Bleed
#> 5    ss_auc  3.25 ± NA                3       Stroke
#> 6    ss_auc  2.13 ± NA                4         Dead
fit_ms_dh$covariate_effects()
#>    covariate        beta target_state_idx target_state
#> 1        age  0.966 ± NA                2        Bleed
#> 2        age -0.415 ± NA                3       Stroke
#> 3        age  0.264 ± NA                4         Dead
#> 4       CrCL -0.775 ± NA                2        Bleed
#> 5       CrCL  0.273 ± NA                3       Stroke
#> 6       CrCL -0.062 ± NA                4         Dead
#> 7     weight -0.858 ± NA                2        Bleed
#> 8     weight  0.301 ± NA                3       Stroke
#> 9     weight -0.034 ± NA                4         Dead
#> 10  dose_amt  1.451 ± NA                2        Bleed
#> 11  dose_amt -0.749 ± NA                3       Stroke
#> 12  dose_amt  0.047 ± NA                4         Dead
fit_death$covariate_effects()
#>   covariate        beta target_state_idx target_state
#> 1       age  0.029 ± NA                2         Dead
#> 2      CrCL -0.016 ± NA                2         Dead
#> 3    weight  0.032 ± NA                2         Dead
#> 4  dose_amt -0.016 ± NA                2         Dead
```

``` r
fit_ms_eh$plot_pk()
#> PK simulation
```

![](exposure-hazard_files/figure-html/unnamed-chunk-10-1.png)
