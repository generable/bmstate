# Exposure-hazard multistate modeling with bmstate

This vignette is work in progress.

## Data simulation

### Setup

``` r
library(bmstate)
#> Attached bmstate 0.2.4. Type ?bmstate to get started.
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
#> maximum absolute value 13.84517, are you sure normalization of covariates is
#> correct?
#> Generating 600 paths
covs_dh <- unique(c(mod_true$data_covs(), "dose_amt"))
simdat_dh <- simdat$paths$subset_covariates(covs_dh, renamed_old = "dose", renamed_new = "dose_amt")
simdat_death <- as_single_event(simdat_dh, "Dead", null_state = "Healthy")

sa <- simdat$paths$subject_df$ss_auc
mod_true$set_auc_normalizers(loc = mean(sa), scale = stats::sd(sa))
#> setting auc normalizers to loc = 266.63989, scale = 197.01093

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
#>  * Covariates = {age, CrCL, weight, dose, t_pre, t_post, conc_pre, conc_post, ss_auc, pk_lloq}
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
#> Shortest time interval (0.100000000000023) is smaller than delta_grid (1.09575). Consider increasing n_grid or decreasing t_max of the model.
#> Using stan file at /home/runner/work/_temp/Library/bmstate/stan/msm.stan
#> setting auc normalizers to loc = 259.5406, scale = 138.62034
#> setting max conc = 15247.12757
#> Warning in max(x[!is.na(x)]): no non-missing arguments to max; returning -Inf
#> Initial log joint probability = -914138 
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#> Error evaluating model log probability: Non-finite gradient. 
#> Error evaluating model log probability: Non-finite gradient. 
#>       99      -9686.62      0.146466       1791.64      0.2262           1      125    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      100      -9673.83     0.0812532       566.721       2.225      0.2225      127    
#> Optimization terminated normally:  
#>   Maximum number of iterations hit, may not be at an optima 
#> Finished in  0.9 seconds.
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
#> Shortest time interval (0.100000000000023) is smaller than delta_grid (1.09575). Consider increasing n_grid or decreasing t_max of the model.
#> Using stan file at /home/runner/work/_temp/Library/bmstate/stan/msm.stan
#> setting auc normalizers to loc = 259.5406, scale = 138.62034
#> setting max conc = 15247.12757
#> Warning in max(x[!is.na(x)]): no non-missing arguments to max; returning -Inf
#> Initial log joint probability = -913988 
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#> Error evaluating model log probability: Non-finite gradient. 
#> Error evaluating model log probability: Non-finite gradient. 
#> Error evaluating model log probability: Non-finite gradient. 
#> Error evaluating model log probability: Non-finite gradient. 
#>       99        -11002     0.0392664       763.081      0.8973      0.8973      131    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      199         -9835     0.0732173       2597.02      0.5619      0.5619      245    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      299      -9131.95      0.202998       1832.69           1           1      359    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      399      -8766.46     0.0433002       1438.55           1           1      476    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      499      -8372.43     0.0280325       1367.53           1           1      591    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      599      -8071.13      0.207503       1727.99           1           1      705    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      699      -7756.65     0.0628896       1386.33           1           1      817    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      799      -7551.84     0.0555508       1747.72           1           1      929    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      899      -7424.55    0.00941451       734.704           1           1     1049    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      999      -7290.49     0.0748073       4297.77           1           1     1169    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1099      -7160.44     0.0526671       2531.31       2.987      0.2987     1285    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1199      -7090.24     0.0149287        1487.4       1.003      0.1003     1407    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1299      -7053.95    0.00649302       1366.36      0.9693      0.9693     1520    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1399      -7005.01    0.00876072       2297.52      0.4472      0.4472     1633    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1499      -6948.53    0.00297488       484.894           1           1     1746    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1599      -6915.76     0.0110758       2892.69      0.1403           1     1859    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1699      -6867.86    0.00243423       2726.25      0.6351      0.6351     1975    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1799      -6790.54     0.0193707       1688.32           1           1     2087    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1899      -6751.63    0.00319797        458.66           1           1     2204    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1999      -6696.34    0.00427482       1917.82      0.3566      0.3566     2319    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     2000      -6695.49     0.0094888       1097.59           1           1     2320    
#> Optimization terminated normally:  
#>   Maximum number of iterations hit, may not be at an optima 
#> Finished in  14.1 seconds.
fit_ms_dh <- fit_stan(mod_ms_dh, simdat_dh, method = "optimize", init = 0)
#> Shortest time interval (0.0999999999999943) is smaller than delta_grid (1.09575). Consider increasing n_grid or decreasing t_max of the model.
#> Using stan file at /home/runner/work/_temp/Library/bmstate/stan/msm.stan
#> Initial log joint probability = -615022 
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>       99      -6469.65     0.0607907       52.5899           1           1      113    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      199      -6443.97        0.1533       60.0679           1           1      223    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      299      -6433.51     0.0764647       29.8971           1           1      334    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      399      -6429.75    0.00358569       6.66979      0.8298      0.8298      441    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      499      -6425.74     0.0430428       27.8584           1           1      549    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      599      -6421.24     0.0355685       18.8165           1           1      661    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      699      -6420.01     0.0325662       17.6514      0.3778           1      773    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      799       -6418.5     0.0353841       11.4806      0.5365      0.5365      884    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      899      -6417.87    0.00306629       2.89229           1           1      999    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      999      -6415.76     0.0183879       19.8563           1           1     1106    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1099       -6415.1    0.00335729       6.52143       0.235           1     1219    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1199       -6414.9    0.00210733       4.94835           1           1     1327    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1299      -6414.09    0.00560676       6.56934      0.4691      0.4691     1431    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1399      -6413.89    0.00308072       2.65389           1           1     1539    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1499      -6413.85   0.000386251       2.06569      0.5557      0.5557     1645    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1599      -6413.29     0.0517785       21.4497      0.4567      0.4567     1756    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1699      -6401.02     0.0244175       13.3994           1           1     1864    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1799      -6394.95     0.0124588       10.8694           1           1     1976    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1899      -6393.59      0.117618       43.7635           1           1     2083    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1999      -6393.34    0.00212374       3.45562      0.4303           1     2192    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     2000      -6393.34    0.00244148       2.78249      0.4346      0.4346     2193    
#> Optimization terminated normally:  
#>   Maximum number of iterations hit, may not be at an optima 
#> Finished in  9.6 seconds.
fit_death <- fit_stan(mod_death, simdat_death, method = "optimize", init = 0)
#> Using stan file at /home/runner/work/_temp/Library/bmstate/stan/msm.stan
#> Initial log joint probability = -220374 
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>       99      -3580.41    0.00838396       1.22793           1           1      130    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      199      -3560.35   0.000689369      0.137888           1           1      258    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      220      -3560.35   0.000154399     0.0902147      0.6747      0.6747      279    
#> Optimization terminated normally:  
#>   Convergence detected: relative gradient magnitude is below tolerance 
#> Finished in  0.1 seconds.
fit_ms_eh$plot_h0()
```

![](exposure-hazard_files/figure-html/unnamed-chunk-9-1.png)

``` r
fit_ms_eh$covariate_effects()
#>   covariate          beta target_state_idx target_state
#> 1       age  0.27878 ± NA                2        Bleed
#> 2       age -0.00012 ± NA                3       Stroke
#> 3       age  0.32507 ± NA                4         Dead
#> 4    ss_auc  0.46823 ± NA                2        Bleed
#> 5    ss_auc -0.44588 ± NA                3       Stroke
#> 6    ss_auc -0.04961 ± NA                4         Dead
fit_ms_dh$covariate_effects()
#>    covariate        beta target_state_idx target_state
#> 1        age  1.256 ± NA                2        Bleed
#> 2        age -0.493 ± NA                3       Stroke
#> 3        age  0.405 ± NA                4         Dead
#> 4       CrCL -1.138 ± NA                2        Bleed
#> 5       CrCL  0.539 ± NA                3       Stroke
#> 6       CrCL -0.039 ± NA                4         Dead
#> 7     weight  0.040 ± NA                2        Bleed
#> 8     weight  0.036 ± NA                3       Stroke
#> 9     weight  0.035 ± NA                4         Dead
#> 10  dose_amt  1.792 ± NA                2        Bleed
#> 11  dose_amt -1.086 ± NA                3       Stroke
#> 12  dose_amt  0.027 ± NA                4         Dead
fit_death$covariate_effects()
#>   covariate        beta target_state_idx target_state
#> 1       age -0.026 ± NA                2         Dead
#> 2      CrCL  0.185 ± NA                2         Dead
#> 3    weight  0.029 ± NA                2         Dead
#> 4  dose_amt -0.270 ± NA                2         Dead
```

``` r
fit_ms_eh$plot_pk()
#> PK simulation
```

![](exposure-hazard_files/figure-html/unnamed-chunk-10-1.png)
