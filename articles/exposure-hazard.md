# Exposure-hazard multistate modeling with bmstate

This vignette is work in progress.

## Data simulation

### Setup

``` r
library(bmstate)
#> Attached bmstate 0.3.0. Type ?bmstate to get started.
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

  beta_xpsr <- similar_fit$get_draws("beta_xpsr")
  beta_xpsr[1, , ] <- t(beta_true$matrix[, 2])

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
    beta_xpsr = beta_xpsr,
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
#> Generating 600 paths
covs_dh <- unique(c(mod_true$data_covs(), "dose_amt"))
simdat_dh <- simdat$paths$subset_covariates(covs_dh, renamed_old = "dose", renamed_new = "dose_amt")
simdat_death <- as_single_event(simdat_dh, "Dead", null_state = "Healthy")

sa <- simdat$paths$subject_df$xpsr
mod_true$set_xpsr_normalizers(loc = mean(sa), scale = stats::sd(sa))
#> setting xpsr normalizers to loc = 7.34488, scale = 0.7649

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
#>  * Covariates = {age, CrCL, weight, dose, t_pre, t_post, conc_pre, conc_post, xpsr, ka, CL, V2, pk_lloq}
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
  hazard_covs = setdiff(mod_true$covs(), "xpsr"),
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
#> Shortest time interval (0.43180367938658) is smaller than delta_grid (1.09575). Consider increasing n_grid or decreasing t_max of the model.
#> Using stan file at /home/runner/work/_temp/Library/bmstate/stan/msm.stan
#> setting xpsr normalizers to loc = 7.40235, scale = 0.56995
#> setting max conc = 11103.35342
#> Initial log joint probability = -957545 
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#> Error evaluating model log probability: Non-finite gradient. 
#> Error evaluating model log probability: Non-finite gradient. 
#> Error evaluating model log probability: Non-finite gradient. 
#>       99       -8871.6      0.274837       777.639           1           1      138    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      100      -8861.85      0.492775       1504.19           1           1      139    
#> Optimization terminated normally:  
#>   Maximum number of iterations hit, may not be at an optima 
#> Finished in  1.0 seconds.
fit_true <- create_oracle_fit(fit_prelim, beta_true, h0_true)
fit_true$covariate_effects()
#>   covariate      beta target_state_idx target_state
#> 1       age  0.0 ± NA                2        Bleed
#> 2       age  0.0 ± NA                3       Stroke
#> 3       age  0.3 ± NA                4         Dead
#> 4      xpsr  1.0 ± NA                2        Bleed
#> 5      xpsr -1.0 ± NA                3       Stroke
#> 6      xpsr  0.0 ± NA                4         Dead
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
#> Shortest time interval (0.43180367938658) is smaller than delta_grid (1.09575). Consider increasing n_grid or decreasing t_max of the model.
#> Using stan file at /home/runner/work/_temp/Library/bmstate/stan/msm.stan
#> setting xpsr normalizers to loc = 7.40235, scale = 0.56995
#> setting max conc = 11103.35342
#> Initial log joint probability = -957395 
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>       99      -13466.9      0.134121       585.054      0.7989      0.7989      160    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      199      -10548.7      0.065471       205.978      0.6478      0.6478      267    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      299      -10287.3     0.0550827       101.984           1           1      374    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      399      -8833.24     0.0117592        498.88      0.5588      0.5588      499    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      499      -8773.69     0.0473253       200.512           1           1      605    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      599      -8762.98    0.00380779       23.4628      0.8102      0.8102      714    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      699      -8757.95     0.0330903       53.5451           1           1      824    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      799      -8755.28     0.0178851       44.4238           1           1      931    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      899       -8752.8    0.00678303       22.1161      0.6257      0.6257     1044    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      999      -8749.74     0.0455133       84.1514      0.4688      0.4688     1155    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1099       -8739.2    0.00686113       52.6303      0.4495      0.4495     1269    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1199      -8736.53    0.00194763       15.3804      0.2488           1     1385    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1299      -8735.71    0.00472893       13.9008           1           1     1494    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1399      -8735.42    0.00246577       24.7813      0.4241      0.4241     1602    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1499         -8735   0.000959264       12.4187           1           1     1711    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1599      -8734.41    0.00118824       41.6362      0.5342      0.5342     1819    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1699      -8724.82     0.0284222       23.9692           1           1     1933    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1799      -8720.56     0.0190098       58.0966           1           1     2042    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1899      -8719.85    0.00163301       14.6072           1           1     2151    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1999      -8719.61   0.000322125       8.37171      0.7877      0.7877     2259    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     2000      -8719.61   0.000604208       10.8696           1           1     2260    
#> Optimization terminated normally:  
#>   Maximum number of iterations hit, may not be at an optima 
#> Finished in  14.4 seconds.
fit_ms_dh <- fit_stan(mod_ms_dh, simdat_dh, method = "optimize", init = 0)
#> Shortest time interval (0.0999999999999943) is smaller than delta_grid (1.09575). Consider increasing n_grid or decreasing t_max of the model.
#> Using stan file at /home/runner/work/_temp/Library/bmstate/stan/msm.stan
#> Initial log joint probability = -721464 
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>       99      -6718.59     0.0843703       96.3895           1           1      120    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      199      -6650.16      0.285527       79.8082           1           1      225    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      299      -6642.57    0.00388817       20.9987      0.6704      0.6704      331    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      399      -6638.63     0.0163856       7.03785           1           1      440    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      499      -6632.07     0.0690188       34.2961           1           1      543    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      599      -6626.16     0.0067006       8.83777      0.9939      0.9939      652    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      699      -6625.32     0.0761852       14.0505           1           1      757    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      799      -6622.02     0.0285054       17.4926       0.286           1      867    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      899      -6614.12    0.00792563       11.0484       0.853       0.853      981    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      999      -6612.06     0.0328372       23.4549           1           1     1089    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1099      -6609.61     0.0766424       10.7358           1           1     1197    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1199      -6606.72     0.0111774       12.9944           1           1     1314    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1299      -6605.14     0.0208145       11.3211           1           1     1430    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1399      -6603.03     0.0191563       31.4591      0.3404      0.3404     1545    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1499      -6601.24     0.0120439       6.81678           1           1     1652    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1599      -6600.98    0.00482129       2.18416       1.656      0.1656     1761    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1699      -6600.87    0.00177176       3.67063           1           1     1866    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1799      -6600.83    0.00287955       1.91997      0.4346      0.4346     1976    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1899       -6600.8    0.00308483      0.885636           1           1     2086    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1999      -6600.78    0.00154565       1.45032           1           1     2192    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     2000      -6600.78   0.000849788      0.513626      0.9207      0.9207     2193    
#> Optimization terminated normally:  
#>   Maximum number of iterations hit, may not be at an optima 
#> Finished in  10.1 seconds.
fit_death <- fit_stan(mod_death, simdat_death, method = "optimize", init = 0)
#> Using stan file at /home/runner/work/_temp/Library/bmstate/stan/msm.stan
#> Initial log joint probability = -255418 
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>       99      -3977.35     0.0183187      0.743335           1           1      121    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      132       -3977.3   0.000260916     0.0623965      0.9149      0.9149      155    
#> Optimization terminated normally:  
#>   Convergence detected: relative gradient magnitude is below tolerance 
#> Finished in  0.1 seconds.
fit_ms_eh$plot_h0()
```

![](exposure-hazard_files/figure-html/unnamed-chunk-9-1.png)

``` r
fit_ms_eh$covariate_effects()
#>   covariate       beta target_state_idx target_state
#> 1       age  0.42 ± NA                2        Bleed
#> 2       age  0.10 ± NA                3       Stroke
#> 3       age  0.31 ± NA                4         Dead
#> 4      xpsr -0.47 ± NA                2        Bleed
#> 5      xpsr -0.80 ± NA                3       Stroke
#> 6      xpsr -0.24 ± NA                4         Dead
fit_ms_dh$covariate_effects()
#>    covariate        beta target_state_idx target_state
#> 1        age  0.582 ± NA                2        Bleed
#> 2        age -0.600 ± NA                3       Stroke
#> 3        age  0.221 ± NA                4         Dead
#> 4       CrCL -0.485 ± NA                2        Bleed
#> 5       CrCL  0.688 ± NA                3       Stroke
#> 6       CrCL  0.054 ± NA                4         Dead
#> 7     weight -0.471 ± NA                2        Bleed
#> 8     weight  0.576 ± NA                3       Stroke
#> 9     weight -0.034 ± NA                4         Dead
#> 10  dose_amt  0.858 ± NA                2        Bleed
#> 11  dose_amt -1.060 ± NA                3       Stroke
#> 12  dose_amt -0.067 ± NA                4         Dead
fit_death$covariate_effects()
#>   covariate        beta target_state_idx target_state
#> 1       age -0.071 ± NA                2         Dead
#> 2      CrCL  0.232 ± NA                2         Dead
#> 3    weight  0.080 ± NA                2         Dead
#> 4  dose_amt -0.179 ± NA                2         Dead
```

``` r
fit_ms_eh$plot_pk()
#> PK simulation
```

![](exposure-hazard_files/figure-html/unnamed-chunk-10-1.png)
