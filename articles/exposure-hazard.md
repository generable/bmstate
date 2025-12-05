# Exposure-hazard multistate modeling with bmstate

This vignette is work in progress.

## Data simulation

### Setup

``` r
library(bmstate)
#> Attached bmstate 0.3.1. Type ?bmstate to get started.
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
#>       99      -8890.48     0.0748033       439.968           1           1      144    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      100      -8884.79      0.164095       1041.58           1           1      145    
#> Optimization terminated normally:  
#>   Maximum number of iterations hit, may not be at an optima 
#> Finished in  0.7 seconds.
#> If these don't roughly match, consider refitting after setting xpsr normalizer loc and scale closer to estimated mean and sd. Otherwise interpret baseline hazards and xpsr effect size accodingly.
#>  - xpsr normalization loc = 7.40235, mean estimated xpsr = 12.36832
#>  - xpsr normalization scale = 0.56995, estimated xpsr sd = 0.70053
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
#>       99      -13219.4       3.70974       1123.64           1           1      160    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      199      -9592.67     0.0380356       206.009           1           1      270    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      299      -9078.43    0.00983121       138.463           1           1      377    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      399      -8936.79      0.117285       143.826           1           1      485    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      499      -8856.68     0.0423536       212.376           1           1      593    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      599      -8802.65     0.0445093       124.903           1           1      703    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      699      -8788.46     0.0245598       87.2478           1           1      818    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      799      -8743.91      0.029015       70.7765           1           1      925    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      899      -8736.77     0.0108436       56.4198       0.367           1     1036    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      999      -8733.21    0.00274683       30.6018           1           1     1146    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1099      -8731.11    0.00307861       22.8986           1           1     1260    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1199      -8726.25     0.0217849       176.386      0.4243           1     1368    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1299      -8721.56     0.0143369       17.6503           1           1     1474    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1399      -8720.08    0.00193987       7.60364      0.8701      0.8701     1580    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1499      -8719.75    0.00569246       25.9519      0.5364      0.5364     1688    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1599      -8719.53    0.00151107       12.8847           1           1     1798    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1699      -8719.41    0.00671789       10.1732           1           1     1909    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1799      -8719.17    0.00245639       3.90602           1           1     2019    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1899      -8719.07   0.000581305       4.54652           1           1     2126    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1999      -8718.97    0.00396648       7.05303      0.8273      0.8273     2236    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     2000      -8718.96    0.00686852       4.71608           1           1     2237    
#> Optimization terminated normally:  
#>   Maximum number of iterations hit, may not be at an optima 
#> Finished in  10.7 seconds.
#> If these don't roughly match, consider refitting after setting xpsr normalizer loc and scale closer to estimated mean and sd. Otherwise interpret baseline hazards and xpsr effect size accodingly.
#>  - xpsr normalization loc = 7.40235, mean estimated xpsr = 12.45435
#>  - xpsr normalization scale = 0.56995, estimated xpsr sd = 0.73747
fit_ms_dh <- fit_stan(mod_ms_dh, simdat_dh, method = "optimize", init = 0)
#> Shortest time interval (0.0999999999999943) is smaller than delta_grid (1.09575). Consider increasing n_grid or decreasing t_max of the model.
#> Using stan file at /home/runner/work/_temp/Library/bmstate/stan/msm.stan
#> Initial log joint probability = -721464 
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>       99      -6710.49     0.0962278       69.9461           1           1      116    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      199       -6645.5     0.0270824       28.7697      0.7753      0.7753      229    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      299      -6640.61      0.022021       19.2347      0.3728      0.3728      340    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      399      -6617.55     0.0537164       25.6929           1           1      446    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      499      -6612.47     0.0248891       30.4679      0.4849           1      553    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      599      -6609.23     0.0240955       18.5197           1           1      664    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      699      -6607.77    0.00351639       10.6477      0.6013      0.6013      770    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      799      -6607.62    0.00860762       3.16114           1           1      873    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      899      -6607.33    0.00500078       4.71301           1           1      982    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      999       -6606.9     0.0163521       5.11164           1           1     1089    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1099      -6603.21     0.0274515       13.5618      0.9937      0.9937     1202    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1199      -6601.94    0.00884304       9.71218      0.8931      0.8931     1311    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1299      -6600.65    0.00188555       4.85025           1           1     1417    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1399      -6599.65    0.00511112       5.37397      0.9658      0.9658     1519    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1499      -6599.27    0.00432644       3.03535           1           1     1629    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1599      -6599.22     0.0013081       1.30039           1           1     1737    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1699      -6599.18   0.000627939       1.95542      0.4945           1     1847    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1799      -6599.09    0.00307922       3.68924      0.2852           1     1957    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1899      -6599.02     0.0179589       1.62049           1           1     2067    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     1999      -6598.97   0.000537596       1.87777      0.2938           1     2181    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>     2000      -6598.97   0.000303309        1.5957      0.2968      0.2968     2182    
#> Optimization terminated normally:  
#>   Maximum number of iterations hit, may not be at an optima 
#> Finished in  6.3 seconds.
fit_death <- fit_stan(mod_death, simdat_death, method = "optimize", init = 0)
#> Using stan file at /home/runner/work/_temp/Library/bmstate/stan/msm.stan
#> Initial log joint probability = -255418 
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>       99      -3738.33     0.0427164       2.19728           1           1      128    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      156      -3738.27    0.00202098      0.080757       1.685      0.1685      198    
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
#> 3       age  0.35 ± NA                4         Dead
#> 4      xpsr -0.48 ± NA                2        Bleed
#> 5      xpsr -0.80 ± NA                3       Stroke
#> 6      xpsr -0.31 ± NA                4         Dead
fit_ms_dh$covariate_effects()
#>    covariate        beta target_state_idx target_state
#> 1        age  0.587 ± NA                2        Bleed
#> 2        age -0.603 ± NA                3       Stroke
#> 3        age  0.222 ± NA                4         Dead
#> 4       CrCL -0.489 ± NA                2        Bleed
#> 5       CrCL  0.694 ± NA                3       Stroke
#> 6       CrCL  0.055 ± NA                4         Dead
#> 7     weight -0.475 ± NA                2        Bleed
#> 8     weight  0.580 ± NA                3       Stroke
#> 9     weight -0.035 ± NA                4         Dead
#> 10  dose_amt  0.863 ± NA                2        Bleed
#> 11  dose_amt -1.070 ± NA                3       Stroke
#> 12  dose_amt -0.067 ± NA                4         Dead
fit_death$covariate_effects()
#>   covariate        beta target_state_idx target_state
#> 1       age -0.079 ± NA                2         Dead
#> 2      CrCL  0.196 ± NA                2         Dead
#> 3    weight  0.094 ± NA                2         Dead
#> 4  dose_amt -0.166 ± NA                2         Dead
```

``` r
fit_ms_eh$plot_pk()
#> PK simulation
```

![](exposure-hazard_files/figure-html/unnamed-chunk-10-1.png)
