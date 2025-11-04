functions {

  // log of hazard multiplier
  matrix compute_log_hazard_multiplier(
      data int N_int,
      array[] vector beta_oth,
      array[] vector beta_auc,
      data array[] vector x_haz,
      array[] vector x_auc,
      real auc_loc,
      real auc_scale,
      data array[] int ttype
  ) {
    int N_trans = size(ttype);
    int nc_haz = size(x_haz);
    matrix[N_int, N_trans] log_C_haz = rep_matrix(0.0, N_int, N_trans);
    for(j in 1:N_trans){
      int h = ttype[j];
      if(nc_haz > 0){
        for(k in 1:nc_haz){
          log_C_haz[,j] += beta_oth[k][h] * x_haz[k];
        }
      }
      if(size(beta_auc)==1){
        log_C_haz[,j] += beta_auc[1][h] * ((x_auc[1] - auc_loc) / auc_scale);
      }
    }
    return(log_C_haz);
  }

  // Sigmoidal curve
  vector sigmoid(vector x, real steepness, real midpoint){
    return 1.0 ./ (1.0 + exp(-steepness*(x - midpoint)));
  }

  // Transform mu, sig, z to actual parameter
  vector raw_params_to_actual(array[] real mu, array[] real sigma,
      vector z, data array[] int tt){
    return(to_vector(mu[tt]) + to_vector(sigma[tt]) .* z);
  }

  // Compute log of baseline hazard
  vector log_basehaz(matrix SBF, vector w, real log_w0){
    return(SBF * w + log_w0);
  }

  // Compute log of instant hazard
  vector log_hazard(vector log_C_haz, matrix SBF, vector w, real log_w0){
    return(log_C_haz + log_basehaz(SBF, w, log_w0));
  }

  // Transform params
  matrix compute_theta_pk(
      array[] vector log_z, vector log_mu, vector log_sig,
      vector beta_ka, vector beta_CL, vector beta_V2,
      array[] vector x_ka, array[] vector x_CL, array[] vector x_V2
  ){
    int N_id = size(log_z);
    matrix[N_id, 3] log_theta;
    for(n in 1:N_id) {
      real v1 = -2 + log_mu[1] + log_z[n][1] * log_sig[1] +
        sum(beta_ka .* x_ka[n]);
      real v2 = -2 + log_mu[2] + log_z[n][2] * log_sig[2] +
        sum(beta_CL .* x_CL[n]);
      real v3 = -2 + log_mu[3] + log_z[n][3] * log_sig[3] +
        sum(beta_V2 .* x_V2[n]);

      // ka
      log_theta[n, 1] = v1;

      // CL
      log_theta[n, 2] = v2;

      // V2
      log_theta[n, 3] = v3;

    }
    return(exp(log_theta));
  }

  // Two-cpt population PK model with multiple doses (steady state)
  array[] vector pop_2cpt_ss(array[] vector t, data vector dose, matrix theta,
      data real tau) {

    int N_id = size(dose);
    int N_t = num_elements(t[1]);
    vector[N_id] tau_ss = rep_vector(tau, N_id);

    vector[N_id] ka = theta[:,1];
    vector[N_id] CL = theta[:,2];
    vector[N_id] V2 = theta[:,3];

    vector[N_id] ke = CL ./ V2;
    vector[N_id] A = (dose ./ V2) .* (ka ./ (ka-ke));

    array[N_id] vector[N_t] conc;
    vector[N_id] ma = A .* inv(-expm1(-ka.*tau_ss)); // 1/(1-exp(-ka*tau))
    vector[N_id] me = A .* inv(-expm1(-ke.*tau_ss));
    vector[N_t] tt;

    for(n in 1:N_id){
      tt = fmod(t[n], tau_ss[n]);
      conc[n] = me[n] * exp(-ke[n]*tt) - ma[n] * exp(-ka[n] * tt);
    }
    return(conc);
  }

  // Two-cpt population PK model (steady state peak time after dose)
  array[] vector pop_2cpt_ss_peak_time(matrix theta, data real tau){

    int N_id = rows(theta);
    vector[N_id] tau_ss = rep_vector(tau, N_id);

    vector[N_id] ka = theta[:,1];
    vector[N_id] CL = theta[:,2];
    vector[N_id] V2 = theta[:,3];

    vector[N_id] ke = CL ./ V2;

    array[N_id] vector[1] t_peak;
    vector[N_id] ma = inv(-expm1(-ka.*tau_ss)); // 1/(1-exp(-ka*tau))
    vector[N_id] me = inv(-expm1(-ke.*tau_ss));
    real r;
    for(n in 1:N_id){
      r = me[n]*ke[n]/(ma[n]*ka[n]);
      t_peak[n] = rep_vector(log(r)/(ke[n]-ka[n]), 1);
    }
    return(t_peak);
  }

  // Analytic solution with general initial condition (A0, C0)
  real two_cpt_central(real t, real ka, real CL, real V2, real A0, real C0){
    if(t < 0){
      reject("t = ", t,", should be non-negative");
    }
    real ke = CL / V2;
    return(A0 * ka/(ka-ke) * (exp(-ke*t) - exp(-ka*t)) + C0 * exp(-ke*t));
  }

  // Analytic solution with general initial condition A0
  real two_cpt_depot(real t, real ka, real A0){
    if(t < 0){
      reject("t = ", t,", should be non-negative");
    }
    return(A0*exp(-t*ka));
  }

  // Two-cpt population PK model with multiple doses (partly steady-state)
  // For each subject,
  // * dose_ss = dose amount in SS
  // * times = the time points after SS
  // * doses = the doses taken after SS
  array[,] vector pop_2cpt_partly_ss_stage1(
    data vector dose_ss,
    data array[] vector times,
    data array[] vector doses,
    matrix theta,
    data real tau
  ) {

    int N_id = num_elements(dose_ss);
    vector[N_id] tau_ss = rep_vector(tau, N_id);
    int N_last = num_elements(times[1]);

    vector[N_id] ka = theta[:,1];
    vector[N_id] CL = theta[:,2];
    vector[N_id] V2 = theta[:,3];
    real A0;
    real C0;
    real t0;

    // Loop over subjects
    array[2, N_id] vector[N_last] amt; // conc[1] = depot, conc[2] = central
    for(n in 1:N_id){

      // Find drug concentration at last SS trough time
      amt[1,n,1] = two_cpt_depot_ss0(tau_ss[n], dose_ss[n], ka[n]) - dose_ss[n];
      amt[2,n,1] = two_cpt_central_ss0(
        tau_ss[n], dose_ss[n], ka[n], CL[n], V2[n]
      );

      // Loop over remaining doses
      for(j in 2:N_last){
        A0 = amt[1,n,j-1];
        C0 = amt[2,n,j-1];
        t0 = times[n][j-1];
        amt[1,n,j] = two_cpt_depot(times[n][j] - t0, ka[n], A0 + doses[n][j-1]);
        amt[2,n,j] = two_cpt_central(
          times[n][j] - t0, ka[n], CL[n], V2[n], A0 + doses[n][j-1], C0
        );
      }
    }
    return(amt);
  }

  // How many elements in increasing vector x are <= y
  int rank_of(vector x, real y){
    int N = num_elements(x);
    int count = 0;
    for(n in 1:N){
      if(x[n] <= y){
        count = count + 1;
      }
    }
    return(count);
  }

  // Two-cpt population PK model with multiple doses (partly steady-state)
  // For each subject,
  // * t = vector of output time points
  // * dose_ss = dose amount in SS
  // * times = the time points after SS
  // * doses = the doses taken after SS
  // * amounts = drug amounts in both compartments
  // at dose occasions after SS (computed in stage1)
  array[] vector pop_2cpt_partly_ss_stage2(
    data array[] vector t,
    data vector dose_ss,
    data array[] vector times,
    data array[] vector doses,
    array[,] vector amounts,
    matrix theta,
    data real tau
  ) {

    int N_id = num_elements(dose_ss);
    vector[N_id] tau_ss = rep_vector(tau, N_id);
    int N_t = num_elements(t[1]);

    vector[N_id] ka = theta[:,1];
    vector[N_id] CL = theta[:,2];
    vector[N_id] V2 = theta[:,3];
    real A0;
    real C0;
    real t0;
    real t_eval;
    real t_ss_offset;
    int idx;

    // Loop over subjects
    array[N_id] vector[N_t] C; // central
    for(n in 1:N_id){

      // How much is last dose time off of a multiple of tau_ss?
      // i.e. which time of day is SS trough?
      t_ss_offset = fmod(times[n][1], tau_ss[n]);

      // Loop trough output times
      for(j in 1:N_t){
        t_eval = t[n][j];
        idx = rank_of(times[n], t_eval);
        if(idx == 0){
          C[n,j] = two_cpt_central_ss(
            t_eval - t_ss_offset, tau_ss[n], dose_ss[n], ka[n], CL[n], V2[n]
          );
        } else {
          t0 = times[n][idx];
          A0 = amounts[1,n,idx];
          C0 = amounts[2,n,idx];
          C[n,j] = two_cpt_central(
            t_eval - t0, ka[n], CL[n], V2[n], A0 + doses[n][idx], C0
          );
        }
      }
      C[n] = C[n]/V2[n];
    }
    return(C);
  }

  // Two-cpt PK model (steady state at given time t)
  real two_cpt_central_ss(real t, real tau, real dose, real ka, real CL, real V2)
  {
    real ke = CL / V2;
    real A = dose * (ka / (ka-ke));
    real tt = fmod(t, tau);
    real ma = A * inv(-expm1(-ka*tau)); // 1/(1-exp(-ka*tau))
    real me = A * inv(-expm1(-ke*tau));
    return(me * exp(-ke*tt) - ma * exp(-ka * tt));
  }

  // Two-cpt PK model (steady state at trough)
  real two_cpt_central_ss0(real tau, real dose, real ka, real CL, real V2) {
    real ke = CL / V2;
    real A = dose * (ka / (ka-ke));
    return(A * (inv(-expm1(-ke*tau)) - inv(-expm1(-ka*tau))));
  }

  // Two-cpt PK model (steady state at trough)
  real two_cpt_depot_ss0(real tau, real dose, real ka) {
    return(dose * inv(-expm1(-ka*tau)));
  }

  // Dummy function
  real STAN_dummy_function(real x){
    return(x + 1.0);
  }

}

data {

  int<lower=1> N_sub;       // number of subjects
  int<lower=1> N_int;       // total number of intervals
  array[N_int] int<lower=1,upper=N_sub> idx_sub; // which subject is interval for
  int<lower=1> N_grid;      // number of integration grid points
  real<lower=0> delta_grid; // grid step size

  // which value in t_grid does each t_start and t_end time correspond to
  array[N_int] int<lower=1,upper=N_grid> t_start_idx_m1;
  array[N_int] int<lower=2,upper=N_grid+1> t_end_idx;

  // how much we need to "correct" the integral for each interval
  vector<lower=0>[N_int] correction_multiplier;

  // Transitions / intervals data
  int<lower=1> N_trans; // number of possible transitions
  int<lower=1> N_trans_types;
  int<lower=1> D_trans; // max num of occurrences for a transition
  int<lower=1> D_risk; // max num of at-risk intervals for a transition

  // Total number of occurred transitions / at-risk intervals
  array[N_trans] int<lower=0> sum_trans;
  array[N_trans] int<lower=1> sum_risk;

  // Indices of intervals where a transition occurred / was at risk to occur
  array[N_trans, D_trans] int<lower=0,upper=N_int> which_trans;
  array[N_trans, D_risk] int<lower=0,upper=N_int> which_risk;

  // Transition type (for example can be same as the target state)
  array[N_trans] int<lower=1,upper=N_trans_types> ttype;

  // Mean h0
  vector[N_trans] mu_w0;

  // Covariates that affect hazard multiplier (excluding exposure)
  int<lower=0> nc_haz; // number of hazard covariates
  array[nc_haz] vector[N_sub] x_haz;

  // Which covariates to include in model
  int<lower=0,upper=1> I_auc;

  // Flags
  int<lower=0,upper=1> do_pk;
  int<lower=0,upper=1> do_haz;

  // Number of basis functions
  int<lower=1> N_sbf;

  // Basis functions evaluated at interval cut points (end points)
  matrix[N_int, N_sbf] SBF;

  // Basis functions evaluated at t_grid
  matrix[N_grid, N_sbf] SBF_grid;

  // Which likelihoods to include
  int<lower=0, upper=1> omit_lik_haz;
  int<lower=0, upper=1> omit_lik_pk;

  // PK model data
  vector<lower=0>[N_sub] pk_lloq; // lower limit of quantification
  array[N_sub] vector<lower=0>[2] conc_pk;
  array[N_sub] vector<lower=0>[2] t_obs_pk;
  vector<lower=0>[N_sub] dose_ss;
  real<lower=0> tau_ss;

  // Used while fitting
  array[N_sub] vector<lower=0>[2] last_two_times; // 1st = last SS trough time
  array[N_sub] vector<lower=0>[2] last_two_doses;

  // PK model covariates
  int<lower=0> nc_ka; // num of predictors for ka
  int<lower=0> nc_CL; // num of predictors for CL
  int<lower=0> nc_V2; // num of predictors for V2
  array[N_sub] vector[nc_ka] x_ka; // covs that affect ka
  array[N_sub] vector[nc_CL] x_CL; // covs that affect CL
  array[N_sub] vector[nc_V2] x_V2; // covs that affect V2

  // AUC normalization
  real<lower=0> auc_loc;
  real<lower=0> auc_scale;

}

transformed data {
  array[N_sub] vector[1] t0_ss;
  for(n in 1:N_sub){
    t0_ss[n] = rep_vector(0.0, 1);
  }
  // Set x corresponding to each interval
  array[nc_haz] vector[N_int] x_haz_long;
  if(nc_haz > 0){
    for(j in 1:nc_haz){
      for(n in 1:N_int){
        x_haz_long[j][n] = x_haz[j][idx_sub[n]];
      }
    }
  }
}

parameters {

  // Aux params ~ N(0,1)
  array[do_haz] matrix[N_sbf, N_trans] z_weights; // SBF weights
  array[do_haz] vector[N_trans] z_w0;

  // Means
  array[do_haz] matrix<lower=-3, upper=3>[N_sbf, N_trans_types] mu_weights;

  // Stds
  array[do_haz] matrix<lower=0, upper=3>[N_sbf, N_trans_types] sig_weights;
  array[do_haz] real<lower=0> sig_w0;

  // Covariate effects
  array[do_haz, nc_haz] vector[N_trans_types] beta_oth;
  array[do_haz, I_auc] vector[N_trans_types] beta_auc;

  // PK model
  array[do_pk, N_sub] vector[3] log_z_pk;
  array[do_pk] vector[3] log_mu_pk;
  array[do_pk] vector<lower=0>[3] log_sig_pk;
  array[do_pk] vector[nc_ka] beta_ka;
  array[do_pk] vector[nc_CL] beta_CL;
  array[do_pk] vector[nc_V2] beta_V2;
  array[do_pk] real<lower=0> sigma_pk;

}

transformed parameters {

  // Baseline hazard
  array[do_haz, N_trans] real log_w0;
  array[do_haz, N_trans] vector[N_sbf] weights;

  // Baseline hazard
  if(do_haz == 1) {
    for(j in 1:N_trans){
      weights[1,j] = mu_weights[1][:,ttype[j]] +
        sig_weights[1][:,ttype[j]] .* z_weights[1][:,j];
      log_w0[1,j] = mu_w0[j] + sig_w0[1] * z_w0[1][j];
    }
  }

  // PK quantities (in-sample)
  array[do_pk] matrix[N_sub, 3] theta_pk;
  array[do_pk, 2, N_sub] vector[2] last_two_amounts_pk;
  array[do_pk, N_sub] vector[2] conc_mu_pk;
  array[do_pk, N_sub] vector[1] t1_ss;
  array[do_pk, N_sub] vector[1] ss_trough;
  array[do_pk, N_sub] vector[1] ss_peak;
  array[do_pk] vector[N_sub] ss_auc; // auc for each subject
  array[do_pk] vector[N_int] x_auc_long; // auc for each interval

  if(do_pk == 1){

    // PK params
    theta_pk[1] = compute_theta_pk(
      log_z_pk[1], log_mu_pk[1], log_sig_pk[1],
      beta_ka[1], beta_CL[1], beta_V2[1], x_ka, x_CL, x_V2
    );

    // Find drug amounts in both compartments at last two dose times
    last_two_amounts_pk[1] = pop_2cpt_partly_ss_stage1(
      dose_ss, last_two_times, last_two_doses, theta_pk[1], tau_ss
    );

    // Drug concentration estimated at t_obs
    conc_mu_pk[1] = pop_2cpt_partly_ss_stage2(
      t_obs_pk, dose_ss, last_two_times, last_two_doses, last_two_amounts_pk[1],
      theta_pk[1], tau_ss
    );

    // Trough, peak, and auc at steady state
    t1_ss[1] = pop_2cpt_ss_peak_time(theta_pk[1], tau_ss);
    ss_trough[1] = pop_2cpt_ss(t0_ss, dose_ss, theta_pk[1], tau_ss);
    ss_peak[1] = pop_2cpt_ss(t1_ss[1], dose_ss, theta_pk[1], tau_ss);
    ss_auc[1] = dose_ss ./ theta_pk[1][:,2]; // D/CL

    // Set AUC corresponding to each interval
    for(n in 1:N_int){
      x_auc_long[1][n] = ss_auc[1][idx_sub[n]];
    }
  }
}

model {

  // HAZARD PARAM PRIOR
  if(do_haz == 1){
    if(nc_haz > 0){
      for(k in 1:nc_haz){
        beta_oth[1, k] ~ normal(0, 1);
      }
    }
    if(I_auc==1){
      beta_auc[1, 1] ~ normal(0, 1); // only if has pk submodel
    }

    // Baseline hazard base level
    sig_w0[1] ~ normal(0, 3);
    to_vector(z_w0[1]) ~ std_normal();

    // Baseline hazard spline weights
    to_vector(mu_weights[1]) ~ normal(0, 1);
    to_vector(sig_weights[1]) ~ normal(0, 0.3);
    to_vector(z_weights[1]) ~ std_normal();
  }
  // END HAZARD PARAM PRIOR

  // PK PARAM PRIOR
  if(do_pk == 1){
    for(n in 1:size(log_z_pk[1])){
      log_z_pk[1, n] ~ normal(0, 1);
    }
    log_mu_pk[1] ~ normal(0, 3);
    log_sig_pk[1] ~ normal(0, 1);
    beta_ka[1] ~ normal(0, 1);
    beta_CL[1] ~ normal(0, 1);
    beta_V2[1] ~ normal(0, 1);
    sigma_pk[1] ~ normal(0, 0.5);
  }
  // END PK PARAM PRIOR

  // Hazard model likelihood
  if (omit_lik_haz == 0 && do_haz == 1) {
    for(h in 1:N_trans){

      // log of hazard multiplier on each interval
      matrix[N_int, N_trans] log_C_haz = compute_log_hazard_multiplier(
        N_int, beta_oth[1], beta_auc[1], x_haz_long, x_auc_long, auc_loc, auc_scale,
        ttype
      );

      // Ragged array access
      array[sum_risk[h]] int idx_atr = which_risk[h, 1:sum_risk[h]];
      array[sum_trans[h]] int idx_occ = which_trans[h, 1:sum_trans[h]];

      // Occurred transitions (log hazard at interval end time)
      target += log_hazard(
        log_C_haz[idx_occ, h], SBF[idx_occ,:], weights[1,h], log_w0[1,h]
      );

      // Evaluate baseline hazard at grid, prepad with zero
      vector[N_grid+1] h0_grid = rep_vector(0.0, N_grid+1);
      h0_grid[2:(N_grid+1)] = exp(log_basehaz(SBF_grid, weights[1,h], log_w0[1,h]));

      // Integrated baseline hazard at grid
      vector[N_grid+1] h0_int = cumulative_sum(h0_grid) * delta_grid;

      // Transitions that were at risk (- integrated hazard over interval)
      // Taking into account that indices are off by one because of
      // The prepadded zero
      target += - exp(log_C_haz[idx_atr, h]) .*
        (h0_int[t_end_idx[idx_atr]] - h0_int[t_start_idx_m1[idx_atr]]) .*
        correction_multiplier[idx_atr];
    }

  }

  // PK model likelihood
  if (omit_lik_pk == 0 && do_pk == 1) {
    for(n in 1:N_sub){
      for(k in 1:2){
        real val = conc_mu_pk[1][n][k] + 1e-8;
        if(conc_pk[n][k] >= pk_lloq[n]){
          target += lognormal_lpdf(conc_pk[n][k] | log(val), sigma_pk[1]);
        } else {
          target += lognormal_lcdf(pk_lloq[n] | log(val), sigma_pk[1]);
        }
      }
    }
  }
}
