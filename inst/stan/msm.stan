functions {

  // log of hazard multiplier
  matrix compute_log_hazard_multiplier(
      array[] matrix beta_sub, 
      array[] matrix beta_cnt, 
      array[] vector beta_oth, 
      array[] vector beta_avg, 
      data array[] int x_sub, 
      data array[] int x_cnt,
      data array[] vector x_oth, 
      array[] vector ss_auc, 
      array[] int ttype
  ) {
    int N = size(x_sub);
    int N_trans = size(ttype);
    int N_oth = size(x_oth);
    matrix[N, N_trans] log_C_haz = rep_matrix(0.0, N, N_trans);
    for(j in 1:N_trans){
      if(size(beta_sub)==1){
        log_C_haz[,j] += to_vector(beta_sub[1][ttype[j], x_sub]);
      }
      if(size(beta_cnt)==1){
        log_C_haz[,j] += to_vector(beta_cnt[1][ttype[j], x_cnt]);
      }
      if(N_oth > 0){
        for(k in 1:N_oth){
          log_C_haz[,j] += beta_oth[k][ttype[j]] * x_oth[k];
        }
      }
      if(size(beta_avg)==1){
        log_C_haz[,j] += beta_avg[1][ttype[j]] * (ss_auc[1][x_sub]/24.0-0.08)/0.04;
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
      
      // ka
      log_theta[n, 1] = -2 + log_mu[1] + log_z[n][1] * log_sig[1] + 
        sum(beta_ka .* x_ka[n]);
        
      // CL
      log_theta[n, 2] = 1 + log_mu[2] + log_z[n][2] * log_sig[2] +
        sum(beta_CL .* x_CL[n]);
      
      // V2
      log_theta[n, 3] = -2 + log_mu[3] + log_z[n][3] * log_sig[3] +
        sum(beta_V2 .* x_V2[n]);
      
    }
    return(exp(log_theta));
  }

  // Two-cpt population PK model with multiple doses (steady state)
  array[] vector pop_2cpt_ss(array[] vector t, data vector dose, matrix theta) {
      
    int N_id = size(dose);
    int N_t = num_elements(t[1]);
    vector[N_id] tau = rep_vector(24.0, N_id);
        
    vector[N_id] ka = theta[:,1];
    vector[N_id] CL = theta[:,2];
    vector[N_id] V2 = theta[:,3];
        
    vector[N_id] ke = CL ./ V2;
    vector[N_id] A = (dose ./ V2) .* (ka ./ (ka-ke));
        
    array[N_id] vector[N_t] conc;
    vector[N_id] ma = A .* inv(-expm1(-ka.*tau)); // 1/(1-exp(-ka*tau))
    vector[N_id] me = A .* inv(-expm1(-ke.*tau));
    vector[N_t] tt;
        
    for(n in 1:N_id){
      tt = fmod(t[n], tau[n]);
      conc[n] = me[n] * exp(-ke[n]*tt) - ma[n] * exp(-ka[n] * tt);
    }
    return(conc);
  }
  
  // Two-cpt population PK model (steady state peak time after dose)
  array[] vector pop_2cpt_ss_peak_time(matrix theta){
    
    int N_id = rows(theta);
    vector[N_id] tau = rep_vector(24.0, N_id);
        
    vector[N_id] ka = theta[:,1];
    vector[N_id] CL = theta[:,2];
    vector[N_id] V2 = theta[:,3];
        
    vector[N_id] ke = CL ./ V2;
        
    array[N_id] vector[1] t_peak;
    vector[N_id] ma = inv(-expm1(-ka.*tau)); // 1/(1-exp(-ka*tau))
    vector[N_id] me = inv(-expm1(-ke.*tau));
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
  // * t_last_ss = last time point when subject is assumed to be in SS trough
  // * times = the time points after SS
  // * doses = the doses taken after SS
  array[,] vector pop_2cpt_partly_ss_stage1(
    data vector dose_ss, 
    data array[] vector times, 
    data array[] vector doses, 
    matrix theta
  ) {
      
    int N_id = num_elements(dose_ss);
    vector[N_id] tau_ss = rep_vector(24.0, N_id);
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
        amt[1,n,j] = two_cpt_depot(times[n][j] - t0, ka[n], A0 + doses[n][j]);
        amt[2,n,j] = two_cpt_central(
          times[n][j] - t0, ka[n], CL[n], V2[n], A0 + doses[n][j], C0
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
    matrix theta
  ) {
      
    int N_id = num_elements(dose_ss);
    vector[N_id] tau_ss = rep_vector(24.0, N_id);
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
  
}

data {

  int<lower=1> N_obs;       // total number of intervals
  int<lower=1> N_grid;      // number of integration grid points
  real<lower=0> delta_grid; // grid step size  
  vector<lower=0>[N_obs] correction_multiplier;
  
  // which value in t_grid does each t_start and t_end time correspond to
  array[N_obs] int<lower=1,upper=N_grid> t_start_idx_m1;  
  array[N_obs] int<lower=2,upper=N_grid+1> t_end_idx;
  
  // Transitions / intervals data
  int<lower=1> N_trans; // number  of possible transitions
  int<lower=1> N_trans_types;
  int<lower=1> D_trans; // max num of occurrences for a transition
  int<lower=1> D_risk; // max num of at-risk intervals for a transition
  
  // Total number of occurred transitions / at-risk intervals
  array[N_trans] int<lower=0> sum_trans;
  array[N_trans] int<lower=1> sum_risk;
  
  // Indices of intervals where a transition occurred / was at risk to occur
  array[N_trans, D_trans] int<lower=0,upper=N_obs> which_trans;
  array[N_trans, D_risk] int<lower=0,upper=N_obs> which_risk;

  // Transition type (for example can be same as the target state)
  array[N_trans] int<lower=1,upper=N_trans_types> ttype;
  
  // Mean h0
  vector[N_trans] mu_w0;
  
  // Subject covariates that affect hazard multiplier
  int<lower=1> N_sub;
  int<lower=0> N_country;
  array[N_obs] int<lower=1,upper=N_sub> x_sub;
  array[N_obs] int<lower=1,upper=N_country> x_cnt;
  int<lower=0> N_oth; // number of other covariates
  array[N_oth] vector[N_obs] x_oth;
  
  // Subject covariates that affect hazard multiplier (out-of-sample)
  int<lower=1> N_sub_oos;
  array[N_sub_oos] int<lower=1,upper=N_sub_oos> x_sub_oos;
  array[N_sub_oos] int<lower=1,upper=N_country> x_cnt_oos;
  array[N_oth] vector[N_sub_oos] x_oth_oos;
  
  // Which covariates to include in model
  int<lower=0,upper=1> I_sub;
  int<lower=0,upper=1> I_cnt;
  int<lower=0,upper=1> I_avg;
  
  // Skips PK evaluation completely if this is 0
  int<lower=0,upper=1> do_pk;
  
  // Number of basis functions
  int<lower=1> N_sbf;
  
  // Basis functions evaluated at interval cut points (end points)
  matrix[N_obs, N_sbf] SBF;
  
  // Basis functions evaluated at t_grid
  matrix[N_grid, N_sbf] SBF_grid;
  
  // Basis functions evaluated at predict times
  int<lower=0> N_pred;
  matrix[N_pred, N_sbf] SBF_pred;
  
  // Which likelihoods to include
  int<lower=0, upper=1> omit_lik_hazard;
  int<lower=0, upper=1> omit_lik_pk;
  
  // PK model data
  int<lower=2> N_last;
  vector<lower=0>[N_sub] pk_lloq; // lower limit of quantification
  array[N_sub] vector<lower=0>[2] conc_pk;
  
  array[N_sub] vector<lower=0>[2] t_obs_pk;
  array[N_sub_oos] vector<lower=0>[2] t_obs_pk_oos;
  vector<lower=0>[N_sub] dose_ss;
  vector<lower=0>[N_sub_oos] dose_ss_oos;
  
  // Used while checking
  array[N_sub] vector<lower=0>[N_last] last_times; // 1st = last SS trough time
  array[N_sub] vector<lower=0>[N_last] last_doses;
  array[N_sub_oos] vector<lower=0>[N_last] last_times_oos;
  array[N_sub_oos] vector<lower=0>[N_last] last_doses_oos; 
  
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
  array[N_sub_oos] vector[nc_ka] x_ka_oos; // covs that affect ka
  array[N_sub_oos] vector[nc_CL] x_CL_oos; // covs that affect CL
  array[N_sub_oos] vector[nc_V2] x_V2_oos; // covs that affect V2
  
  // Concentration prediction in a dense grid
  int<lower=0> N_pred_pk;
  int<lower=0,upper=1> pred_pk;
  array[N_sub] vector<lower=0>[N_pred_pk] t_pred_pk;
  array[N_sub_oos] vector<lower=0>[N_pred_pk] t_pred_pk_oos;
  
  // Options for beta_computing beta_sub in GQ
  // 0: set z to 0
  // 1: draw z from prior
  // 2: use actual z fitted for the subject (not possible for out-of-sample)
  int<lower=0, upper=2> z_mode_is;
  int<lower=0, upper=1> z_mode_oos;
}

transformed data {
  array[N_sub] vector[1] t0_ss;
  array[N_sub_oos] vector[1] t0_ss_oos;
  array[N_sub_oos] vector[3] z0;
  for(n in 1:N_sub){
    t0_ss[n] = rep_vector(0.0, 1);
  }
  for(n in 1:N_sub_oos){
    t0_ss_oos[n] = rep_vector(0.0, 1);
    z0[n] = rep_vector(0.0, 3);
  }
}
parameters {
  
  // Aux params ~ N(0,1)
  matrix[N_sbf, N_trans] z_weights; // SBF weights
  vector[N_trans] z_w0;
  
  // Means
  matrix<lower=-2, upper=2>[N_sbf, N_trans_types] mu_weights;
  
  // Stds
  matrix<lower=0, upper=2>[N_sbf, N_trans_types] sig_weights;
  real<lower=0> sig_w0;
  
  // Covariate effects
  array[I_sub] matrix[N_trans_types, N_sub] beta_sub;
  array[I_cnt] matrix[N_trans_types, N_country] beta_cnt;
  array[N_oth] vector[N_trans_types] beta_oth;
  array[I_avg] vector[N_trans_types] beta_avg;
  
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
  array[N_trans] real log_w0;
  array[N_trans] vector[N_sbf] weights;
  
  // Baseline hazard
  for(j in 1:N_trans){
    weights[j,] = mu_weights[:,ttype[j]] + 
      sig_weights[:,ttype[j]] .* z_weights[:,j];

    log_w0[j] = mu_w0[j] + sig_w0 * z_w0[j];
  }
  
  // PK quantities (in-sample)
  array[do_pk] matrix[N_sub, 3] theta_pk;
  array[do_pk, 2, N_sub] vector[2] last_two_amounts_pk;
  array[do_pk, N_sub] vector[2] conc_mu_pk;
  array[do_pk, N_sub] vector[1] t1_ss;
  array[do_pk, N_sub] vector[1] ss_trough;
  array[do_pk, N_sub] vector[1] ss_peak;
  array[do_pk] vector[N_sub] ss_auc;
    
  if(do_pk == 1){
    
    // PK params
    theta_pk[1] = compute_theta_pk(
      log_z_pk[1], log_mu_pk[1], log_sig_pk[1], 
      beta_ka[1], beta_CL[1], beta_V2[1], x_ka, x_CL, x_V2
    );
    
    // Find drug amounts in both compartments at last two dose times
    last_two_amounts_pk[1] = pop_2cpt_partly_ss_stage1(
      dose_ss, last_two_times, last_two_doses, theta_pk[1]
    );
  
    // Drug concentration estimated at t_obs
    conc_mu_pk[1] = pop_2cpt_partly_ss_stage2(
      t_obs_pk, dose_ss, last_two_times, last_two_doses, last_two_amounts_pk[1], 
      theta_pk[1]
    );
  
    // Trough, peak, and auc at steady state
    t1_ss[1] = pop_2cpt_ss_peak_time(theta_pk[1]);
    ss_trough[1] = pop_2cpt_ss(t0_ss, dose_ss, theta_pk[1]);
    ss_peak[1] = pop_2cpt_ss(t1_ss[1], dose_ss, theta_pk[1]);
    ss_auc[1] = dose_ss ./ theta_pk[1][:,2]; // D/CL
  }
  
  // log of hazard multiplier on each interval
  matrix[N_obs, N_trans] log_C_haz = compute_log_hazard_multiplier(
    beta_sub, beta_cnt, beta_oth, beta_avg, x_sub, x_cnt, x_oth, ss_auc, ttype
  );
}

model {
  
  // Prior
  if(I_sub==1){
    to_vector(beta_sub[1]) ~ normal(0, 1);
  }
  if(I_cnt==1){
    to_vector(beta_cnt[1]) ~ normal(0, 1);
  }
  if(N_oth > 0){
    for(k in 1:N_oth){
      beta_oth[k] ~ normal(0, 1);
    }
  }
  if(I_avg==1){
    beta_avg[1] ~ normal(0, 1);
  }
  
  // PK model
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
  
  // Weights
  sig_w0 ~ normal(0, 5);
  to_vector(z_w0) ~ std_normal();
  to_vector(mu_weights) ~ normal(0, 1);
  to_vector(sig_weights) ~ normal(0, 0.3);
  to_vector(z_weights) ~ std_normal();

  // Hazard model likelihood
  if (omit_lik_hazard == 0) {
    for(h in 1:N_trans){

      // Ragged array access
      array[sum_risk[h]] int idx_atr = which_risk[h, 1:sum_risk[h]];
      array[sum_trans[h]] int idx_occ = which_trans[h, 1:sum_trans[h]];
      
      // Occurred transitions (log hazard at interval end time)
      target += log_hazard(
        log_C_haz[idx_occ, h], SBF[idx_occ,:], weights[h], log_w0[h]
      );
        
      // Evaluate baseline hazard at grid, prepad with zero
      vector[N_grid+1] h0_grid = rep_vector(0.0, N_grid+1);
      h0_grid[2:(N_grid+1)] = exp(log_basehaz(SBF_grid, weights[h], log_w0[h]));
        
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
        real val = conc_mu_pk[1][n][k] + 1e-7;
        if(conc_pk[n][k] >= pk_lloq[n]){
          target += lognormal_lpdf(conc_pk[n][k] | log(val), sigma_pk[1]);
        } else {
          target += lognormal_lcdf(pk_lloq[n] | log(val), sigma_pk[1]);
        }
      }
    }
  }
}

generated quantities {
  
  array[N_trans] vector[N_pred] log_h0;
  
  // Log baseline hazard for each transition
  if(N_pred > 0){
    for(h in 1:N_trans){
      log_h0[h] = log_basehaz(SBF_pred, weights[h], log_w0[h]); // transition h
    }
  }
  
  // PK params (OOS subjects)
  array[do_pk, N_sub_oos] vector[3] log_z_pk_oos;
  if(do_pk == 1){
    for(n in 1:N_sub_oos){
      log_z_pk_oos[1, n] = rep_vector(0.0, 3);
      for(d in 1:3){
        log_z_pk_oos[1, n][d] = std_normal_rng();
      }
    }
  }

  
  // PK quantities (out-of-sample)
  array[do_pk] matrix[N_sub_oos, 3] theta_pk_oos;
  array[do_pk, 2, N_sub_oos] vector[N_last] amounts_pk_oos;
  array[do_pk, N_sub_oos] vector[2] conc_mu_pk_oos;
  array[do_pk, N_sub_oos] vector[1] t1_ss_oos;
  array[do_pk, N_sub_oos] vector[1] ss_trough_oos;
  array[do_pk, N_sub_oos] vector[1] ss_peak_oos;
  array[do_pk] vector[N_sub_oos] ss_auc_oos;
  array[do_pk, N_sub] vector[2] conc_mu_pk_obs;
  
  // Predict
  array[do_pk, 2, N_sub] vector[N_last] amounts_pk;
  array[pred_pk*do_pk*N_sub] vector[N_pred_pk] conc_mu_pk_pred;
  array[pred_pk*do_pk*N_sub] vector[N_pred_pk] conc_mu_pk_pred_lasttwo;
  array[pred_pk*do_pk*N_sub_oos] vector[N_pred_pk] conc_mu_pk_pred_oos;
    
  if(do_pk == 1){
    
    // PK params
    theta_pk_oos[1] = compute_theta_pk(
      log_z_pk_oos[1], log_mu_pk[1], log_sig_pk[1], 
      beta_ka[1], beta_CL[1], beta_V2[1], x_ka_oos, x_CL_oos, x_V2_oos
    );
    
    // Find drug amounts in both compartments at last_n dose times
    amounts_pk[1] = pop_2cpt_partly_ss_stage1(
      dose_ss, last_times, last_doses, theta_pk[1]
    );
    amounts_pk_oos[1] = pop_2cpt_partly_ss_stage1(
      dose_ss_oos, last_times_oos, last_doses_oos, theta_pk_oos[1]
    );
  
    // Copied so can be extracted from GQ
    conc_mu_pk_obs[1] = conc_mu_pk[1];
  
    // Drug concentration estimated at t_obs (OOS)
    conc_mu_pk_oos[1] = pop_2cpt_partly_ss_stage2(
      t_obs_pk_oos, dose_ss_oos, last_times_oos, last_doses_oos, 
      amounts_pk_oos[1], theta_pk_oos[1]
    );
  
    // Trough, peak, and auc at steady state
    t1_ss_oos[1] = pop_2cpt_ss_peak_time(theta_pk_oos[1]);
    ss_trough_oos[1] = pop_2cpt_ss(t0_ss_oos, dose_ss_oos, theta_pk_oos[1]);
    ss_peak_oos[1]   = pop_2cpt_ss(t1_ss_oos[1], dose_ss_oos, theta_pk_oos[1]);
    ss_auc_oos[1] = dose_ss_oos ./ theta_pk_oos[1][:,2]; // D/CL
  }
  
  // Drug concentration estimated at t_pred_pk
  if(pred_pk == 1 && do_pk == 1){
    conc_mu_pk_pred = pop_2cpt_partly_ss_stage2(
      t_pred_pk, dose_ss, last_times, last_doses, amounts_pk[1], theta_pk[1]
    );
    conc_mu_pk_pred_lasttwo = pop_2cpt_partly_ss_stage2(
      t_pred_pk, dose_ss, last_two_times, last_two_doses, 
      last_two_amounts_pk[1], theta_pk[1]
    );
    conc_mu_pk_pred_oos = pop_2cpt_partly_ss_stage2(
      t_pred_pk_oos, dose_ss_oos, last_times_oos, last_doses_oos, 
      amounts_pk_oos[1], theta_pk_oos[1]
    );
  }
  
  // Frailty (IS and OOS subjects)
  array[I_sub] matrix[N_trans, N_sub_oos] beta_sub_oos;
  array[I_sub] matrix[N_trans, N_sub] beta_sub_is;
  if(I_sub==1){
    beta_sub_oos[1] = rep_matrix(0.0, N_trans, N_sub_oos);
    beta_sub_is[1] = rep_matrix(0.0, N_trans, N_sub);
    real zzz;
    for(j in 1:N_trans){
      for(k in 1:N_sub){
        if(z_mode_is == 0){
          beta_sub_is[1][j, k] = 0.0;
        } else if (z_mode_is == 1){
          beta_sub_is[1][j, k] = std_normal_rng();
        } else {
          beta_sub_is[1][j, k] = beta_sub[1][j, k];
        }
      }
      for(k in 1:N_sub_oos){
        if(z_mode_oos == 0){
          beta_sub_oos[1][j, k] = 0.0;
        } else {
          beta_sub_oos[1][j, k] = std_normal_rng();
        }
      }
    }
  }
  
  // log of hazard multiplier for IS subjects (possibly different than in
  // transformed params)
  matrix[N_obs, N_trans] log_C_haz_is = compute_log_hazard_multiplier(
    beta_sub_is, beta_cnt, beta_oth, beta_avg, x_sub, x_cnt, x_oth, ss_auc, ttype
  );
  
  // log of hazard multiplier for OOS subjects
  matrix[N_sub_oos, N_trans] log_C_haz_oos = compute_log_hazard_multiplier(
    beta_sub_oos, beta_cnt, beta_oth, beta_avg, x_sub_oos, x_cnt_oos, 
    x_oth_oos, ss_auc_oos, ttype
  );
  
}
