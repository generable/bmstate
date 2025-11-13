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
      log_theta[n, 2] = -2 + log_mu[2] + log_z[n][2] * log_sig[2] +
        sum(beta_CL .* x_CL[n]);

      // V2
      log_theta[n, 3] = -2 + log_mu[3] + log_z[n][3] * log_sig[3] +
        sum(beta_V2 .* x_V2[n]);

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

  // Two-cpt PK model (partly steady-state)
  // For each subject,
  // * t = vector of output time points
  // * dose_ss = dose amount in SS
  // * times = the time points after SS
  // * doses = the doses taken after SS
  // * tau = dosing interval
  array[] vector pop_2cpt_partly_ss(
    data array[] vector t,
    data vector dose_ss,
    data array[] vector times,
    data array[] vector doses,
    matrix theta,
    data real tau
  ){

    // Find drug amounts in both compartments at last two dose times
    int N_sub = num_elements(dose_ss);
    array[2, N_sub] vector[2] amounts = pop_2cpt_partly_ss_stage1(
      dose_ss, times, doses, theta, tau
    );

    // Drug concentration estimated at t
    return(pop_2cpt_partly_ss_stage2(
      t, dose_ss, times, doses, amounts, theta, tau
    ));
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
