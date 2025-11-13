
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
