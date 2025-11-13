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
  // * dose_ss = dose amount
  // * dose_times = the time points after SS
  // * doses = the doses taken after SS
  // * tau = dosing interval
  array[] vector pop_2cpt_partly_ss(
    data array[] vector t,
    data vector dose_ss,
    data array[] vector dose_times,
    data array[] vector doses,
    matrix theta,
    data real tau
  ){

    // Find drug amounts in both compartments at dose_times
    int N_t = num_elements(t[1]);
    int N_sub = num_elements(dose_ss);
    int D = num_elements(dose_times[1]);
    array[2, N_sub] vector[D] amounts = pop_2cpt_partly_ss_stage1(
      dose_ss, dose_times, doses, theta, tau
    );

    array[N_sub] vector[N_t] conc = pop_2cpt_partly_ss_stage2(
      t, dose_ss, dose_times, doses, amounts, theta, tau
    );

    // Drug concentration at t
    return(conc);
  }
