# Package index

## All functions

- [`DosingData`](https://generable.github.io/bmstate/reference/DosingData.md)
  : Dosing data class (R6 class)
- [`JointData`](https://generable.github.io/bmstate/reference/JointData.md)
  : Joint data class (R6 class)
- [`MultistateModel`](https://generable.github.io/bmstate/reference/MultiStateModel.md)
  : Main model class
- [`MultistateModelFit`](https://generable.github.io/bmstate/reference/MultistateModelFit.md)
  : Minimal fit class
- [`MultistateSystem`](https://generable.github.io/bmstate/reference/MultistateSystem.md)
  : A multistate system with proportional spline hazards
- [`PKModel`](https://generable.github.io/bmstate/reference/PKModel.md)
  : One-compartment PK Model
- [`PSSDosingData`](https://generable.github.io/bmstate/reference/PSSDosingData.md)
  : Partially steady-state dosing data class (R6 class)
- [`PathData`](https://generable.github.io/bmstate/reference/PathData.md)
  : Path data class (R6 class)
- [`TransitionMatrix`](https://generable.github.io/bmstate/reference/TransitionMatrix.md)
  : Defines states and possible transitions
- [`as_any_event()`](https://generable.github.io/bmstate/reference/as_any_event.md)
  : PathData to time-to-event data format for any state other than null
  state
- [`as_single_event()`](https://generable.github.io/bmstate/reference/as_single_event.md)
  : PathData to time-to-event data format with a single event
- [`as_survival()`](https://generable.github.io/bmstate/reference/as_survival.md)
  : PathData to event-free survival format
- [`bmstate-package`](https://generable.github.io/bmstate/reference/bmstate-package.md)
  [`bmstate`](https://generable.github.io/bmstate/reference/bmstate-package.md)
  : The 'bmstate' package.
- [`c_index()`](https://generable.github.io/bmstate/reference/c_index.md)
  : Compute concordance index
- [`create_msm()`](https://generable.github.io/bmstate/reference/create_msm.md)
  : Create a multistate model
- [`create_scoring_df()`](https://generable.github.io/bmstate/reference/create_scoring_df.md)
  : Create a data frame to evaluate state visit probability predictions
- [`create_stan_data()`](https://generable.github.io/bmstate/reference/create_stan_data.md)
  : Creating Stan data list
- [`create_stan_model()`](https://generable.github.io/bmstate/reference/create_stan_model.md)
  : Create the main 'Stan' model
- [`default_stan_filepath()`](https://generable.github.io/bmstate/reference/default_stan_filepath.md)
  : Get path to default Stan file
- [`df_to_pathdata()`](https://generable.github.io/bmstate/reference/df_to_pathdata.md)
  : Create 'PathData' from a data frame of one observed path per subject
- [`ensure_exposed_stan_functions()`](https://generable.github.io/bmstate/reference/ensure_exposed_stan_functions.md)
  : Expose 'Stan' functions if they are not yet exposed
- [`example_sim_setup_illnessdeath()`](https://generable.github.io/bmstate/reference/example_sim_setup_illnessdeath.md)
  : Example simulation setup
- [`fit_stan()`](https://generable.github.io/bmstate/reference/fit_stan.md)
  : Fit a model using 'Stan'
- [`generate_paths()`](https://generable.github.io/bmstate/reference/generate_paths.md)
  : Path generation for 'MultistateModelFit'
- [`msmfit_check_xpsr_norm()`](https://generable.github.io/bmstate/reference/msmfit_check_xpsr_norm.md)
  : Check exposure normalization
- [`msmfit_exposure()`](https://generable.github.io/bmstate/reference/msmfit_exposure.md)
  : Compute exposure
- [`msmfit_exposure_df()`](https://generable.github.io/bmstate/reference/msmfit_exposure_df.md)
  : Compute exposure and return as rvar in df
- [`msmfit_inst_hazard_param_draws()`](https://generable.github.io/bmstate/reference/msmfit_inst_hazard_param_draws.md)
  : Extract and reshape draws of instant hazard related parameters
- [`msmfit_log_hazard_multipliers()`](https://generable.github.io/bmstate/reference/msmfit_log_hazard_multipliers.md)
  : Compute log_hazard multipliers
- [`msmfit_pk_params()`](https://generable.github.io/bmstate/reference/msmfit_pk_params.md)
  : Evaluate PK parameters
- [`p_state_occupancy()`](https://generable.github.io/bmstate/reference/p_state_occupancy.md)
  : Solve state occupancy probabilities for each subject and draw in
  'MultistateModelFit'
- [`p_state_visit()`](https://generable.github.io/bmstate/reference/p_state_visit.md)
  : For each non-source state, compute probability of visiting it at
  least once before given time
- [`p_state_visit_per_subject()`](https://generable.github.io/bmstate/reference/p_state_visit_per_subject.md)
  : For each subject, compute probability of visiting a given state at
  least once before given time
- [`pk_2cpt_pss()`](https://generable.github.io/bmstate/reference/pk_2cpt_pss.md)
  : Partially steady-state PK model
- [`plot_stan_data_integral()`](https://generable.github.io/bmstate/reference/plot_stan_data_integral.md)
  : Visualize numerical integration of hazard
- [`plot_stan_data_matrix()`](https://generable.github.io/bmstate/reference/plot_stan_data_matrix.md)
  : Visualize parts of the 'Stan' data
- [`plot_state_occupancy()`](https://generable.github.io/bmstate/reference/plot_state_occupancy.md)
  : Plot mean state occupancy probabilities over time for each subject
- [`potential_covariates()`](https://generable.github.io/bmstate/reference/potential_covariates.md)
  : Look for potential covariates that affect transitions
- [`simulate_dosing()`](https://generable.github.io/bmstate/reference/simulate_dosing.md)
  : Simulate dosing data
- [`simulate_example_data()`](https://generable.github.io/bmstate/reference/simulate_example_data.md)
  : Simulate example data
- [`solve_time_evolution()`](https://generable.github.io/bmstate/reference/solve_time_evolution.md)
  : Solve the Kolmogorov forward equation of a Markovian multistate
  system
- [`solve_trans_prob_matrix()`](https://generable.github.io/bmstate/reference/solve_trans_prob_matrix.md)
  : Solve the transition probability matrix for each given output time
- [`split_data()`](https://generable.github.io/bmstate/reference/split_data.md)
  : Train-test split
- [`transmat_full()`](https://generable.github.io/bmstate/reference/transmat.md)
  [`transmat_comprisk()`](https://generable.github.io/bmstate/reference/transmat.md)
  [`transmat_survival()`](https://generable.github.io/bmstate/reference/transmat.md)
  [`transmat_illnessdeath()`](https://generable.github.io/bmstate/reference/transmat.md)
  [`transmat_progression()`](https://generable.github.io/bmstate/reference/transmat.md)
  [`transmat_diamond()`](https://generable.github.io/bmstate/reference/transmat.md)
  : Create a common transition matrix
