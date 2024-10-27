DSC:
  run: simulate * fit_esmr
  exec_path: R
  output: simulation_results/2024-06-11_ma_wib_3node

simulate:  renv.R + utils.R + ma_wib_3node/sim_ma_wib_3node.R
  N: 30000, 2e5
  J: 5000
  pi_J: 0.1
  effect_scale_factor: 0.2
  h2: 0.3
  graph_type: "correlated", "mediation"
  $true_beta: true_beta
  $dat: dat
  $seed: DSC_SEED

fit_esmr: renv.R + utils.R + ma_wib_3node/fit_esmr_ma_wib_3node.R
  eta: c(0.5, 1)
  alpha: 1, 1e-3, 1e-5, 1e-8, 1e-10
  dat: $dat
  true_beta: $true_beta
  $sim_results: sim_results
  $cursed_mod_results: cursed_mod_results
  $true_mod_results: true_mod_results
  $ma_mod_results: ma_mod_results
  $cursed_model: cursed_model
  $true_model: true_model
  $ma_models: ma_models
