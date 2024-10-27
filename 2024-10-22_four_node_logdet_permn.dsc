DSC:
  run: simulate * process
  exec_path: R
  output: simulation_results/2024-10-22_4node_logdet_permn

process: R(
  library(esmr); library(NESMR.Sims); library(tidyr); library(plyr); library(dplyr);
    summary_sim_results <- logdet_process_sim_results(
      sim_results,
      G = G))
  G: $G
  sim_results: $sim_results
  $summary_sim_results: summary_sim_results
  $seed: DSC_SEED

simulate:  R(
  library(esmr); library(NESMR.Sims);
  new_G <- G1 * effect_modifier;
  sim_results <- logdet_perm_fas_sim(
    G = new_G,
    alpha = alpha,
    include_all_perms = ncol(G1) <= 4,
    gwas_params = list(
        h2 = h2,
        J = J,
        N = N,
        pi = pi_J,
        sporadic_pleiotropy = TRUE
        )
    ))
  N: 20000
  J: 5000
  pi_J: 0.1
  h2: 0.3
  alpha: 5e-8
  effect_modifier: 1, 0.75
  G1: matrix(c(0, 0, 0, 0,
    0.1, 0, 0, 0,
    0, -0.05, 0, 0,
    -0.05, 0, 0.1, 0),
    nrow = 4,
    byrow = TRUE
  )
  $G: new_G
  $seed: DSC_SEED
  $sim_results: sim_results
