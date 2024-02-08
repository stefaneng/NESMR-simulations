#!/usr/bin/env dsc

GWAS_sim_five_node: five_node_simulate.R
  $GWAS_dat: dat
  $seed: DSC_SEED

discovery_algo: five_node_simulate_backward.R
  dat: $GWAS_dat
  $loglik_results: loglik_results
  $results_df: results_df
  $backward_results: backward_results
  $backward_select_edges: backward_select_edges
  $backward_select_pvals: backward_select_pvals
  $backward_mod_results: backward_mod_results
  $discovery_model: discovery_model
  $last_backward_mod: last_backward_mod
  $backward_select_adj_mat: backward_select_adj_mat
  $seed: DSC_SEED

DSC:
  define:
    simulate: GWAS_sim_five_node
    analyze: discovery_algo
  run: simulate * analyze
  exec_path: R
  output: five_node_discovery_no_LD_winners_curse
