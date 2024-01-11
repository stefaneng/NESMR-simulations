#!/usr/bin/env dsc

GWAS_sim_five_node: five_node_simulate.R
  $GWAS_dat: dat
  $seed: DSC_SEED

one_edge_test: five_node_one_extra.R
  dat: $GWAS_dat
  $true_model: true_model
  $true_ll: true_ll
  $extra_edge_df: extra_edge_df
  $seed: DSC_SEED

DSC:
  define:
    simulate: GWAS_sim_five_node
    analyze: one_edge_test
  run: simulate * analyze
  exec_path: R
  output: five_node_five_one
