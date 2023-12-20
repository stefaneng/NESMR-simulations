sim_three_node_two_edge: three_node_all_configs.R
  $loglik: model_loglik
  $dm_log_pval: pval_dm
  $seed: DSC_SEED

DSC:
  define:
    simulate: sim_three_node_two_edge
  run: simulate
  exec_path: R
  output: three_node_output_two_edge
