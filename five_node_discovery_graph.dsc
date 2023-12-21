sim_five_node: five_node_simulate.R
  $loglik_results: loglik_results
  $results_df: results_df
  $seed: DSC_SEED

DSC:
  define:
    simulate: sim_five_node
  run: simulate
  exec_path: R
  output: five_node_discovery
