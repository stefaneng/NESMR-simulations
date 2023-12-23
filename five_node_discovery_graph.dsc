sim_five_node: five_node_simulate.R
  $loglik_results: loglik_results
  $results_df: results_df
  $seed: DSC_SEED

sim_five_node_backward: five_node_simulate_backward.R
  $loglik_results: loglik_results
  $results_df: results_df
  $backward_select_edges: backward_select_edges
  $seed: DSC_SEED

DSC:
  define:
    simulate: sim_five_node_backward
  run: simulate
  exec_path: R
  output: five_node_discovery_backward
