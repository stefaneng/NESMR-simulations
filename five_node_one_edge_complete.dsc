#!/usr/bin/env dsc

five_node_one_edge_complete2: five_node_complete_sim.R
  $model_ll: model_ll
  $results_df: results_df

DSC:
  define:
    simulate: five_node_one_edge_complete2
  run: simulate
  exec_path: R
  output: five_node_one_edge_complete
