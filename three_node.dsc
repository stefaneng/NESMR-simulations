sim_three_node: three_node_simulate.R
  $true_ll: true_ll
  $full_ll: full_ll
  $lrt_pvalue: lrt_pvalue
  $pval_dm: pval_dm
  $seed: DSC_SEED

DSC:
  define:
    simulate: sim_three_node
  run: simulate
  exec_path: R
  output: three_node_output