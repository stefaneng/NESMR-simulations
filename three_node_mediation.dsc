sim_three_node_mediation_no_pleiotropy: three_node_mediation.R
  $true_model: true_model
  $incorrect_model: incorrect_model
  $dm_pvalue: dm_pvalue
  $lrt_pvalue: lrt_pvalue

DSC:
  define:
    simulate: sim_three_node_mediation_no_pleiotropy
  run: simulate
  exec_path: R
  output: three_node_mediation_no_LD
