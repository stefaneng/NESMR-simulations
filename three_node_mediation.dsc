no_pleiotropy_pval_threshold: three_node_mediation.R
  $true_model: true_model
  $incorrect_model: incorrect_model
  $dm_pvalue: dm_pvalue
  $lrt_pvalue: lrt_pvalue
  $min_pval: min_pval

DSC:
  define:
    simulate: no_pleiotropy_pval_threshold
  run: simulate
  exec_path: R
  output: three_node_mediation_no_LD_test
