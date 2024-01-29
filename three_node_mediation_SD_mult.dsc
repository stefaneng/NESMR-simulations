SD_mult_no_pleiotropy_pval_threshold: three_node_mediation_SD_multiplier.R
  $true_model: true_model
  $incorrect_model: incorrect_model
  $true_ll: true_ll
  $incorrect_ll: incorrect_ll
  $dm_pvalue: dm_pvalue
  $lrt_pvalue: lrt_pvalue
  $SD_multiplier: SD_multiplier

DSC:
  define:
    simulate: SD_mult_no_pleiotropy_pval_threshold
  run: simulate
  exec_path: R
  output: three_node_mediation_no_LD_test
