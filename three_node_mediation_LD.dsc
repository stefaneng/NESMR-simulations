LD_r2_pval_threshold: three_node_mediation_LD.R
  $true_model: true_model
  $incorrect_model: incorrect_model
  $dm_pvalue: dm_pvalue
  $lrt_pvalue: lrt_pvalue
  $min_pval: min_pval

DSC:
  define:
    simulate: LD_r2_pval_threshold
  run: simulate
  exec_path: R
  output: three_node_mediation_LD
