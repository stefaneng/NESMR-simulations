r2_threshold_collider: three_node_collider.R
  $true_model: true_model
  $incorrect_model: incorrect_model
  $dm_pvalue: dm_pvalue
  $lrt_pvalue: lrt_pvalue
  $log10_lrt_pvalue: log10_lrt_pvalue
  $log10_dm_pvalue: log10_dm_pvalue

DSC:
  define:
    simulate: r2_threshold_collider
  run: simulate
  exec_path: R
  output: three_node_collider
