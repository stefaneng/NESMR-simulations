#nohup dsc --replicate 100  --host config.yml -c 4 xxxxx.dsc > xxxxx.out &



DSC:
  define:
    analyze: nesmr_true, genomic_sem, nesmr_discovery, mvmr_true
  run: simulate * analyze
  exec_path: R
  output: 2024-04-29_five_node

simulate: five_node_GSEM_compare/five_node_simulate_LD_2.R
  n: 40000
  J: 500000
  pi: 0.001
  x: 0.05, 0.1, 0.5, 1
  h2: c(0.5, 0.3, 0.25, 0.4, 0.3)
  $N: n
  $G: G
  $dat: dat
  $seed: DSC_SEED

genomic_sem: ldsc.R + five_node_GSEM_compare/compute_ldsc.R + five_node_GSEM_compare/genomic_sem_fit_2.R
  dat: $dat
  N: $N
  $genomic_sem_fit: genomic_sem_fit
  $gsem_result: gsem_result
  $seed: DSC_SEED

nesmr_true: five_node_GSEM_compare/nesmr_true.R
  dat: $dat
  G: $G
  $nesmr_model: nesmr_model
  $nesmr_result: nesmr_result
  $seed: DSC_SEED

mvmr_true: five_node_GSEM_compare/mvmr_true.R
  dat: $dat
  $mvmr_result: mvmr_result
  $seed: DSC_SEED

nesmr_discovery: five_node_GSEM_compare/nesmr_full_alg.R
  dat: $dat
  threshold1: 0.05, 0.1, 1
  $discovery_res: mvmr_res_all
  $seed: DSC_SEED
