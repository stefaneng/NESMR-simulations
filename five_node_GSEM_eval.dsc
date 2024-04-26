#!/usr/bin/env dsc

GWAS_sim_five_node: five_node_GSEM_compare/five_node_simulate.R
  n: 10000, 40000
  J: 500000
  pi: 0.001
  $GWAS_dat: dat
  $seed: DSC_SEED

comptute_ldsc: ldsc.R + five_node_GSEM_compare/compute_ldsc.R
  dat: $GWAS_dat
  $ldsc_res: ldsc_res
  $seed: DSC_SEED

genomic_sem: five_node_GSEM_compare/genomic_sem_fit.R
  ldsc_res: ldsc_res
  $sem_mod_DWLS: sem_mod_DWLS
  $sem_mod_ML: sem_mod_ML
  $genomic_sem_fit: genomic_sem_fit
  $seed: DSC_SEED

DSC:
  define:
    simulate: GWAS_sim_five_node
    ldsc: comptute_ldsc
    analyze: nesmr, genomic_sem
  run: simulate * analyze
  exec_path: R
  output: five_node_GSEM_eval
