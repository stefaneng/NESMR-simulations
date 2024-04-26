#!/usr/bin/env dsc

GWAS_sim_five_node: five_node_GSEM_compare/five_node_simulate_LD.R
  n: 10000, 40000
  J: 500000
  pi: 0.001
  $N: n
  $G: G
  h2: c(0.5, 0.3, 0.25, 0.4, 0.3)
  $GWAS_dat: dat
  $seed: DSC_SEED

compute_ldsc: ldsc.R + five_node_GSEM_compare/compute_ldsc.R
  dat: $GWAS_dat
  N: $N
  $ldsc_res: ldsc_res
  $seed: DSC_SEED

genomic_sem: five_node_GSEM_compare/genomic_sem_fit.R
  ldsc_res: $ldsc_res
  $sem_mod_DWLS: sem_mod_DWLS
  $sem_mod_ML: sem_mod_ML
  $genomic_sem_fit: genomic_sem_fit
  $seed: DSC_SEED

nesmr: five_node_GSEM_compare/nesmr.R
  dat: $GWAS_dat
  G: $G
  $nesmr_model: nesmr_model
  $seed: DSC_SEED

DSC:
  define:
    simulate: GWAS_sim_five_node
    ldsc: compute_ldsc
    analyze: nesmr, genomic_sem
  run: simulate * ldsc * analyze
  exec_path: R
  output: five_node_GSEM_eval
