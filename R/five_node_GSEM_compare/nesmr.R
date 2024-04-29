library(GWASBrewer)
devtools::load_all('/nfs/turbo/sph-jvmorr/NESMR/esmr')
B_lower <- lower.tri(G) + 0

# Error in inverting total to direct
nesmr_model <- with(
  dat,
  esmr(
    beta_hat_X = beta_hat,
    se_X = s_estimate,
    variant_ix = ix1,
    G = diag(5), # required for network problem
    direct_effect_template = B_lower,
    max_iter = 300))