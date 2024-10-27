library(GWASBrewer)
library(esmr)
#devtools::load_all('/nfs/turbo/sph-jvmorr/NESMR/esmr')
B_lower <- lower.tri(G) + 0


minp <- apply(dat$pval_true[dat$ld_list_external, ], 1, min)
ix1 <- dat$ld_list_external[minp < 5e-8]

nesmr_model <- with(
  dat,
  esmr(
    beta_hat_X = beta_hat,
    se_X = s_estimate,
    variant_ix = ix1,
    G = diag(5), # required for network problem
    direct_effect_template = B_lower,
    max_iter = 300))
