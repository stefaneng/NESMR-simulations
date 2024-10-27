library(GWASBrewer)
devtools::load_all('/nfs/turbo/sph-jvmorr/NESMR/esmr')
B_lower <- lower.tri(G) + 0

# Calculate LD scores
# dat <- GWASBrewer:::calc_ld_scores(dat, R_LD = ld_mat_list)

# Select on observed p-values (Winner's curse)
Zobs <- with(dat, beta_hat/s_estimate)
pval_obs <- 2*pnorm(-abs(Zobs))
minp <- apply(pval_obs, 1, min)

# ld pruning
dat$ld_list_minp <- GWASBrewer::sim_ld_prune(dat, R_LD = GWASBrewer::ld_mat_list, pvalue = minp)

ix <- dat$ld_list_minp
minp <- apply(pval_obs[ix, ], 1, min)
dat$ix1 <- dat$ix1_wc <- ix[which(minp < 5e-8)]

nesmr_model <- with(
  dat,
  esmr(
    beta_hat_X = beta_hat,
    se_X = s_estimate,
    variant_ix = ix1_wc,
    G = diag(5), # required for network problem
    direct_effect_template = B_lower,
    max_iter = 300))