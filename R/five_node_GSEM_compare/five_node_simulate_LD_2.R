#renv::load('/nfs/turbo/sph-jvmorr/NESMR/simulations')
#devtools::load_all('/nfs/turbo/sph-jvmorr/NESMR/esmr')
library(GWASBrewer)

G <- matrix(
  c(0, 0, 0, 0, 0,
    sqrt(0.3), 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, -1 * sqrt(0.1), 0, 0, 0,
    0, -1 * sqrt(0.1), sqrt(0.2), sqrt(0.25), 0),
  nrow = 5,
  byrow = 5
)

G <- x*G

data(ld_mat_list)
data(AF)
## simulate summary statistics
# Parameters from dsc
# h2 <- c(0.5, 0.3, 0.25, 0.4, 0.3)
# n <- 40000
# J <- 5e5
#pi <- 500 / J
dat <- sim_mv(
  G = G,
  N = n,
  J = J,
  h2 = h2,
  pi = pi,
  R_LD = GWASBrewer::ld_mat_list,
  af = AF,
  sporadic_pleiotropy = FALSE,
  est_s = TRUE
)

# Calculate LD scores
dat <- GWASBrewer:::calc_ld_scores(dat, R_LD = ld_mat_list)

# Select on true variants
Ztrue <- with(dat, beta_marg/se_beta_hat)
Zhat <- with(dat, beta_hat/s_estimate)
dat$pval_true <- 2*pnorm(-abs(Ztrue))
dat$pval_hat <- 2*pnorm(-abs(Zhat))
minp <- apply(dat$pval_true, 1, min)
minphat <- apply(dat$pval_hat, 1, min)

# ld pruning
dat$ld_list_external <- GWASBrewer::sim_ld_prune(dat, R_LD = GWASBrewer::ld_mat_list, pvalue = minp,
                                             r2_thresh = 0.001)
dat$ld_list_internal <- GWASBrewer::sim_ld_prune(dat, R_LD = GWASBrewer::ld_mat_list, pvalue = minp,
                                             r2_thresh = 0.001)
