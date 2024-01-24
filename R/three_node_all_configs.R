renv::load('/nfs/turbo/sph-jvmorr/NESMR/simulations')
devtools::load_all('/nfs/turbo/sph-jvmorr/NESMR/esmr')

library(GWASBrewer)

G <- matrix(
  c(0, 0, 0,
    sqrt(0.4), 0, 0,
    sqrt(0.2), 0, 0),
  nrow = 3,
  byrow = TRUE
)

B_true <- G
B_true[!B_true == 0] <- 1

h2 <- c(0.5, 0.3, 0.25)
## simulate summary statistics
data(ld_mat_list)
data(AF)
    dat <- GWASBrewer::sim_mv(
    G = G,
    N = 40000,
    J = 5e5,
    h2 = h2,
    pi = 500/5e5,
    R_LD = GWASBrewer::ld_mat_list,
    af = AF,
    est_s = TRUE
)

Z <- with(dat, beta_hat/s_estimate)
dat$pval <- 2*pnorm(-abs(Z))
minp <- apply(dat$pval, 1, min)

# ld pruning

dat$ld_list_minp <- GWASBrewer::sim_ld_prune(dat, R_LD = GWASBrewer::ld_mat_list, pvalue = minp,  r2_thresh = )

minp <- apply(dat$pval, 1, min)

ix <- dat$ld_list_minp
minp <- apply(dat$pval[ix, ], 1, min)
ix1 <- which(minp < 5e-8)

# Incorrect configurations
B_alt1 <- matrix(
  c(0, 0, 0,
    1, 0, 0,
    0, 1, 0),
  nrow = 3,
  byrow = TRUE
)

B_alt2 <- matrix(
  c(0, 0, 0,
    0, 0, 0,
    1, 1, 0),
  nrow = 3,
  byrow = TRUE
)

all_B <- list(
    B_true, B_alt1, B_alt2, t(B_true), t(B_alt1), t(B_alt2)
)

all_models <- lapply(all_B, function(B) {
  with(dat,
        esmr(
        beta_hat_X = beta_hat[ix,],
        se_X = s_estimate[ix,],
        ix1 = ix1,
        G = diag(3),
        direct_effect_template = B,
        max_iter = 300))
  })

rm(all_models)

model_loglik <- lapply(all_models, logLik.esmr)
pval_dm <- lapply(all_models, function(x) {
    x$pvals_dm[3,2]
  })
