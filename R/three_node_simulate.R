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

# Need to trim the cycles
B_alt1 <- matrix(
  c(0, 0, 0,
    1, 0, 0,
    1, 1, 0),
  nrow = 3,
  byrow = TRUE
)

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
dat$ld_list_minp <- GWASBrewer::sim_ld_prune(dat, R_LD = GWASBrewer::ld_mat_list, pvalue = minp)

minp <- apply(dat$pval, 1, min)

ix <- dat$ld_list_minp
minp <- apply(dat$pval[ix, ], 1, min)
ix1 <- which(minp < 5e-8)

true_model <- with(dat,
                        esmr(
                        beta_hat_X = beta_hat[ix,],
                        se_X = s_estimate[ix,],
                        ix1 = ix1,
                        G = diag(3), # required for network problem
                        direct_effect_template = B_true,
                        max_iter = 300))

full_model <- with(dat,
                        esmr(
                        beta_hat_X = beta_hat[ix,],
                        se_X = s_estimate[ix,],
                        ix1 = ix1,
                        G = diag(3), # required for network problem
                        direct_effect_template = B_alt1,
                        max_iter = 300
                        )
)

true_ll <- logLik(true_model)
full_ll <- logLik(full_model)

lrt_pvalue <- pchisq(
    - 2 * (true_ll - full_ll),
    df = sum(B_alt1) - sum(B_true),
    lower.tail = FALSE
)

pval_dm <- full_ll$pvals_dm[3,2]
cat('lrt_pvalue', lrt_pvalue, 'pval_dm', pval_dm, '\n')

true_fbar <- true_model$f$fbar
full_fbar <- full_model$f$fbar
