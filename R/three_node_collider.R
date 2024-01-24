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

# Incorrect configurations
B_inc <- matrix(
  c(0, 0, 0,
    1, 0, 0,
    1, 1, 0),
  nrow = 3,
  byrow = TRUE
)

test_idx <- which(B_inc - B_true > 0, arr.ind = TRUE)

h2 <- c(0.5, 0.3, 0.25)
## simulate summary statistics
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

true_model <- list()
incorrect_model <- list()
lrt_pvalue <- list()
dm_pvalue <- list()
log10_lrt_pvalue <- list()
log10_dm_pvalue <- list()

r2_thresholds <- c(0.1, 0.01, 0.001, 0.0001)
for (i in seq_along(r2_thresholds)) {
  r2_thresh <- r2_thresholds[[i]]

  Z <- with(dat, beta_hat/s_estimate)
  dat$pval <- 2*pnorm(-abs(Z))
  minp <- apply(dat$pval, 1, min)
  # ld pruning

  dat$ld_list_minp <- GWASBrewer::sim_ld_prune(
    dat, R_LD = GWASBrewer::ld_mat_list, pvalue = minp, r2_thresh = r2_thresh)

  minp <- apply(dat$pval, 1, min)

  ix <- dat$ld_list_minp
  minp <- apply(dat$pval[ix, ], 1, min)
  ix1 <- which(minp < 5e-8)

  # TODO: Want to run with beta_marg as well
  true_model[[i]] <- with(dat,
                          esmr(
                            beta_hat_X = beta_hat[ix, ][ix1,], # TODO: Double check [ix, ][ix1, ]
                            se_X = s_estimate[ix, ][ix1,],
                            pval_thresh = 1,
                            G = diag(3),
                            direct_effect_template = B_true,
                            max_iter = 300))

  incorrect_model[[i]] <- with(dat,
                               esmr(
                                 beta_hat_X = beta_hat[ix, ][ix1,],
                                 se_X = s_estimate[ix, ][ix1,],
                                 pval_thresh = 1,
                                 G = diag(3),
                                 direct_effect_template = B_inc,
                                 max_iter = 300))

  true_ll <- logLik.esmr(true_model[[i]])
  incorrect_ll <- logLik.esmr(incorrect_model[[i]])

  lrt_pvalue[[i]] <- pchisq(
    - 2 * (true_ll - incorrect_ll),
    df = 1,
    lower.tail = FALSE
  )

  log10_lrt_pvalue[[i]] <- pchisq(
    - 2 * (true_ll - incorrect_ll),
    df = 1,
    lower.tail = FALSE,
    log.p = TRUE
    ) / log(10)

  log10_dm_pvalue[[i]] <- incorrect_model[[i]]$pvals_dm[test_idx] / log(10)
  dm_pvalue[[i]] <- exp(incorrect_model[[i]]$pvals_dm[test_idx])
  print(dm_pvalue)
  print(lrt_pvalue)
}

print(lrt_pvalue)
print(dm_pvalue)
