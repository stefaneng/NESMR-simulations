renv::load('/nfs/turbo/sph-jvmorr/NESMR/simulations')
devtools::load_all('/nfs/turbo/sph-jvmorr/NESMR/esmr')

library(GWASBrewer)

G <- matrix(
  c(0, 0, 0,
    sqrt(0.4), 0, 0,
    0, sqrt(0.2), 0),
  nrow = 3,
  byrow = TRUE
)

B_true <- G
B_true[!B_true == 0] <- 1

# Incorrect configurations
B_inc_mediation <- matrix(
  c(0, 0, 0,
    1, 0, 0,
    1, 1, 0),
  nrow = 3,
  byrow = TRUE
)

test_idx <- which(B_inc_mediation - B_true > 0, arr.ind = TRUE)

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

Z <- with(dat, beta_hat/s_estimate)
dat$pval <- 2*pnorm(-abs(Z))

pval_thresholds <- c(5e-8, 5e-9, 5e-10)
r2_thresholds <- c(0.1, 0.01, 0.001, 0.0001)

init_list <- replicate(
  length(r2_thresholds), vector("list", length(pval_thresholds)),
  simplify = FALSE)

true_model <- init_list
incorrect_model <- init_list
lrt_pvalue <- init_list
dm_pvalue <- init_list
log10_lrt_pvalue <-init_list
log10_dm_pvalue <- init_list

min_pval <- min(dat$pval)

for (i in seq_along(r2_thresholds)) {
  for (j in seq_along(pval_thresholds)) {
    pval_thres <- pval_thresholds[[j]]
    r2_thres <- r2_thresholds[[i]]

    if (min_pval > pval_thres) break;

    Z <- with(dat, beta_hat/s_estimate)
    dat$pval <- 2*pnorm(-abs(Z))
    minp <- apply(dat$pval, 1, min)
    # ld pruning

    dat$ld_list_minp <- GWASBrewer::sim_ld_prune(
      dat, R_LD = GWASBrewer::ld_mat_list, pvalue = minp, r2_thresh = r2_thres, pval_thresh = pval_thres)

    minp <- apply(dat$pval, 1, min)

    ix1 <- dat$ld_list_minp

    true_model[[i]][[j]] <- with(dat,
                                 esmr(
                                   beta_hat_X = beta_hat[ix1,],
                                   se_X = s_estimate[ix1,],
                                   pval_thresh = 1,
                                   G = diag(3),
                                   direct_effect_template = B_true,
                                   max_iter = 300))

    incorrect_model[[i]][[j]] <- with(dat,
                                      esmr(
                                        beta_hat_X = beta_hat[ix1,],
                                        se_X = s_estimate[ix1,],
                                        pval_thresh = 1,
                                        G = diag(3),
                                        direct_effect_template = B_inc_mediation,
                                        max_iter = 300))

    true_ll <- logLik.esmr(true_model[[i]][[j]])
    incorrect_ll <- logLik.esmr(incorrect_model[[i]][[j]])

    lrt_pvalue[[i]][[j]] <- pchisq(
      - 2 * (true_ll - incorrect_ll),
      df = 1,
      lower.tail = FALSE
    )

    log10_lrt_pvalue[[i]][[j]] <- pchisq(
      - 2 * (true_ll - incorrect_ll),
      df = 1,
      lower.tail = FALSE,
      log.p = TRUE
    ) / log(10)

    log10_dm_pvalue[[i]][[j]] <- incorrect_model[[i]][[j]]$pvals_dm[test_idx] / log(10)
    dm_pvalue[[i]][[j]] <- exp(incorrect_model[[i]][[j]]$pvals_dm[test_idx])
  }
}
print(lrt_pvalue)
print(dm_pvalue)