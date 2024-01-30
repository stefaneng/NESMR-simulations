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

h2 <- c(0.5, 0.3, 0.25)
## simulate summary statistics
dat <- GWASBrewer::sim_mv(
  G = G,
  N = 40000,
  J = 5e5,
  h2 = h2,
  pi = 500/5e5,
  sporadic_pleiotropy = FALSE,
  est_s = TRUE
)

Z <- with(dat, beta_hat/s_estimate)
dat$pval <- 2*pnorm(-abs(Z))

pval_thresholds <- c(5e-8, 5e-9, 5e-10)

min_pval <- min(dat$pval)

for (i in seq_along(pval_thresholds)) {
  pval <- pval_thresholds[[i]]

  if (min_pval > pval) break;

  minp <- apply(dat$pval, 1, min)
  ix1 <- which(minp < 5e-8)

  beta_var <- 'beta_hat'

  # TODO: Want to run with beta_marg as well
  true_model[[i]] <- with(dat,
                          esmr(
                            beta_hat_X = get(beta_var)[ix1,],
                            se_X = s_estimate[ix1,],
                            pval_thresh = 1,
                            G = diag(3),
                            direct_effect_template = B_true,
                            max_iter = 300))

  incorrect_model[[i]] <- with(dat,
                               esmr(
                                 beta_hat_X = beta_hat[ix1,],
                                 se_X = s_estimate[ix1,],
                                 pval_thresh = 1,
                                 G = diag(3),
                                 direct_effect_template = B_inc_mediation,
                                 max_iter = 300))

  true_ll <- logLik.esmr(true_model)
  incorrect_ll <- logLik.esmr(incorrect_model)

  lrt_pvalue[[i]] <- pchisq(
    - 2 * (true_ll - incorrect_ll),
    df = 1,
    lower.tail = FALSE
  )

  dm_pvalue[[i]] <- exp(incorrect_model$pvals_dm[3, 1])
}
