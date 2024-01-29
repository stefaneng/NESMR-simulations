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
    sporadic_pleiotropy = FALSE,
    est_s = TRUE
)

Z <- with(dat, beta_hat/s_estimate)
dat$pval <- 2*pnorm(-abs(Z))

SD_multiplier <- c(1, 1.5, 2, 2.5, 3)

true_ll <- list()
incorrect_ll <- list()
true_model <- list()
incorrect_model <- list()
lrt_pvalue <- list()
dm_pvalue <- list()

# min_pval <- min(dat$pval)

for (i in seq_along(SD_multiplier)) {
    SD_mult <- SD_multiplier[[i]]

    minp <- apply(dat$pval, 1, min)
    ix1 <- which(minp < 5e-8)

    true_model[[i]] <- with(dat,
            esmr(
            beta_hat_X = beta_hat[ix1,],
            se_X = s_estimate[ix1,] * SD_mult,
            pval_thresh = 1,
            G = diag(3),
            direct_effect_template = B_true,
            max_iter = 300))

    incorrect_model[[i]] <- with(dat,
            esmr(
              beta_hat_X = beta_hat[ix1,],
              se_X = s_estimate[ix1,] * SD_mult,
              pval_thresh = 1,
              G = diag(3),
              direct_effect_template = B_inc_mediation,
              max_iter = 300))

    true_ll[[i]] <- logLik.esmr(true_model[[i]])
    incorrect_ll[[i]] <- logLik.esmr(incorrect_model[[i]])

    lrt_pvalue[[i]] <- pchisq(
        - 2 * (true_ll[[i]] - incorrect_ll[[i]]),
        df = 1,
        lower.tail = FALSE
      )

    dm_pvalue[[i]] <- exp(incorrect_model[[i]]$pvals_dm[test_idx])
    }
print(lrt_pvalue)
print(dm_pvalue)
