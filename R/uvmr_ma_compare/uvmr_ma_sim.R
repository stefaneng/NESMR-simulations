library(esmr)

ma_rivw_est <- function(beta_x, se_x, beta_y, se_y) {
  denom <- sum((beta_x^2 - se_x^2) / se_y^2)
  beta_rivw <- sum(beta_x * beta_y / se_y^2) / denom
  se_num <- sum(
    (beta_x * beta_y - beta_rivw * (beta_x^2 - se_x^2))^2 / se_y^4
  )
  se_rivw <- se_num / denom^2
  return(list(
    beta_rivw = beta_rivw, se_rivw = se_rivw
  ))
}

ivw_est <- function(beta_x, se_x, beta_y, se_y) {
  beta_ivw <- sum(beta_x * beta_y / se_y^2) / sum(beta_x^2 / se_y^2)
  se_ivw <- 1 / sqrt(sum(beta_x^2 / se_y^2))

  # wls_res <- lm(beta_y ~ beta_x - 1, weights = 1 / se_y^2)
  # mr_res <- mr_ivw(
  #   mr_input(bx = beta_x, bxse = se_x, by = beta_y, byse = se_y)
  #   )

  return(
    list(
      beta_ivw = beta_ivw, se_ivw = se_ivw
    ))
  # mr.divw:::ivw(
  #   beta.exposure = beta_x, se.exposure = se_x, beta.outcome = beta_y, se.outcome = se_y
  #   )[c('beta.hat', 'beta.se')]
  # ))
}

uvmr_rand_wc_sim <- function(
    eta = 1,
    alpha = 5e-8,
    lambda = qnorm(1 - alpha / 2),
    true_beta = 0.1,
    n = 10000,
    J = 5e5,
    h2 = 0.3,
    pi_J = 500 / J
) {

  G <- matrix(
    c(0, 0,
      true_beta, 0
    ),
    nrow = 2,
    byrow = TRUE
  )

  # G <- scale_effects * G
  B_true <- G
  B_true[!B_true == 0] <- 1
  B_lower <- lower.tri(B_true) + 0

  dat <- GWASBrewer::sim_mv(
    G = G,
    N = n,
    J = J,
    h2 = h2,
    pi = pi_J,
    sporadic_pleiotropy = FALSE,
    est_s = TRUE
  )

  B_lower <- lower.tri(G) + 0

  # Winner's cursed estimates
  # TODO: is this actually an issue in how we select the variants?
  # In the UVMR case we only case about G -> X but not G -> Y
  # But in NESMR we select variants on and G -> X_j...
  Z_cursed_XY <- with(dat, beta_hat / s_estimate)
  Z_cursed_X <- Z_cursed_XY[, 2]
  pval_cursed <- 2 * pnorm(-abs(Z_cursed_XY))
  minp_cursed <- apply(pval_cursed, 1, min)
  cursed_ix_XY <- which(minp_cursed < alpha)
  cursed_ix_X <- which(pval_cursed[, 2] < alpha)

  # Select only on G -> X variants
  cursed_esmr_X <- try(with(dat, esmr::esmr(
    beta_hat_X = beta_hat,
    se_X = s_estimate,
    variant_ix = cursed_ix_X,
    G = diag(2),
    direct_effect_template = B_lower,
    max_iter = 300)
  ), silent = TRUE)

  # Select on either G -> X or G -> Y variants (incorrect for UVMR case..)
  cursed_esmr_XY <- try(with(dat, esmr::esmr(
    beta_hat_X = beta_hat,
    se_X = s_estimate,
    variant_ix = cursed_ix_XY,
    G = diag(2),
    direct_effect_template = B_lower,
    max_iter = 300)
  ), silent = TRUE)

  ## IVW Estimator:
  ivw_cursed <- ivw_est(
    beta_x = dat$beta_hat[cursed_ix_X, 2],
    se_x = dat$s_estimate[cursed_ix_X, 2],
    beta_y = dat$beta_hat[cursed_ix_X, 1],
    se_y = dat$s_estimate[cursed_ix_X, 1]
  )

  # Step 1: Randomized instrument selection
  Z_noise <- rnorm(J, 0, sd = eta)
  snp_select <- abs(Z_cursed_X + Z_noise)

  rand_ix <- which(snp_select > lambda)

  # Create a new data object with the same structure as dat
  ini_scale_factor <- 1 + 1/eta^2
  ini_est_dat <- dat
  # This uses the true SE(beta_hat)
  # TODO: Should this only be for the selected SNPs?
  ini_est_dat$beta_hat[, 2] <- dat$beta_hat[, 2] - Z_noise * dat$se_beta_hat[, 2] / eta^2
  ini_est_dat$s_estimate[, 2] <- dat$s_estimate[, 2] * sqrt(ini_scale_factor)
  ini_est_dat$se_beta_hat[, 2] <- dat$se_beta_hat[, 2] * sqrt(ini_scale_factor)

  ini_N_rescale <- n / ini_scale_factor
  # Use resample_sumstats to make fair comparison
  ini_rescale_dat <- GWASBrewer::resample_sumstats(
    dat,
    N = ini_N_rescale,
    est_s = TRUE
  )

  ini_cursed_ix_X <- which(
    abs(ini_rescale_dat$beta_hat[, 2] / ini_rescale_dat$s_estimate[, 2]) > lambda
  )

  ini_ivw_cursed <- ivw_est(
    beta_x = ini_rescale_dat$beta_hat[ini_cursed_ix_X, 2],
    se_x = ini_rescale_dat$s_estimate[ini_cursed_ix_X, 2],
    beta_y = ini_rescale_dat$beta_hat[ini_cursed_ix_X, 1],
    se_y = ini_rescale_dat$s_estimate[ini_cursed_ix_X, 1]
  )

  ini_reduced_N_esmr <- try(with(ini_rescale_dat, esmr::esmr(
    beta_hat_X = beta_hat,
    se_X = s_estimate,
    variant_ix = ini_cursed_ix_X,
    G = diag(2),
    direct_effect_template = B_lower,
    max_iter = 300)
  ), silent = TRUE)

  ini_esmr <- try(with(ini_est_dat, esmr::esmr(
    beta_hat_X = beta_hat,
    se_X = s_estimate,
    variant_ix = rand_ix,
    G = diag(2),
    direct_effect_template = B_lower,
    max_iter = 300)
  ), silent = TRUE)

  ### Ma's winner curse adjustment

  # Only adjust the significant SNP values
  unbias_SNPs <- esmr::snp_beta_rb(
    beta = dat$beta_hat[rand_ix, 2],
    se_beta = dat$s_estimate[rand_ix, 2],
    lambda = lambda,
    eta = eta
  )

  ma_adj_dat <- list()
  ma_adj_dat$beta_hat <- dat$beta_hat
  ma_adj_dat$s_estimate <- dat$s_estimate
  # Only update the significant SNPs beta and SE
  ma_adj_dat$beta_hat[rand_ix, 2] <- unbias_SNPs$beta_rb
  ma_adj_dat$s_estimate[rand_ix, 2] <- unbias_SNPs$se_rb

  ma_scale_factor <- mean(ma_adj_dat$s_estimate[rand_ix, 2] / dat$s_estimate[rand_ix, 2], na.rm = TRUE)

  # TODO: Can we resample with different Ns?
  ma_N_rescale <- n / ma_scale_factor
  ma_rescale_dat <- GWASBrewer::resample_sumstats(
    dat,
    N = ma_N_rescale,
    est_s = TRUE
  )

  ma_cursed_ix_X <- which(
    abs(ma_rescale_dat$beta_hat[, 2] / ma_rescale_dat$s_estimate[, 2]) > lambda
  )

  ma_ivw_cursed <- ivw_est(
    beta_x = ma_rescale_dat$beta_hat[ma_cursed_ix_X, 2],
    se_x = ma_rescale_dat$s_estimate[ma_cursed_ix_X, 2],
    beta_y = ma_rescale_dat$beta_hat[ma_cursed_ix_X, 1],
    se_y = ma_rescale_dat$s_estimate[ma_cursed_ix_X, 1]
  )

  # names(ma_ivw_cursed) <- paste0('ma_', names(ivw_cursed))

  # Fit model with Ma variants

  # Use the randomized selected variants in esmr
  ma_esmr <- try(with(ma_adj_dat, esmr::esmr(
    beta_hat_X = beta_hat,
    se_X = s_estimate,
    variant_ix = rand_ix,
    G = diag(2),
    direct_effect_template = B_lower,
    max_iter = 300)
  ), silent = TRUE)

  ma_reduced_N_esmr <- try(with(ma_rescale_dat, esmr::esmr(
    beta_hat_X = beta_hat,
    se_X = s_estimate,
    variant_ix = ma_cursed_ix_X,
    G = diag(2),
    direct_effect_template = B_lower,
    max_iter = 300)
  ), silent = TRUE)

  ## Fit the Ma RIVR estimator
  ma_rivw <- ma_rivw_est(
    beta_x = ma_adj_dat$beta_hat[rand_ix, 2],
    se_x = ma_adj_dat$s_estimate[rand_ix, 2],
    beta_y = ma_adj_dat$beta_hat[rand_ix, 1],
    se_y = ma_adj_dat$s_estimate[rand_ix, 1]
  )

  extract_results <- function(m) {
    if (class(m) != 'try-error') {
      return(list(
        beta = m$direct_effects[2,1],
        se_beta_log10 = - m$pvals_dm[2,1] / log(10)
      ))
    } else {
      return(list(
        beta = NA,
        se_beta_log10 = NA
      ))
    }
  }

  ## Return results
  res <- lapply(
    list(
      cursed_esmr_X = cursed_esmr_X,
      cursed_esmr_XY = cursed_esmr_XY,
      ini_reduced_N_esmr = ini_reduced_N_esmr,
      ini_esmr = ini_esmr,
      ma_esmr = ma_esmr,
      ma_reduced_N_esmr = ma_reduced_N_esmr
    ),
    extract_results
  )

  res$ivw_cursed <- ivw_cursed
  res$ini_ivw_cursed <- ini_ivw_cursed
  res$ma_ivw_cursed <- ma_ivw_cursed
  res$ma_rivw <- ma_rivw
  return(res)
}
