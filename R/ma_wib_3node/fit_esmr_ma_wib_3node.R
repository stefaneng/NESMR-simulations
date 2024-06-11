library(esmr)
library(dplyr)

extract_edge_info <- function(model, sim_beta, name = NULL) {
  if (is.null(name)) {
    model_name <- as.character(as.list(match.call()[-1])$model)
  } else {
    model_name <- name
  }

  if (inherits(model, "try-error")) {
    write("Cannot extract edge info from model", stderr())
    write(model, stderr())
    return(NULL)
  }
  beta <- extract_lower_triangular(model$direct_effects)
  bias <- beta - sim_beta
  log10_pvals <- - extract_lower_triangular(model$pvals_dm) / log(10)

  res <- data.frame(
    beta,
    bias,
    log10_pvals,
    model = model_name)
  res$edge <- rownames(res)
  rownames(res) <- NULL
  res
}

B_lower <- lower.tri(diag(3)) + 0

lambda <- qnorm(1 - alpha / 2)
## Winner's cursed variant selection
valid_trait_ix <- -1
# valid_trait_ix <
Z_cursed <- with(dat, beta_hat / s_estimate)
pval_cursed <- 2 * pnorm(-abs(Z_cursed))
minp_cursed <- apply(pval_cursed[, valid_trait_ix, drop = FALSE], 1, min)
cursed_ix <- which(minp_cursed < alpha)

## True variant selection
Z_true <- with(dat, beta_marg / se_beta_hat)
pval_true <- 2 * pnorm(-abs(Z_true))
minp_true <- apply(pval_true[, -1, drop = FALSE], 1, min)
true_ix <- which(minp_true < alpha)

cursed_model <- try(with(dat, esmr::esmr(
  beta_hat_X = beta_hat,
  se_X = s_estimate,
  variant_ix = cursed_ix,
  G = diag(3),
  direct_effect_template = B_lower,
  max_iter = 300)
), silent = TRUE)

true_model <- try(with(dat, esmr::esmr(
  beta_hat_X = beta_hat,
  se_X = s_estimate,
  variant_ix = true_ix,
  G = diag(3),
  direct_effect_template = B_lower,
  max_iter = 300)
), silent = TRUE)

## Ma adjustment
ma_models <- lapply(eta, function(.eta) {
  noise_N <- Reduce('*', dim(dat$beta_hat[, -1, drop = FALSE]))
  W <- rnorm(noise_N, 0, sd = .eta)
  snp_select <- abs(Z_cursed[, -1, drop = FALSE] + W)

  rand_ix <- which(snp_select > lambda, arr.ind = TRUE)

  # Add one for the first column
  rand_ix[, 2] <- rand_ix[, 2] + 1

  # Create a new data object with the same structure as dat

  ma_adj_dat <- list()
  ma_adj_dat$beta_hat <- dat$beta_hat
  ma_adj_dat$s_estimate <- dat$s_estimate

  # Only adjust the significant SNP values
  unbias_SNPs <- esmr::snp_beta_rb(
    beta = dat$beta_hat[rand_ix],
    se_beta = dat$s_estimate[rand_ix],
    alpha = alpha,
    eta = .eta
  )

  ma_adj_dat$beta_hat[rand_ix] <- unbias_SNPs$beta_rb
  ma_adj_dat$s_estimate[rand_ix] <- unbias_SNPs$se_rb

  # Check that the variance is not NA
  valid_ix <- !is.na(unbias_SNPs$se_rb)
  rand_ix <- rand_ix[valid_ix, ]

  rand_ix_flat <- unique(rand_ix[, 1])

  ## Fit NESMR on the Ma adjusted snps
  try(with(ma_adj_dat, esmr::esmr(
    beta_hat_X = beta_hat,
    se_X = s_estimate,
    variant_ix = rand_ix_flat,
    G = diag(3),
    direct_effect_template = B_lower,
    max_iter = 300)
  ), silent = TRUE)
})

cursed_mod_results <- extract_edge_info(cursed_model, true_beta)
true_mod_results <- extract_edge_info(true_model, true_beta)

ma_mod_results <- lapply(1:seq_along(eta), function(i) {
  cbind.data.frame(
    extract_edge_info(ma_models[[i]], true_beta, name = paste0("ma_", i)),
    eta = eta[i]
  )
})

sim_results <- bind_rows(
  cursed_mod_results,
  true_mod_results,
  ma_mod_results
)

# n_variants <- lapply(
#   list(true_ix, cursed_ix, rand_ix_flat), length
# )
