B_lower <- lower.tri(G) + 0

lambda <- qnorm(1 - alpha / 2)

# Winner's cursed estimates
Z_cursed <- with(dat, beta_hat / s_estimate)
pval_cursed <- 2 * pnorm(-abs(Z_cursed))
minp_cursed <- apply(pval_cursed, 1, min)
cursed_ix <- which(minp_cursed < alpha)

cursed_model <- try(with(dat, esmr::esmr(
  beta_hat_X = beta_hat,
  se_X = s_estimate,
  pval_thresh = alpha,
  G = diag(3),
  direct_effect_template = B_lower,
  max_iter = 300)
), silent = TRUE)

# Step 1: Randomized instrument selection
noise <- Reduce('*', dim(dat$beta_hat))
snp_select <- abs(Z_cursed + rnorm(noise, 0, sd = eta))

rand_ix <- which(snp_select > lambda, arr.ind = TRUE)

# Create a new data object with the same structure as dat

ma_adj_dat <- list()
ma_adj_dat$beta_hat <- dat$beta_hat
ma_adj_dat$s_estimate <- dat$s_estimate

# Only adjust the significant SNP values
unbias_SNPs <- esmr::snp_beta_rb(
  beta = dat$beta_hat[rand_ix],
  se_beta = dat$s_estimate[rand_ix]
)

# Only update the significant SNPs beta and SE
ma_adj_dat$beta_hat[rand_ix] <- unbias_SNPs$beta_rb
ma_adj_dat$s_estimate[rand_ix] <- unbias_SNPs$se_rb

# Fit model with Ma variants

# Use the randomized selected variants in esmr
ma_ix <- unique(rand_ix[, 1])

ma_SNP_model <- try(with(ma_adj_dat, esmr::esmr(
  beta_hat_X = beta_hat,
  se_X = s_estimate,
  variant_ix = ma_ix,
  G = diag(3),
  direct_effect_template = B_lower,
  max_iter = 300)
), silent = TRUE)
