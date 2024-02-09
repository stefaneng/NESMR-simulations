renv::load('/nfs/turbo/sph-jvmorr/NESMR/simulations')
devtools::load_all('/nfs/turbo/sph-jvmorr/NESMR/esmr')
library(GWASBrewer)
# library(igraph)
library(dplyr)

threshold <- 0.05

G <- matrix(
  c(0, 0, 0, 0, 0,
    sqrt(0.3), 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, -1 * sqrt(0.1), 0, 0, 0,
    0, -1 * sqrt(0.1), sqrt(0.2), sqrt(0.25), 0),
  nrow = 5,
  byrow = 5
)

B_true <- 1 * (abs(G) > 0)

h2 <- c(0.5, 0.3, 0.25, 0.4, 0.3)
## simulate summary statistics
dat <- sim_mv(
  G = G,
  N = 40000,
  J = 5e5,
  h2 = h2,
  pi = 500/5e5,
  sporadic_pleiotropy = FALSE,
  est_s = TRUE
)

# Remove winner's curse bias
Ztrue <- with(dat, beta_marg/se_beta_hat)
pval_true <- 2*pnorm(-abs(Ztrue))
minp <- apply(pval_true, 1, min)
ix <- which(minp < 5e-8)

true_edges <- which(G != 0 & lower.tri(G), arr.ind = TRUE)
extra_edges <- which(G == 0 & lower.tri(G), arr.ind = TRUE)

 # TODO: Attach the true beta value to the results_df
results_df <- rbind.data.frame(
    cbind(as.data.frame(true_edges), correct_edge = TRUE),
    cbind(as.data.frame(extra_edges), correct_edge = FALSE)
)

results_df <- merge(results_df,
    data.frame(true_edges, beta = G[true_edges]),
    by = c('row', 'col'),
    all.x = TRUE)

results_df$beta[is.na(results_df$beta)] <- 0

# Fit true model
true_model <- with(
    dat,
    esmr(
        beta_hat_X = beta_hat[ix,],
        se_X = s_estimate[ix,],
        pval_thresh = 1,
        G = diag(5), # required for network problem
        direct_effect_template = B_true,
        max_iter = 300))

model_ll <- list()
model_ll$true_ll <- logLik.esmr(true_model)

# Add the true model p-values from true_model and effects to the results_df
true_model_edges <- which(B_true != 0, arr.ind = TRUE)
true_model_results <- cbind(
    as.data.frame(true_model_edges),
    true_model_log10_pvals = (true_model$pvals_dm / log(10))[true_model_edges],
    true_model_beta = true_model$direct_effects[true_model_edges]
)

results_df <- merge(
    results_df,
    true_model_results,
    by = c('row', 'col'),
    all.x = TRUE
)

# Fit complete model
complete_model <- with(
    dat,
    esmr(
        beta_hat_X = beta_hat[ix,],
        se_X = s_estimate[ix,],
        pval_thresh = 1,
        G = diag(5), # required for network problem
        direct_effect_template = 1 * lower.tri(G),
        max_iter = 300))

model_ll$complete_ll <- logLik.esmr(complete_model)

complete_model_edges <- which(lower.tri(G), arr.ind = TRUE)
complete_model_results <- cbind(
    as.data.frame(complete_model_edges),
    complete_model_log10_pvals = (complete_model$pvals_dm / log(10))[complete_model_edges],
    complete_model_beta = complete_model$direct_effects[complete_model_edges]
)

results_df <- merge(
    results_df,
    complete_model_results,
    by = c('row', 'col'),
    all.x = TRUE
)

# TODO: Need to add each of the one-edges to the results_df
# TODO: Build up the loglik for each of the models in separate result list

# 5 edges that we need to test
k <- nrow(extra_edges)
one_extra_models <- vector("list", k)
for(i in seq_len(k)) {
    B_one_edge <- B_true
    B_one_edge[extra_edges[i,, drop = F]] <- 1
    edge_name <- paste(extra_edges[i,], collapse = '-')

    m <- with(
        dat,
        esmr(
            beta_hat_X = beta_hat[ix,],
            se_X = s_estimate[ix,],
            pval_thresh = 1,
            G = diag(5), # required for network problem
            direct_effect_template = B_one_edge,
            max_iter = 300))

    model_ll[[edge_name]] <- logLik.esmr(m)

    model_edges <- which(B_one_edge != 0, arr.ind = TRUE)
    model_results <- cbind(
        as.data.frame(model_edges),
        model_log10_pvals = (m$pvals_dm / log(10))[model_edges],
        model_beta = m$direct_effects[model_edges]
    ) %>%
    setNames(c('row', 'col', paste0(edge_name, '_log10_pvals'), paste0(edge_name, '_beta')))

    results_df <- merge(
        results_df,
        model_results,
        by = c('row', 'col'),
        all.x = TRUE
    )
}

print(model_ll)
print(results_df)
