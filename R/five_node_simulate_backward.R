renv::load('/nfs/turbo/sph-jvmorr/NESMR/simulations')
devtools::load_all('/nfs/turbo/sph-jvmorr/NESMR/esmr')
library(GWASBrewer)
library(dplyr)
library(igraph)

threshold <- 0.05

is_acyclic <- function(g) {
  tryCatch(
    !is.null(topo_sort(g)),
    error = function(e) FALSE)
}

# Select p-values we are estimating

p.adjust.log <- function(
    p,
    method = c('bonferroni', 'fdr', 'BH'),
    n = length(p),
    p.log = TRUE,
    log.base = exp(1)) {
  if (!p.log) return(p.adjust(p, method))

  method <- match.arg(method)
  if (method == 'fdr') {
    method <- 'BH'
  }

  if (method == 'BH') {
    lp <- length(p)
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    res <- logb(n, base = log.base) - logb(i, base = log.base) + p[o]
    pmin(1, cummin(res))[ro]
  } else if (method == 'bonferroni') {
    logb(n, base = log.base) + p
  }
}

G <- matrix(
  c(0, 0, 0, 0, 0,
    sqrt(0.3), 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, -1 * sqrt(0.1), 0, 0, 0,
    0, -1 * sqrt(0.1), sqrt(0.2), sqrt(0.25), 0),
  nrow = 5,
  byrow = 5
)

# Remove winner's curse bias
Ztrue <- with(dat, beta_marg/se_beta_hat)
pval_true <- 2*pnorm(-abs(Ztrue))
minp <- apply(pval_true, 1, min)
ix <- which(minp < 5e-8)

## 0. Fit true model

B_true <- dat$direct_trait_effects
B_true[!B_true == 0] <- 1

five_node_g <- graph_from_adjacency_matrix(B_true)

true_model <- with(dat,
                   esmr(
                     beta_hat_X = beta_hat,
                     se_X = s_estimate,
                     variant_ix = ix,
                     G = diag(5), # required for network problem
                     direct_effect_template = B_true))

## 1. esmr to generate super graph
MVMR_models <- lapply(seq_len(nrow(G)), function(i) {
  mvmr_minp <- apply(pval_true[,-i], 1, min)
  mvmr_ix <- which(mvmr_minp < 5e-8)

  # For now use true G, eventual switch to estimating G

  trait_order <- c(i, (1:5)[-i])
  true_G_total <- dat$total_trait_effects[trait_order,trait_order]
  true_G_total[1, ] <- true_G_total[, 1] <- 0
  diag(true_G_total) <- 1

  with(dat,
       esmr(beta_hat_Y = beta_hat[,i],
            se_Y = s_estimate[,i],
            beta_hat_X = beta_hat[,-i],
            se_X = s_estimate[,-i],
            variant_ix = mvmr_ix,
            G = t(true_G_total),
            beta_joint = TRUE)
  )
})

mvmr_beta_df <- do.call(
  'rbind.data.frame',
  lapply(1:5, function(i) {
    x <- MVMR_models[[i]]

    res <- x$beta[c('beta_m', 'beta_s')]
    res$to <- rep(i, 4)
    res$from <- setdiff(1:5, i)
    res
  })
)

## Build adj matrix from MVMR results
mvmr_beta_df$pval <- 2 * pnorm(-abs(mvmr_beta_df$beta_m / mvmr_beta_df$beta_s))
mvmr_beta_df$pval_weights <- pmin(-log10(mvmr_beta_df$pval), 20)

mvmr_beta_edgelist <- mvmr_beta_df %>%
  filter(pval < threshold) %>%
  select(from, to, pval_weights)

discovery_G <- mvmr_G <- graph_from_edgelist(
  as.matrix(select(mvmr_beta_edgelist, -pval_weights)))

print('Discovery graph: ')
print(mvmr_G)

# ## 1.5: Check if super graph is acyclic; Otherwise 2.
if (! is_acyclic(mvmr_g)) {
  ## 2. Perform feedback arc set with weights on p-values to reduce super graph to "best" acyclic subgraph
  fas <- feedback_arc_set(discovery_G, weights = mvmr_beta_edgelist$pval_weights, algo = "exact")
  print('Removing edges: ')
  print(fas)
  discovery_G <- discovery_G - fas
  print('New graph: ')
  print(discovery_G)
}

missed_edges <- setNames(
  as.data.frame(as_edgelist(five_node_g - discovery_G)),
  c('from', 'to'))
extra_edges <- setNames(
  as.data.frame(as_edgelist(discovery_G - five_node_g)),
  c('from', 'to'))
if (nrow(missed_edges) > 0) missed_edges$missing_edge <- TRUE
if (nrow(extra_edges) > 0) extra_edges$extra_edge <- TRUE

## 3. Fit discovery model

## Fit model on discovery set
discovery_adj_mat <- as.matrix(as_adjacency_matrix(discovery_G))

discovery_model <- with(dat,
                        esmr(
                          beta_hat_X = beta_hat,
                          se_X = s_estimate,
                          variant_ix = ix,
                          G = diag(5), # required for network problem
                          direct_effect_template = discovery_adj_mat,
                          max_iter = 300))

discovery_idx <- which(discovery_adj_mat != 0, arr.ind = T)
discovery_log_pvals <- discovery_model$pvals_dm[discovery_idx]

discovery_log10 <- pmin(- discovery_log_pvals / log(10), 20)
# Bonferroni correct on log p-vals
n <- length(discovery_log10)
# Bonferroni correct: log10(n * p) = log(n * p) / log(10) = (log(n) + log(p)) / log(10)
discovery_bf_adjust_log10 <- - p.adjust.log(
  - discovery_log10, method = 'bonferroni', log.base = 10)
discovery_fdr_adjust_log10 <- - p.adjust.log(
  - discovery_log10, method = 'fdr', log.base = 10
)

keep_no_adjust <- discovery_log10 > -log10(threshold)
keep_bonferroni <- discovery_bf_adjust_log10 > -log10(threshold)
keep_fdr <- discovery_fdr_adjust_log10 > -log10(threshold)

results_df <- cbind.data.frame(
  setNames(data.frame(discovery_idx), c('from', 'to')),
  keep_no_adjust,
  keep_bonferroni,
  keep_fdr,
  discovery_log10,
  discovery_bf_adjust_log10,
  discovery_fdr_adjust_log10)

adj_mat <- matrix(0, nrow = nrow(G), ncol = ncol(G))
adj_mat[
  as.matrix(results_df[results_df$keep_no_adjust, c('from', 'to')])
  ] <- 1

no_adj_model <- with(dat,
                     esmr(
                       beta_hat_X = beta_hat,
                       se_X = s_estimate,
                       variant_ix = ix,
                       G = diag(5), # required for network problem
                       direct_effect_template = adj_mat,
                       max_iter = 300))

if (all(results_df$keep_no_adjust == results_df$keep_fdr)) {
  fdr_model <- no_adj_model
} else {
  # Create new templates for nesmr model
  fdr_adj_mat <- matrix(0, nrow = nrow(G), ncol = ncol(G))
  fdr_adj_mat[
    as.matrix(results_df[results_df$keep_fdr, c('from', 'to')])
  ] <- 1

  fdr_model <- with(dat,
                    esmr(
                      beta_hat_X = beta_hat,
                      se_X = s_estimate,
                      variant_ix = ix,
                      G = diag(5), # required for network problem
                      direct_effect_template = fdr_adj_mat,
                      max_iter = 300))
}

if (all(results_df$keep_fdr == results_df$keep_bonferroni)) {
  bonferroni_model <- fdr_model
} else {
  # Create new templates for nesmr model
  bf_adj_mat <- matrix(0, nrow = nrow(G), ncol = ncol(G))
  bf_adj_mat[
    as.matrix(results_df[results_df$keep_bonferroni, c('from', 'to')])
  ] <- 1

  bonferroni_model <- with(dat,
                           esmr(
                             beta_hat_X = beta_hat,
                             se_X = s_estimate,
                             variant_ix = ix,
                             G = diag(5), # required for network problem
                             direct_effect_template = bf_adj_mat,
                             max_iter = 300))
}

# Add additional information about incorrect edges
results_df <- merge(results_df, extra_edges, by = c('from', 'to'), all = TRUE)
results_df <- merge(results_df, missed_edges, by = c('from', 'to'), all = TRUE)

# True backwards selection algorithm..
# 1. start with discovery model
# 2. remove the edge with the biggest p-value if bigger than x (e.g. your bonferroni etc cutoff)
# 3. refit the model missing the removed edge
# 4. repeat until there are no p-values bigger than x
# Just do for BF for now
backward_df <- results_df[, c('from', 'to', 'discovery_log10')]
idx_log10_pval <- which.min(results_df$discovery_log10)
min_log10_pval <- results_df$discovery_bf_adjust_log10[idx_log10_pval]
# Copy discovery matrix
backward_select_adj_mat <- discovery_adj_mat

backward_select_edges <- list()
backward_select_models <- list()
backward_mod_results <- list()
backward_select_pvals <- list()
backward_results <- list()
last_backward_mod <- NULL
iter <- 1
cat('Starting backward selection...\n')
# While we still have large p-values...
while(min_log10_pval <= -log10(threshold) && iter < n) {
  cat('Backward select: ', iter, '\n')
  backward_select_edges[[iter]] <- as.matrix(backward_df[idx_log10_pval, 1:2])

  backward_results[[iter]] <- cbind.data.frame(
    backward_select_edges[[iter]], min_log10_pval
  )

  # Remove highest p-value from adj matrix
  backward_select_adj_mat[
    backward_select_edges[[iter]]
  ] <- 0

  # Fit the new model
  last_backward_mod <- backward_select_models[[iter]] <- with(
    dat,
    esmr(
      beta_hat_X = beta_hat,
      se_X = s_estimate,
      variant_ix = ix,
      G = diag(5), # required for network problem
      direct_effect_template = backward_select_adj_mat,
      max_iter = 300))


  # TODO Keep track of discovery_bf_adjust_log10?
  backward_idx <- which(backward_select_adj_mat != 0, arr.ind = TRUE)
  backward_log_pvals <- backward_select_models[[iter]]$pvals_dm[backward_idx]
  backward_log10 <- backward_select_pvals[[iter]] <- pmin(- backward_log_pvals / log(10), 20)
  backward_df <- backward_mod_results[[iter]] <- data.frame(
    setNames(data.frame(backward_idx), c('from', 'to')),
    backward_log10
  )
  # backward_bf_adjust_log10 <- - p.adjust.log(
  #   - backward_log10, method = 'bonferroni', log.base = 10, n = n)

  idx_log10_pval <- which.min(backward_log10)
  min_log10_pval <- backward_log10[idx_log10_pval]

  merge_backward_df <- backward_df
  names(merge_backward_df) <- c('from', 'to', paste0('backward_log10_', iter))
  results_df <- merge(
    results_df, merge_backward_df,
    by = c('from', 'to'),
    all.x = TRUE
  )

  iter <- iter + 1
}

cat('Selected', iter - 1, 'edges via backward selection...\n')

# Log Likelihoods
loglik_results <- lapply(
  list(true_model, discovery_model, no_adj_model, fdr_model, bonferroni_model),
  logLik.esmr)
names(loglik_results) <- c(
  "true", "discovery", "no_adj", "fdr", "bonferroni")

if (length(backward_select_models) > 0) {
  backward_loglik <- setNames(
    lapply(backward_select_models, logLik.esmr),
    paste0('backward_select_', seq_along(backward_select_models)))

  loglik_results <- c(loglik_results, backward_loglik)
  print(loglik_results)
}
