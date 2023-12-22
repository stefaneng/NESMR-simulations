renv::load('/nfs/turbo/sph-jvmorr/NESMR/simulations')
devtools::load_all('/nfs/turbo/sph-jvmorr/NESMR/esmr')
library(GWASBrewer)
library(igraph)
library(dplyr)

threshold <- 0.05

is_acyclic <- function(g) {
  tryCatch(
    !is.null(topo_sort(g)),
    error = function(e) FALSE)
}

G <- matrix(
  c(0, 0, 0, 0, 0,
    sqrt(0.3), 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, -1*sqrt(0.1), 0, 0, 0,
    0, -1*sqrt(0.1), sqrt(0.2), sqrt(0.25), 0),
  nrow = 5,
  byrow = 5
)

h2 <- c(0.5, 0.3, 0.25, 0.4, 0.3)
## simulate summary statistics
data(ld_mat_list)
data(AF)
dat <- sim_mv(
  G = G,
  N = 40000,
  J = 5e5,
  h2 = h2,
  pi = 500/5e5,
  R_LD = ld_mat_list,
  af = AF,
  est_s = TRUE
)

Z <- with(dat, beta_hat/s_estimate);
dat$pval <- 2*pnorm(-abs(Z));
minp <- apply(dat$pval, 1, min)

# ld pruning
dat$ld_list_minp <- sim_ld_prune(dat, R_LD = ld_mat_list, pvalue = minp)

minp <- apply(dat$pval, 1, min)

ix <- dat$ld_list_minp
minp <- apply(dat$pval[ix,], 1, min)
ix1 <- which(minp < 5e-8)

## 0. Fit true model

B_true <- dat$direct_trait_effects
B_true[!B_true == 0] <- 1

five_node_g <- graph_from_adjacency_matrix(
  B_true
)

# which_beta gives the indices of the F matrix that should be estimated
true_model <- with(dat,
                   esmr(
                     beta_hat_X = beta_hat[ix,],
                     se_X = s_estimate[ix,],
                     ix1 = ix1,
                     G = diag(5), # required for network problem
                     direct_effect_template = B_true))

## 1. esmr to generate super graph

MVMR_models <- lapply(seq_len(nrow(G)), function(i) {
  with(dat,
    esmr(beta_hat_Y = beta_hat[ix,i],
    se_Y = s_estimate[ix,i],
    beta_hat_X = beta_hat[ix,-i],
    se_X = s_estimate[ix,-i],
    ix1 = ix1, # the subset of variants esmr will fit with
    augment_G = TRUE,
    beta_joint = TRUE)
    )
})

mvmr_beta_df <- do.call(
  'rbind.data.frame',
  lapply(1:5, function(i) {
    x <- MVMR_models[[i]]

    res <- x$beta[c('beta_m', 'beta_s')]
    res$from <- rep(i, 4)
    res$to <- setdiff(1:5, i)
    res
  })
)

## Build adj matrix from MVMR results
mvmr_beta_df$pval <- 2*pnorm(-abs(mvmr_beta_df$beta_m / mvmr_beta_df$beta_s))
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
  # Missed edges
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
                     beta_hat_X = beta_hat[ix,],
                     se_X = s_estimate[ix,],
                     ix1 = ix1,
                     G = diag(5), # required for network problem
                     direct_effect_template = discovery_adj_mat,
                     max_iter = 300))

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

adj_mat <- as.matrix(as_adjacency_matrix(
  graph_from_edgelist(
    as.matrix(results_df[results_df$keep_no_adjust, c('from', 'to')])
  )
))

no_adj_model <- with(dat,
                  esmr(
                    beta_hat_X = beta_hat[ix,],
                    se_X = s_estimate[ix,],
                    ix1 = ix1,
                    G = diag(5), # required for network problem
                    direct_effect_template = adj_mat,
                    max_iter = 300))

if (all(results_df$keep_no_adjust == results_df$keep_fdr)) {
  fdr_model <- no_adj_model
} else {
  # Create new templates for nesmr model
  fdr_adj_mat <- as.matrix(as_adjacency_matrix(
    graph_from_edgelist(
      as.matrix(results_df[results_df$keep_fdr, c('from', 'to')])
    )
  ))

  fdr_model <- with(dat,
                    esmr(
                      beta_hat_X = beta_hat[ix,],
                      se_X = s_estimate[ix,],
                      ix1 = ix1,
                      G = diag(5), # required for network problem
                      direct_effect_template = fdr_adj_mat,
                      max_iter = 300))
}

if (all(results_df$keep_fdr == results_df$keep_bonferroni)) {
  bonferroni_model <- fdr_model
} else {
  bf_adj_mat <- as.matrix(as_adjacency_matrix(
    graph_from_edgelist(
      as.matrix(results_df[results_df$keep_bonferroni, c('from', 'to')])
    )
  ))

  bonferroni_model <- with(dat,
                    esmr(
                      beta_hat_X = beta_hat[ix,],
                      se_X = s_estimate[ix,],
                      ix1 = ix1,
                      G = diag(5), # required for network problem
                      direct_effect_template = bf_adj_mat,
                      max_iter = 300))
}

# Add addition information about incorrect edges
results_df <- merge(results_df, extra_edges, by = c('from', 'to'), all = TRUE)
results_df <- merge(results_df, missed_edges, by = c('from', 'to'), all = TRUE)

# Log Likelihoods
loglik_results <- lapply(
  list(true_model, discovery_model, no_adj_model, fdr_model, bonferroni_model),
  logLik.esmr)
names(loglik_results) <- c(
  "true", "discovery", "no_adj", "fdr", "bonferroni")
