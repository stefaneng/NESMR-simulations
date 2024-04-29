renv::load('/nfs/turbo/sph-jvmorr/NESMR/simulations')
devtools::load_all('/nfs/turbo/sph-jvmorr/NESMR/esmr')

library(igraph)

# Note: This has in-selection bias
minp <- apply(dat$pval, 1, min)
ix <- dat$ld_list_minp
minp <- apply(dat$pval[ix,], 1, min)
ix1 <- which(minp < 5e-8)

B_true <- dat$direct_trait_effects
B_true[!B_true == 0] <- 1

# five_node_g <- graph_from_adjacency_matrix(B_true)

# which_beta gives the indices of the F matrix that should be estimated
true_model <- with(dat,
                   esmr(
                     beta_hat_X = beta_hat[ix,],
                     se_X = s_estimate[ix,],
                     ix1 = ix1,
                     G = diag(5), # required for network problem
                     direct_effect_template = B_true))

true_ll <- logLik.esmr(true_model)

# Extra edges
extra_edges <- which(lower.tri(B_true) & ! B_true, arr.ind = TRUE)

extra_edge_df <- do.call('rbind.data.frame', apply(extra_edges, 1, function(idx) {
  idx <- t(as.matrix(idx))
  B <- B_true
  B[idx] <- 1
  extra_model <- with(dat,
                      esmr(
                        beta_hat_X = beta_hat[ix,],
                        se_X = s_estimate[ix,],
                        ix1 = ix1,
                        G = diag(5), # required for network problem
                        direct_effect_template = B))
  extra_ll <- logLik.esmr(extra_model)

  lrt_log_pvalue <- pchisq(
    - 2 * (true_ll - extra_ll),
    df = 1,
    lower.tail = FALSE,
    log.p = TRUE
  )

  dm_log_pvalue <- extra_model$pvals_dm[idx]
  data.frame(
    idx,
    lrt_log_pvalue = lrt_log_pvalue,
    dm_log_pvalue = dm_log_pvalue
  )
}))
