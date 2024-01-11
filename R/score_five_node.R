
library(ggplot2)
library(dscrutils)
library(tidyr)
library(plyr)
library(dplyr)

reticulate::use_condaenv('dsc')

G <- matrix(
  c(0, 0, 0, 0, 0,
    sqrt(0.3), 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, -1*sqrt(0.1), 0, 0, 0,
    0, -1*sqrt(0.1), sqrt(0.2), sqrt(0.25), 0),
  nrow = 5,
  byrow = 5
)

# Load your data
dscout <- dscquery(dsc.outdir = "/nfs/turbo/sph-jvmorr/NESMR/simulations/five_node_discovery",
                   targets    = c(
                    "discovery_algo.loglik_results", "discovery_algo.results_df", "discovery_algo.DSC_TIME",
                    "discovery_algo.backward_select_edges", "discovery_algo.discovery_model",
                    "discovery_algo.last_backward_mod", "discovery_algo.backward_select_adj_mat",
                    "discovery_algo.backward_results", "discovery_algo.backward_select_pvals",
                    "discovery_algo.backward_select_edges", "discovery_algo.backward_mod_results",
                    "discovery_algo.results_df"
                    ))

#

threshold <- 0.05


results_df <- dscout$discovery_algo.results_df

five_one_pvals_discovery <- lapply(results_df, function(x) {
  10^(-x$discovery_log10[x$from == 5 & x$to == 1])
})


# lapply(results_df, function(x) {
#   x[x$from == 5 & x$to == 1, ]
# })
# head(results_df)


# backward_select_df <- lapply(dscout$discovery_algo.backward_mod_results, function(x) {
#   do.call('rbind.data.frame', lapply(seq_along(x), function(i) {
#     y <- x[[i]]
#     res <- y %>%
#       filter(backward_log10 <= -log10(threshold)) %>%
#       slice(which.min(backward_log10))
#     if (nrow(res) > 0) cbind(res, select_idx = i)
#     else res
#     }))
#   })

## Discovery model edges here: dscout$discovery_algo.results_df
# results_backward <- lapply(seq_along(dscout$discovery_algo.results_df), function(i) {
#   x <- dscout$discovery_algo.results_df[[i]]
#   y <- backward_select_df[[i]]
#   merge(x, y, by = c('from', 'to'), all.x = TRUE)
# })

# results_backward <- lapply(dscout$discovery_algo.backward_select_edges, function(x) {
#   do.call('rbind.data.frame', lapply(seq_along(x), function(i) {
#     y <- x[[i]]
#     cbind(y, backward_idx = i)
#   }))
# })

# merged_results_backward <- lapply(seq_along(dscout$discovery_algo.results_df), function(i) {
#   x <- dscout$discovery_algo.results_df[[i]]
#   y <- results_backward[[i]]
#   merge(x, y, by = c('from', 'to'), all.x = TRUE)
# })

# First get extra edges for
n <- length(dscout$discovery_algo.results_df)
# discovery_extra_edges <- rep(0, n)
discovery_extra_edges <- list()
backward_extra_edges <- list()
results_back_select <- dscout$discovery_algo.results_df
for (i in seq_len(n)) {
  back_select_edges <- dscout$discovery_algo.backward_select_edges[[i]]
  merge_backward_df <- do.call('rbind.data.frame', lapply(seq_along(back_select_edges), function(j) {
    y <- back_select_edges[[j]]
    cbind(y, back_select_idx = j)
  }))
  results_back_select[[i]] <- merge(
    results_back_select[[i]],
    merge_backward_df,
    by = c('from', 'to'),
    all.x = TRUE
  )
  results_back_select[[i]]$keep_back_select <- is.na(results_back_select[[i]]$back_select_idx)

  discovery_extra_edges[[i]] <- results_back_select[[i]] %>%
    filter(extra_edge) %>%
    mutate(edge = sprintf('%s->%s', from, to)) %>%
    pull(edge) %>%
    as.list()

  backward_extra_edges[[i]] <- results_back_select[[i]] %>%
    filter(extra_edge & keep_back_select) %>%
    mutate(edge = sprintf('%s->%s', from, to)) %>%
    pull(edge) %>%
    as.list()


}

# Overall extra edges in discovery
as.data.frame(table(unlist(discovery_extra_edges)))

as.data.frame(table(unlist(backward_extra_edges)))
