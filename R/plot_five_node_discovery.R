
library(ggplot2)
library(dscrutils)
library(dplyr)
library(tidyr)

reticulate::use_condaenv('dsc')

# Load your data
dscout <- dscquery(dsc.outdir = "/nfs/turbo/sph-jvmorr/NESMR/simulations/five_node_discovery",
                   targets    = c("simulate.loglik_results", "simulate.results_df", "simulate.DSC_TIME"),
                   ignore.missing.files = TRUE)

loglik_results <- do.call("rbind.data.frame", lapply(dscout$simulate.loglik_results, as.data.frame))

true_edges <- as.data.frame(list(from = c(2, 4, 5, 5, 5), to = c(1, 2, 2, 3, 4)))
true_edges$correct_edge <- TRUE

pval_adj_methods <- c('keep_no_adjust', 'keep_bonferroni', 'keep_fdr')

# results <- lapply(seq_along(dscout$simulate.results_df), function(i) {
#         x <- dscout$simulate.results_df[[i]]
#         x <- merge(x, true_edges, by = c("from", "to"), all = TRUE)
#         x$correct_edge[is.na(x$correct_edge)] <- FALSE
#         x
#     })

edge_results <- do.call('rbind.data.frame',
    lapply(seq_along(dscout$simulate.results_df), function(i) {
        x <- dscout$simulate.results_df[[i]]
        x <- merge(x, true_edges, by = c("from", "to"), all = TRUE)
        x$correct_edge[is.na(x$correct_edge)] <- FALSE
        fct_edges <- lapply(
            x[, c(pval_adj_methods, 'correct_edge')], factor, levels = c('TRUE', 'FALSE'))

        res <- do.call('rbind.data.frame', lapply(pval_adj_methods, function(m) {
            cbind(
                as.data.frame(table(interaction(fct_edges[[m]], fct_edges$correct_edge))),
                method = m)
        }))

        res$replicate <- i
        res
    })
)

edge_results %>%
    group_by(Var1, method) %>%
    rename(KeepEdge.Actual = Var1) %>%
    mutate(method = gsub('keep_', '', method)) %>%
    summarise(
        median = median(Freq),
        min = min(Freq),
        max = max(Freq),
        mean = mean(Freq)
    ) %>%
    write.csv(file = "simulations/five_node_discovery/results/correct_edges.csv", row.names = FALSE)
