
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
dscout <- dscquery(dsc.outdir = "/nfs/turbo/sph-jvmorr/NESMR/simulations/five_node_discovery_backward",
                   targets    = c("simulate.loglik_results", "simulate.results_df", "simulate.DSC_TIME", "simulate.backward_select_edges"),
                   ignore.missing.files = TRUE)

loglik_results <- do.call("rbind.fill", lapply(dscout$simulate.loglik_results, as.data.frame))

B_true <- G
B_true[abs(G) > 0] <- 1

incorrect_edges <- cbind(
    setNames(as.data.frame(which(B_true == 0, arr.ind = TRUE)), c('from', 'to')),
    correct_edge = FALSE)

correct_edges <- cbind(
    setNames(as.data.frame(which(B_true > 0, arr.ind = TRUE)), c('from', 'to')),
    correct_edge = TRUE)

pval_adj_methods <- c('keep_no_adjust', 'keep_bonferroni', 'keep_fdr', 'keep_backward')

head(dscout$simulate.results_df)

bse <- lapply(
    seq_along(dscout$simulate.backward_select_edges),
    function(i) {
        x <- dscout$simulate.backward_select_edges[[i]]
        tryCatch({
            removed_edges <- cbind(
                do.call('rbind.data.frame', lapply(x, `colnames<-`, c('from', 'to'))),
                removed_backward = TRUE,
                i = i)
                left_join(dscout$simulate.results_df[[i]], removed_edges, by = c('from', 'to'))
        },
        error = function(e) NULL
        )
    })

# TODO: Check why 60 is NA?
bse <- bse[sapply(bse, Negate(is.null))]

# TODO: Need to add the removed edges from backward select here too
edge_results <- do.call('rbind.data.frame',
    lapply(seq_along(bse), function(i) {
        x <- bse[[i]]
        x <- merge(x, true_edges, by = c("from", "to"), all = TRUE)
        x$correct_edge[is.na(x$correct_edge)] <- FALSE
        x$removed_backward[is.na(x$removed_backward)] <- FALSE
        x$keep_backward <- ! x$removed_backward

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
