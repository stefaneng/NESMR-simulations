library(dscrutils)
library(dplyr)
library(tidyr)
reticulate::use_condaenv('dsc')

dscout <- dscquery(dsc.outdir = "/nfs/turbo/sph-jvmorr/NESMR/simulations/five_node_GSEM_eval",
                   targets    = c(
                    "nesmr.nesmr_model",
                    "genomic_sem.genomic_sem_fit"),
                   ignore.missing.files = TRUE)

dsc_n <- dscquery(dsc.outdir = "/nfs/turbo/sph-jvmorr/NESMR/simulations/five_node_GSEM_eval",
                   targets    = c("simulate.n"),
                   ignore.missing.files = TRUE)
#dscout.gsem <- dscquery(dsc.outdir = "/nfs/turbo/sph-jvmorr/NESMR/simulations/five_node_GSEM_eval",
#                   targets    = c("genomic_sem.seed", "genomic_sem.genomic_sem_fit"),
#                   ignore.missing.files = TRUE)

replicates <- length(dscout$DSC) / 2
dscout$nesmr.nesmr_model <- dscout$nesmr.nesmr_model[1:replicates]
dscout$genomic_sem.genomic_sem_fit <- dscout$genomic_sem.genomic_sem_fit[(replicates + 1):length(dscout$genomic_sem.genomic_sem_fit)]

G <- matrix(
  c(0, 0, 0, 0, 0,
    sqrt(0.3), 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, -1 * sqrt(0.1), 0, 0, 0,
    0, -1 * sqrt(0.1), sqrt(0.2), sqrt(0.25), 0),
  nrow = 5,
  byrow = 5
)

extract_lower_triangular <- function(mat) {
    if (is.null(colnames(mat))) {
        col_names <- paste0("V", seq_len(ncol(mat)))
    } else {
        col_names <- colnames(mat)
    }

    low_tri_idx <- which(lower.tri(mat), arr.ind = TRUE)
    low_tri_arrow <- paste0(col_names[low_tri_idx[, 1]], '->', col_names[low_tri_idx[, 2]])
    setNames(mat[lower.tri(mat)], low_tri_arrow)
}

nesmr_results <- lapply(seq_along(dscout$nesmr.nesmr_model), function(i) {
    x <- dscout$nesmr.nesmr_model[[i]]
    if (typeof(x) != 'list') return(NULL)
    de <- extract_lower_triangular(x$direct_effects)

    # Convert log_e to log10 p-values
    pvals_log10 <- extract_lower_triangular(x$pvals_dm) / log(10)

    data.frame(
        model = 'NESMR',
        beta = de,
        pvals_log10 = pvals_log10,
        edge = names(de),
        n = dsc_n$simulate.n[i]
    )
})

genomic_sem_results <- lapply(seq_along(dscout$genomic_sem.genomic_sem_fit), function(i) {
    x <- dscout$genomic_sem.genomic_sem_fit[[i]]
    if (typeof(x) != 'list') return(NULL)
    x <- x$results
    edges <- paste0(x$rhs, '->', x$lhs)
    res <- data.frame(
        model = 'Genomic_SEM',
        beta = x$Unstand_Est,
        pvals_log10 = log10(x$p_value),
        edge = edges,
        n = dsc_n$simulate.n[i]
    )
    res[x$rhs != x$lhs, ]
})

true_edges <- data.frame(
    true_beta = extract_lower_triangular(G)
)
true_edges$edge <- rownames(true_edges)

merged_results <- bind_rows(nesmr_results, genomic_sem_results, .id = 'id') %>%
    `rownames<-`(NULL) %>%
    left_join(
        true_edges,
        by = 'edge'
    ) %>%
    replace_na(list(true_beta = 0))

# Write merged_results to file with todays date
write.csv(merged_results, file = "simulations/five_node_GSEM_eval/nesmr_vs_genomic_sem_results.csv", row.names = FALSE)