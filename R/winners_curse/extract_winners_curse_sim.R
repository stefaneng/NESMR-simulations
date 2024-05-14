library(dscrutils)
library(dplyr)
library(tidyr)
reticulate::use_condaenv('dsc')

extract_lower_triangular <- function(mat, prefix = "") {
    if (is.null(colnames(mat))) {
        col_names <- paste0("V", seq_len(ncol(mat)))
    } else {
        col_names <- colnames(mat)
    }

    low_tri_idx <- which(lower.tri(mat), arr.ind = TRUE)
    low_tri_arrow <- paste0(prefix, col_names[low_tri_idx[, 1]], '_', col_names[low_tri_idx[, 2]])
    setNames(mat[lower.tri(mat)], low_tri_arrow)
}

dir <- "/nfs/turbo/sph-jvmorr/NESMR/simulations/2024-05-13_winners_curse"
dscout1 <- dscquery(dsc.outdir = dir,
                   targets    = c(
                    "simulate.scale_effects",
                    "winners_curse_nesmr.cursed_model",
                    "winners_curse_nesmr.cursed_ix",
                    "winners_curse_nesmr.ma_ix",
                    "winners_curse_nesmr.ma_SNP_model",
                    "winners_curse_nesmr.eta",
                    "winners_curse_nesmr.alpha"
                    ),
                    return.type = "list",
                   ignore.missing.files = TRUE)

valid_sims <- sapply(seq_along(dscout1$cursed_model), function(i) {
    return(class(dscout1$cursed_model[[i]]) == "list" &&
        class(dscout1$ma_SNP_model[[i]]) == "list")
})

res <- lapply(seq_along(dscout1$DSC)[valid_sims], function(i) {
    cursed_mod <- dscout1$cursed_model[[i]]
    ma_mod <- dscout1$ma_SNP_model[[i]]
    cursed_n_variants <- length(dscout1$cursed_ix[[i]])
    ma_n_variants <- length(dscout1$ma_ix[[i]])


    data.frame(
        t(extract_lower_triangular(cursed_mod$direct_effects, prefix = "cursed_beta_")),
        t(extract_lower_triangular(ma_mod$direct_effects, prefix = "ma_beta_")),
        t(extract_lower_triangular(-cursed_mod$pvals_dm / log(10), prefix = "cursed_log10_pvals_")),
        t(extract_lower_triangular(-ma_mod$pvals_dm / log(10), prefix = "ma_log10_pvals_")),
        cursed_n_variants,
        ma_n_variants,
        eta = dscout1$eta[[i]],
        alpha = dscout1$alpha[[i]],
        scale_effects = dscout1$scale_effects[[i]],
        replicate = dscout1$DSC[i]
    )
})

res_merged <- do.call(rbind, res)

today <- format(Sys.Date(), "%Y-%m-%d")
saveRDS(res_merged, file = print(file.path(dir, paste0(today, "_winners_curse_results.rds"))))