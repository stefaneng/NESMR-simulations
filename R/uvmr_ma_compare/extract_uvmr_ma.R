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

dir <- "/nfs/turbo/sph-jvmorr/NESMR/simulations/2024-06-04_uvmr_ma_compare"
dscout1 <- dscquery(dsc.outdir = dir,
                   targets    = c(
                    "simulate.res"
                    ),
                    return.type = "list",
                   ignore.missing.files = TRUE)

all_res <- do.call('rbind.data.frame', lapply(dscout1$simulate.res, as.data.frame))

today <- format(Sys.Date(), "%Y-%m-%d")
saveRDS(all_res, file = print(file.path(dir, paste0(today, "_uvmr_ma_compare.rds"))))
