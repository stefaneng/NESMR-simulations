library(dscrutils)
library(dplyr)
library(tidyr)
reticulate::use_condaenv('dsc')

dir <- "/nfs/turbo/sph-jvmorr/NESMR/simulations/2024-06-05_uvmr_ma"
dscout1 <- dscquery(dsc.outdir = dir,
                   targets    = c(
                    "simulate.res"
                    ),
                    return.type = "list",
                   ignore.missing.files = TRUE)

dscout_rest <- dscquery(dsc.outdir = dir,
                   targets    = c(
                    "simulate.n",
                    "simulate.J",
                    "simulate.pi_J",
                    "simulate.h2",
                    "simulate.true_beta"
                    ),
                   ignore.missing.files = TRUE)


all_res <- do.call('rbind.data.frame', lapply(dscout1$simulate.res, as.data.frame))

all_res <- cbind(dscout_rest, all_res)

today <- format(Sys.Date(), "%Y-%m-%d")
saveRDS(all_res, file = print(file.path(dir, paste0(today, "_uvmr_ma_compare.rds"))))
