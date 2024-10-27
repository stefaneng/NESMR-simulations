library(dscrutils)
library(dplyr)
library(tidyr)
library(reticulate)

reticulate::use_condaenv('dsc')

dir <- "/nfs/turbo/sph-jvmorr/NESMR/simulations/simulation_results/2024-06-11_ma_wib_3node"
dsc_params <- dscquery(dsc.outdir = dir,
                    targets    = c(
                      "simulate.N",
                      "simulate.J",
                      "simulate.pi_J",
                      "simulate.effect_scale_factor",
                      "simulate.h2",
                      "simulate.graph_type",
                      "fit_esmr.alpha"
                    ),
                    ignore.missing.files = TRUE)
dsc_params$i <- as.character(seq_len(nrow(dsc_params)))

sim_results <- dscquery(dsc.outdir = dir,
                   targets    = c(
                     "fit_esmr.sim_results"
                   ),
                   return.type = "list",
                   ignore.missing.files = TRUE)

# error_results <- which(! unlist(lapply(sim_results$fit_esmr.sim_results, is.data.frame)))

#all_results <- bind_rows(sim_results$fit_esmr.sim_results[[error_results]], .id = "i") %>%
all_results <- bind_rows(sim_results$fit_esmr.sim_results, .id = "i") %>%
  right_join(dsc_params, by = "i") %>%
  select(-i)
today <- format(Sys.Date(), "%Y-%m-%d")
output_path <- print(file.path(dir, paste0(today, "_ma_wib_3node_no_errors.rds")))
saveRDS(all_results, file = output_path)