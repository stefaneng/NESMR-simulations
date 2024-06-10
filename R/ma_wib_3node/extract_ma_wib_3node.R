library(dscrutils)
library(dplyr)
library(tidyr)
library(reticulate)
library(here)

reticulate::use_condaenv('dsc')

dir <- here("simulations/simulation_results/2024-06-10_ma_wib_3node")
dsc_params <- dscquery(dsc.outdir = dir,
                    targets    = c(
                      "simulate.N",
                      "simulate.J",
                      "simulate.pi_J",
                      "simulate.effect_scale_factor",
                      "simulate.h2",
                      "simulate.graph_type",
                      "fit_esmr.eta",
                      "fit_esmr.alpha"
                    ),
                    ignore.missing.files = TRUE)
dsc_params$i <- as.character(seq_len(nrow(dsc_params)))

sim_results <- dscquery(dsc.outdir = dir,
                   targets    = c(
                     "fit_esmr.sim_results",
                     "fit_esmr.n_variants",
                     "fit_esmr.cursed_mod_results",
                     "fit_esmr.true_mod_results",
                     "fit_esmr.ma_mod_results"
                   ),
                   return.type = "list",
                   ignore.missing.files = TRUE)

all_results <- bind_rows(sim_results$fit_esmr.sim_results, .id = "i") %>%
  left_join(dsc_params, by = "i") %>%
  select(-i)
today <- format(Sys.Date(), "%Y-%m-%d")
saveRDS(all_results, file = print(file.path(dir, paste0(today, "_ma_wib_3node.rds"))))