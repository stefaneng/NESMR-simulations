library(dscrutils)
library(dplyr)
library(tidyr)
library(NESMR.Sims)
reticulate::use_condaenv('dsc')

dir <- "/nfs/turbo/sph-jvmorr/NESMR/simulations/simulation_results/2024-10-22_4node_logdet_permn"
dscout1 <- dscquery(dsc.outdir = dir,
                   targets    = c(
                    "process.summary_sim_results",
                    "simulate.effect_modifier"
                    ),
                   ignore.missing.files = TRUE)

sim_data_join <- bind_rows(dscout1$process.summary_sim_result, .id = "sim_id") %>%
    left_join(
        data.frame(sim_id = as.character(seq_along(dscout1$simulate.effect_modifier)), effect_modifier = dscout1$simulate.effect_modifier)
    )

dscout2 <- dscquery(dsc.outdir = dir,
                   targets    = c(
                    "simulate.effect_modifier",
                    "simulate.N",
                    "simulate.J",
                    "simulate.pi_J",
                    "simulate.h2",
                    "simulate.alpha"
                    ),
                   ignore.missing.files = TRUE)

dscout2$sim_id <- as.character(seq_along(dscout2$DSC))

sim_data_join <- left_join(
    sim_data_join,
    dscout2
)

save_results(
    sim_data_join, "all_permn_4node_sim", ext = "rds", output_dir = "/nfs/turbo/sph-jvmorr/NESMR/simulations/simulation_results/2024-10-22_4node_logdet_permn/final_outputs")

sim_data_join %>%
  group_by(sim_id, effect_modifier) %>%
  arrange(AIC) %>%
  summarize(
    aic_rank_correct = which(correct_config)
  ) %>%
  View()

View(head(sim_data_join))

# Compute the "coverage" of the likelihood space

sim_data_join %>%
    group_by(sim_id, effect_modifier) %>%
    mutate(
        aic_rank = rank(AIC)
    ) %>%
    summarize(
        graph_search_coverage = sum(posterior_probs * !is.na(backselect_post_probs)),
        aic_rank_correct = max(correct_config * aic_rank),
        aic_correct = min(AIC * correct_config),
        min_aic = min(AIC)
    ) %>%
    mutate(
        aic_diff = aic_correct - min_aic
    ) %>%
    View()


sim_data_join$backselect_post_probs
