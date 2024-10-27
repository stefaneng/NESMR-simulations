
renv::load('/nfs/turbo/sph-jvmorr/NESMR/simulations')
# devtools::load_all('/nfs/turbo/sph-jvmorr/NESMR/esmr')
library(ggplot2)
library(dscrutils)
library(reshape2)
library(forcats)

# source("/nfs/turbo/sph-jvmorr/NESMR/simulations/R/pvalue_plots.R")

reticulate::use_condaenv('dsc')

# Define the output directory
output_collider_dir <- "/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node_collider/figures"

dscout.collider.LD <- dscquery(dsc.outdir = "/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node_collider",
                   targets    = c("simulate.dm_pvalue", "simulate.lrt_pvalue",
                   "simulate.true_model", "simulate.incorrect_model"))

dscout.collider.LD$direct_effects <- lapply(dscout.collider.LD$simulate.incorrect_model, function(x) {
    lapply(x, function(z) {
            list(beta_1 = z$direct_effects[2,1], beta_2 = z$direct_effects[3,1], beta_3 = z$direct_effects[3,2])
        })
})

# pval_thresholds <- c(5e-8, 5e-9, 5e-10)
r2_thresholds <- c(0.1, 0.01, 0.001, 0.0001)

collider.melt <- melt(dscout.collider.LD[c('simulate.lrt_pvalue', 'simulate.dm_pvalue', 'direct_effects')]) %>%
    rename(
        beta = L4,
        r2_threshold = L3,
        sim_number = L2,
        key = L1
    ) %>%
    mutate(
        r2_threshold = r2_thresholds[r2_threshold]
    )

collider.melt$pval_threshold <- 5e-8

true_beta <- data.frame(
    beta = c('beta_1', 'beta_2', 'beta_3'),
    true_value = c(sqrt(0.4), sqrt(0.2), 0)
    )

beta_df <- left_join(
    collider.melt %>% filter(key == 'direct_effects'),
    true_beta,
    by = 'beta'
    ) %>%
    mutate(
        bias = value - true_value,
        beta = fct_recode(
            factor(beta),
            `beta["2,1"] == sqrt(0.4)` = "beta_1",
            `beta["3,1"] == sqrt(0.2)` = "beta_2",
            `beta["3,2"] == 0` = "beta_3"
            ),
        r2_threshold = fct_rev(as.factor(r2_threshold))
    )

# levels(beta_df$r2_threshold) <-
r2_labeller <- sapply(levels(beta_df$r2_threshold), function(x) as.expression(bquote(r^2 < .(as.numeric(x)))))

levels(beta_df$r2_threshold) <- r2_labeller

beta_boxplot <- beta_df %>%
    mutate(
        pval_threshold = fct_rev(as.factor(pval_threshold)),
        ) %>%
ggplot(aes(x = pval_threshold, y = bias)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, lty = 2, col = 'red') +
    facet_grid(cols = vars(r2_threshold), rows = vars(beta), labeller = label_parsed) +
    xlab("p-value threshold") +
    ylab(expression(hat(beta) - beta)) +
    theme_minimal() +
    theme(text = element_text(size = 24), plot.title = element_text(size = 24), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle("3 Node Collider with LD varying r2 thresholds")

ggsave(
    filename = print(
        file.path(output_collider_dir, sprintf('%s_beta_boxplot_3node_r2_thresholds_collider.jpeg', Sys.Date()))),
    plot = beta_boxplot,
    units = "in", width = 12, height = 14
    )
