
renv::load('/nfs/turbo/sph-jvmorr/NESMR/simulations')
devtools::load_all('/nfs/turbo/sph-jvmorr/NESMR/esmr')
library(ggplot2)
library(dscrutils)
library(reshape2)
library(forcats)

reticulate::use_condaenv('dsc')

# Define the output directory
output_mediation_LD_dir <- "/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node_mediation_LD/figures"

# dscout.mediation <- dscquery(dsc.outdir = "/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node_mediation_no_LD_test",
#                    targets    = c("simulate.dm_pvalue", "simulate.lrt_pvalue",
#                    "simulate.true_model", "simulate.incorrect_model"))

dscout.mediation.LD <- dscquery(dsc.outdir = "/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node_mediation_LD",
                   targets    = c("simulate.dm_pvalue", "simulate.lrt_pvalue",
                   "simulate.true_model", "simulate.incorrect_model"))


dscout.mediation.LD$direct_effects <- lapply(dscout.mediation.LD$simulate.incorrect_model, function(x) {
    lapply(x, function(y) {
        lapply(y, function(z) {
            list(beta_1 = z$direct_effects[2,1], beta_2 = z$direct_effects[3,2], beta_3 = z$direct_effects[3,1])
        })
    })
})

pval_thresholds <- c(5e-8, 5e-9, 5e-10)
r2_thresholds <- c(0.1, 0.01, 0.001, 0.0001)
# L1 = pvalue/beta
# L2 = sim number
# L3 = p-value threshold
# L4 = r2 threshold
res <- melt(dscout.mediation.LD[c('simulate.lrt_pvalue', 'simulate.dm_pvalue', 'direct_effects')]) %>%
    rename(
        beta = L5,
        pval_threshold = L4,
        r2_threshold = L3,
        sim_number = L2,
        key = L1
    ) %>%
    mutate(
        pval_threshold = pval_thresholds[pval_threshold],
        r2_threshold = r2_thresholds[r2_threshold]
    )


true_beta <- data.frame(
    beta = c('beta_1', 'beta_2', 'beta_3'),
    true_value = c(sqrt(0.4), sqrt(0.2), 0)
    )

beta_df <- left_join(
    res %>% filter(key == 'direct_effects'),
    true_beta,
    by = 'beta'
    ) %>%
    mutate(
        bias = value - true_value,
        beta = fct_recode(
            factor(beta),
            `beta["2,1"] == sqrt(0.4)` = "beta_1",
            `beta["3,2"] == sqrt(0.2)` = "beta_2",
            `beta["3,1"] == 0` = "beta_3"
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
    ggtitle("3 Node with LD varying r2 and p-value thresholds")

ggsave(
    filename = print(
        file.path(output_mediation_LD_dir, sprintf('%s_beta_boxplot_3node_pval_r2_thresholds.jpeg', Sys.Date()))),
    plot = beta_boxplot,
    units = "in", width = 12, height = 10
    )
