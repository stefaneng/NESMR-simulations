
library(ggplot2)
library(dscrutils)

# reticulate::use_condaenv('dsc')

dscout <- dscquery(dsc.outdir = "/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node_output",
                   targets    = c("simulate.lrt_pvalue", "simulate.pval_dm"),
                   ignore.missing.files = TRUE)

sim_pvalues <- dscout$simulate.lrt_pvalue
dm_pvalues <- exp(dscout$simulate.pval_dm)
pvals_df <- data.frame(lrt_pvalues = sim_pvalues, dm_pvalues = dm_pvalues)

pvalue_hist <- ggplot(pvals_df) +
    geom_histogram(aes(x = lrt_pvalues, after_stat(density)), breaks = seq(0, 1, by = 0.1)) +
    theme_minimal() +
    xlab("LRT P-values") +
    theme(text = element_text(size = 24))

ggsave(
    filename = sprintf("/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node/results/%s_lrt_pvalue_hist_dsc.jpg", Sys.Date()),
    units = "in", width = 5, height = 5,
    plot = pvalue_hist
    )

pvalue_dm_hist <- ggplot(pvals_df) +
    geom_histogram(aes(x = dm_pvalues, after_stat(density)), breaks = seq(0, 1, by = 0.1)) +
    theme_minimal() +
    xlab("Delta Method P-values") +
    theme(text = element_text(size = 24))

ggsave(
    filename = print(sprintf("/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node/results/%s_delta_method_pvalue_hist_dsc.jpg", Sys.Date())),
    units = "in", width = 5, height = 5,
    plot = pvalue_dm_hist
    )

lrt_pvalue_qqplot <- ggplot(pvals_df, aes(sample = lrt_pvalues)) +
    geom_qq(distribution = stats::qunif) +
    geom_qq_line(color = "darkgrey", distribution = stats::qunif) +
    coord_fixed(ratio = 1) +
    xlab("theoretical") +
    ylab("sample") +
    theme_minimal() +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
    theme(text = element_text(size = 24))

ggsave(
    filename = sprintf("/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node/results/%s_lrt_pvalue_qqplot_dsc.jpg", Sys.Date()),
    plot = lrt_pvalue_qqplot,
    units = "in", width = 5, height = 5
)

# https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
qlog10unif <- function(p) {
    -log10(qunif(1 - p))
}

lrt_pvalue_qqplot_log <- ggplot(pvals_df, aes(sample = -log10(lrt_pvalues))) +
    geom_qq(distribution = qlog10unif) +
    geom_qq_line(color = "darkgrey", distribution = qlog10unif) +
    xlab(expression(paste("Expected (",-log[10], " p-value)"))) +
    ylab(expression(paste("LRT (",-log[10], " p-value)"))) +
    theme_minimal() +
    theme(text = element_text(size = 24))

ggsave(
    filename = sprintf("/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node/results/%s_lrt_pvalue_qqplot_log_dsc.jpg", Sys.Date()),
    plot = lrt_pvalue_qqplot_log,
    units = "in", width = 5, height = 5
)

dm_pvalue_qqplot_log <- ggplot(pvals_df, aes(sample = -log10(dm_pvalues))) +
    geom_qq(distribution = qlog10unif) +
    geom_qq_line(color = "darkgrey", distribution = qlog10unif) +
    xlab(expression(atop("Expected", paste("(", -log[10], " p-value)")))) +
    ylab(expression(atop("Delta Method", paste("(", -log[10], " p-value)")))) +
    theme_minimal() +
    theme(text = element_text(size = 24))

ggsave(
    filename = print(sprintf("/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node/results/%s_delta_method_pvalue_qqplot_log_dsc.jpg", Sys.Date())),
    plot = dm_pvalue_qqplot_log,
    units = "in", width = 5, height = 5
)

## Compare LRT and delta method

scatter_pvals_log <- ggplot(pvals_df, aes(x = -log10(dm_pvalues), y = -log10(lrt_pvalues))) +
    geom_point() +
    theme_minimal() +
    xlab(expression(atop("Delta Method", paste("(", -log[10], " p-value)")))) +
    ylab(expression(atop("LRT", paste("(", -log[10], " p-value)")))) +
    geom_abline(slope = 1, color = "red") +
    theme(text = element_text(size = 24))

 ggsave(
    filename = print(sprintf("/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node/results/%s_scatter_pvalues_dm_lrt_log10.jpg", Sys.Date())),
    plot = scatter_pvals_log,
    units = "in", width = 5, height = 5
)

scatter_pvals <- ggplot(pvals_df, aes(x = dm_pvalues, y = lrt_pvalues)) +
    geom_point() +
    theme_minimal() +
    xlab("Delta Method p-value") +
    ylab("LRT p-value") +
    geom_abline(slope = 1, color = "red") +
    theme(text = element_text(size = 24))

 ggsave(
    filename = print(sprintf("/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node/results/%s_scatter_pvalues_dm_lrt.jpg", Sys.Date())),
    plot = scatter_pvals,
    units = "in", width = 5, height = 5
)

mean(sim_pvalues < 0.05, na.rm =TRUE)
