
library(ggplot2)
library(dplyr)
library(dscrutils)

reticulate::use_condaenv('dsc')

# https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
qlog10unif <- function(p) {
    -log10(qunif(1 - p))
}

# Function to create and save plots
create_and_save_plots <- function(
    pvals_df, pvalue_col, method_name = c('delta_method', 'lrt'), key = "", output_dir
    ) {
    method_name <- match.arg(method_name)
    nice_name <- if (method_name == "delta_method") "Delta Method" else "LRT"
    file_key <- paste0(gsub(' ', '_', key), collapse = '_')
    title_key <- paste0(gsub(' ', '_', key), collapse = '-')
    # Histogram
    pvalue_hist <- ggplot(pvals_df) +
        geom_histogram(aes(x = get(pvalue_col), after_stat(density)), breaks = seq(0, 1, by = 0.1)) +
        theme_minimal() +
        ggtitle(title_key) +
        xlab(sprintf("%s P-values", nice_name)) +
        theme(text = element_text(size = 24))

    ggsave(
        filename = file.path(
            output_dir,
            sprintf("%s_%s_%s_pvalue_hist_dsc.jpg", Sys.Date(), method_name, file_key)),
        units = "in", width = 5, height = 5,
        plot = pvalue_hist
    )

    # QQ Plot Log
    pvalue_qqplot_log <- ggplot(pvals_df, aes(sample = -log10(get(pvalue_col)))) +
        geom_qq(distribution = qlog10unif) +
        geom_qq_line(color = "darkgrey", distribution = qlog10unif) +
        xlab(expression(atop("Expected", paste("(", -log[10], " p-value)")))) +
        ylab(bquote(atop(.(nice_name), paste("(", -log[10], " p-value)")))) +
        theme_minimal() +
        ggtitle(title_key) +
        theme(text = element_text(size = 24))

    ggsave(
        filename = file.path(
            output_dir,
            sprintf("%s_%s_%s_pvalue_qqplot_log_dsc.jpg", Sys.Date(), method_name, file_key)),
        plot = pvalue_qqplot_log,
        units = "in", width = 5, height = 5
    )
}

# Define the output directory
output_dir <- "/nfs/turbo/sph-jvmorr/NESMR/simulations/five_node_five_one/figures"

# Load your data
dscout <- dscquery(dsc.outdir = "/nfs/turbo/sph-jvmorr/NESMR/simulations/five_node_five_one",
                   targets    = c("one_edge_test.extra_edge_df"),
                   ignore.missing.files = TRUE)

extra_edge_df <- bind_rows(dscout$one_edge_test.extra_edge_df, .id = "id")

# Calibration
extra_edge_df %>%
    mutate(
        lrt_pvalue = exp(lrt_log_pvalue),
        dm_pvalue = exp(dm_log_pvalue)
    ) %>%
    group_by(row, col) %>%
    summarise(
        lrt_cal = mean(lrt_pvalue < 0.05),
        dm_cal = mean(dm_pvalue < 0.05)
    )

extra_edge_df %>%
    mutate(
        lrt_pvalue = exp(lrt_log_pvalue),
        dm_pvalue = exp(dm_log_pvalue)
    ) %>%
    group_by(row, col) %>%
    group_walk(~ create_and_save_plots(
        .x, pvalue_col = "lrt_pvalue", method_name = "lrt", key = .y, output_dir)
        ) %>%
    group_walk(~ create_and_save_plots(
        .x, pvalue_col = "dm_pvalue", method_name = "delta_method", key = .y, output_dir)
        )


# What is going on here?

stop()

lrt_log10_pvalue <- -log10(lrt_pvalue)

sim_pvalues <- dscout$simulate.lrt_pvalue
dm_pvalues <- exp(dscout$simulate.pval_dm)
pvals_df <- data.frame(lrt_pvalues = sim_pvalues, dm_pvalues = dm_pvalues)

# Call the function for LRT P-values
create_and_save_plots(pvals_df, "lrt_pvalues", "lrt", output_dir)

# Call the function for Delta Method P-values
create_and_save_plots(pvals_df, "dm_pvalues", "delta_method", output_dir)

# Scatter plots of p-values

scatter_pvals <- ggplot(pvals_df, aes(x = dm_pvalues, y = lrt_pvalues)) +
    geom_point() +
    theme_minimal() +
    xlab("Delta Method p-value") +
    ylab("LRT p-value") +
    geom_abline(slope = 1, color = "red") +
    theme(text = element_text(size = 24))

ggsave(
    filename = file.path(
        output_dir,
        sprintf("%s_scatter_pvalues_dm_lrt.jpg", Sys.Date()
            )
    ),
    plot = scatter_pvals,
    units = "in",
    width = 5, height = 5
)

cat('LRT p-values < 0.05', mean(dscout$lrt_pvalue  < 0.05, na.rm =TRUE), '\n')
cat('Delta method p-values < 0.05', mean(dscout$dm_pvalue < 0.05, na.rm =TRUE), '\n')
