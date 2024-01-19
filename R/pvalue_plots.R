library(ggplot2)
library(ggpubr)

# Function to create and save plots
create_and_save_plots <- function(pvals_df, pvalue_col, method_name = c('delta_method', 'lrt'), suffix = '', output_dir, main = "") {
    method_name <- match.arg(method_name)
    nice_name <- if (method_name == "delta_method") "Delta Method" else "LRT"
    if (! startsWith(suffix, "_")) {
        suffix <- paste0("_", suffix)
    }
    # Histogram
    pvalue_hist <- ggplot(pvals_df) +
        geom_histogram(aes(x = get(pvalue_col), after_stat(density)), breaks = seq(0, 1, by = 0.1)) +
        theme_minimal() +
        xlab(sprintf("%s P-values", nice_name)) +
        theme(text = element_text(size = 24), plot.title = element_text(size = 20)) +
        ggtitle(main)

    pvals <- pvals_df[[pvalue_col]]

    alpha <- 0.05
    lower_bound <- qbeta(alpha / 2, seq_along(pvals), rev(seq_along(pvals)))
    upper_bound <- qbeta(1 - alpha / 2, seq_along(pvals), rev(seq_along(pvals)))

    data <- data.frame(
        pvals = pvals
    )

    pvalue_col <- 'pvals'

    pvalue_qqplot_log <- data %>%
        arrange(pvals) %>%
        mutate(
            observed = -log10(pvals),
            expected = -log10(ppoints(nrow(data))),
            lower_bound_log = -log10(lower_bound),
            upper_bound_log = -log10(upper_bound)) %>%
    ggplot(aes(x = expected, y = observed)) +
        geom_point() +
        geom_ribbon(aes(ymin = lower_bound_log, ymax = upper_bound_log), alpha = 0.2) +
        geom_abline(slope = 1, intercept = 0, color = "red") +
        xlab(expression(atop("Expected", paste("(", -log[10], " p-value)")))) +
        ylab(bquote(atop(.(nice_name), paste("(", -log[10], " p-value)")))) +
        theme_minimal() +
        theme(text = element_text(size = 24), plot.title = element_text(size = 20)) +
        ggtitle(main)

    final_plot <- ggarrange(pvalue_hist, pvalue_qqplot_log, ncol = 2)

    ggsave(
        filename = file.path(
            output_dir,
            sprintf("%s_%s_pvalue_hist_qqplot_log_dsc%s.jpg", Sys.Date(), method_name, suffix)),
        plot = final_plot,
        units = "in", width = 12, height = 6
    )
}