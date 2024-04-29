library(dplyr)
library(tidyr)
library(ggplot2)

results_dir <- "/nfs/turbo/sph-jvmorr/NESMR/simulations/five_node_GSEM_eval"
merged_results <- read.csv(file.path(results_dir, "nesmr_vs_genomic_sem_results_2.csv")) %>%
    mutate(
        bias = beta - true_beta,
        true_null = abs(true_beta) < 1e-8
    )

merged_results %>%
    group_by(from, to, model, n) %>%
    summarise(
        bias = mean(bias),
        max_bias = max(bias),
        min = min(bias)
    ) %>%
    arrange(desc(abs(max_bias)))

true_edges <- select(merged_results, from, to, true_beta) %>%
    distinct()

merged_results %>%
    group_by(model, n) %>%
    group_map(~{
        # TODO: Filter out the most extreme points
        # Get mean bias per edge
        mean_bias <- .x %>%
            group_by(from, to) %>%
            summarise(
                mean_bias = mean(bias)
            )
        p <- ggplot(.) +
            facet_grid(rows = vars(from), cols = vars(to),
                labeller = labeller(from = label_both, to = label_both))
        if (.y$model == "Genomic_SEM") {
            p <- p + xlim(c(-1, 1))
        }
        p <- p +
            geom_histogram(aes(x = bias, fill = true_null), alpha = 0.5, position="identity", na.rm = T) +
            scale_fill_manual(name = "True Null", values = c("TRUE" = "grey", "FALSE" = "orange")) +
            geom_vline(data = true_edges, aes(xintercept = 0), color = "red", lty = 2) +
            geom_vline(data = mean_bias, aes(xintercept = mean_bias), color = "black") +
            ggtitle(sprintf('N = %s, Model = %s', .y$n, .y$model)) +
            theme_minimal(base_size = 24) +
            theme(
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

        ggsave(print(file.path(results_dir, sprintf("hist_bias_%s_%s.pdf", .y$n, .y$model))),
            p, width = 14, height = 14, dpi = 300)
    })
