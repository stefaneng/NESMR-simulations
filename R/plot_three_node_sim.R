
renv::load('/nfs/turbo/sph-jvmorr/NESMR/simulations')
devtools::load_all('/nfs/turbo/sph-jvmorr/NESMR/esmr')
library(ggplot2)
library(dscrutils)
library(forcats)

# source("/nfs/turbo/sph-jvmorr/NESMR/simulations/R/pvalue_plots.R")

reticulate::use_condaenv('dsc')

# Define the output directory
output_collider_dir <- "/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node_collider/figures"
output_mediation_dir <- "/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node_mediation_no_LD_test/figures"
output_mediation_LD_dir <- "/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node_mediation_LD/figures"

dscout.collider.LD <- dscquery(dsc.outdir = "/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node_collider",
                   targets    = c("simulate.dm_pvalue", "simulate.lrt_pvalue",
                   "simulate.true_model", "simulate.incorrect_model"))

dscout.mediation <- dscquery(dsc.outdir = "/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node_mediation_no_LD_test",
                   targets    = c("simulate.dm_pvalue", "simulate.lrt_pvalue",
                   "simulate.true_model", "simulate.incorrect_model"))


dscout.mediation.LD <- dscquery(dsc.outdir = "/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node_mediation_LD",
                   targets    = c("simulate.dm_pvalue", "simulate.lrt_pvalue",
                   "simulate.true_model", "simulate.incorrect_model"))

r2_thresholds <- c(0.1, 0.01, 0.001, 0.0001)
pval_thresholds <- c(5e-8, 5e-9, 5e-10)

mediation_LD_pvals_df <- rbind.data.frame(
    bind_rows(lapply(dscout.mediation.LD$simulate.dm_pvalue, function(x) {
        do.call('rbind.data.frame', lapply(seq_along(x), function(i) {
            y <- x[[i]]
            data.frame(
                pvalue = unlist(y),
                pval_threshold = pval_thresholds,
                r2_threshold = r2_thresholds[[i]],
                pvalue_type = 'dm')
        }))
    })),
    bind_rows(lapply(dscout.mediation.LD$simulate.lrt_pvalue, function(x) {
            do.call('rbind.data.frame', lapply(seq_along(x), function(i) {
                y <- x[[i]]
                data.frame(
                    pvalue = unlist(y),
                    pval_threshold = pval_thresholds,
                    r2_threshold = r2_thresholds[[i]],
                    pvalue_type = 'lrt')
            }))
        }))
    )


collider_pvals_df <- do.call('rbind.data.frame',
        lapply(seq_along(dscout.collider.LD$simulate.dm_pvalue), function(i) {
        dm_pvals <- unlist(dscout.collider.LD$simulate.dm_pvalue[[i]])
        lrt_pvals <- unlist(dscout.collider.LD$simulate.lrt_pvalue[[i]])
        rbind.data.frame(
            data.frame(
                pvalue = unlist(dm_pvals), r2_threshold = r2_thresholds, pvalue_type = 'dm'),
            data.frame(
                pvalue = unlist(lrt_pvals), r2_threshold = r2_thresholds, pvalue_type = 'lrt'
                )
        )
    }))

mediation_pvals_df <- do.call('rbind.data.frame',
        lapply(seq_along(dscout.mediation$simulate.dm_pvalue), function(i) {
        dm_pvals <- unlist(dscout.mediation$simulate.dm_pvalue[[i]])
        lrt_pvals <- unlist(dscout.mediation$simulate.lrt_pvalue[[i]])
        rbind.data.frame(
            data.frame(
                pvalue = unlist(dm_pvals), pval_threshold = pval_thresholds, pvalue_type = 'dm',
                beta_var = c(rep("beta_hat", 3), rep("beta_marg", 3))),
            data.frame(
                pvalue = unlist(lrt_pvals), pval_threshold = pval_thresholds, pvalue_type = 'lrt',
                beta_var = c(rep("beta_hat", 3), rep("beta_marg", 3)))
                )
    }))

n_pvals_collider <- nrow(collider_pvals_df) / (length(r2_thresholds) * 2)
alpha <- 0.05
lower_bound_collider <- qbeta(alpha / 2, seq_len(n_pvals_collider), rev(seq_len(n_pvals_collider)))
upper_bound_collider <- qbeta(1 - alpha / 2, seq_len(n_pvals_collider), rev(seq_len(n_pvals_collider)))

grouped_collider_plot <- collider_pvals_df %>%
    group_by(pvalue_type, r2_threshold) %>%
    arrange(pvalue) %>%
    mutate(
        r2_threshold = as.factor(r2_threshold),
        observed = -log10(pvalue),
        expected = -log10(ppoints(n_pvals_collider)),
        lower_bound_log = -log10(lower_bound_collider),
        upper_bound_log = -log10(upper_bound_collider)) %>%
    ggplot() +
    facet_grid(cols =  vars(pvalue_type)) +
        geom_point(aes(x = expected, y = observed, color = r2_threshold)) +
        geom_ribbon(aes(x = expected, y = observed, ymin = lower_bound_log, ymax = upper_bound_log), alpha = 0.2) +
        geom_abline(slope = 1, intercept = 0, color = "red") +
        xlab(expression(atop("Expected", paste("(", -log[10], " p-value)")))) +
        ylab(expression(atop("Observed", paste("(", -log[10], " p-value)")))) +
        theme_minimal() +
        ylim(c(0, 6)) +
        theme(text = element_text(size = 24), plot.title = element_text(size = 24)) +
        ggtitle("3 Node Collider with varying r^2 threshold")

ggsave(
    filename = print(
        file.path(output_collider_dir, sprintf('%s_QQ_collider_3node_r2_thresholds.jpeg', Sys.Date()))),
    plot = grouped_collider_plot,
    units = "in", width = 12, height = 8
    )

n_pvals_mediation <- nrow(mediation_pvals_df) / (2 * 2 * 3)
lower_bound_mediation <- qbeta(alpha / 2, seq_len(n_pvals_mediation), rev(seq_len(n_pvals_mediation)))
upper_bound_mediation <- qbeta(1 - alpha / 2, seq_len(n_pvals_mediation), rev(seq_len(n_pvals_mediation)))

grouped_mediation_plot <- mediation_pvals_df %>%
    group_by(pvalue_type, pval_threshold, beta_var) %>%
    arrange(pvalue) %>%
    mutate(
        pval_threshold = as.factor(pval_threshold),
        observed = -log10(pvalue),
        expected = -log10(ppoints(n_pvals_mediation)),
        lower_bound_log = -log10(lower_bound_mediation),
        upper_bound_log = -log10(upper_bound_mediation)) %>%
    ggplot() +
    facet_grid(cols =  vars(pvalue_type), rows = vars(beta_var)) +
        geom_point(aes(x = expected, y = observed, color = pval_threshold)) +
        geom_ribbon(aes(x = expected, y = observed, ymin = lower_bound_log, ymax = upper_bound_log), alpha = 0.2) +
        geom_abline(slope = 1, intercept = 0, color = "red") +
        xlab(expression(atop("Expected", paste("(", -log[10], " p-value)")))) +
        ylab(expression(atop("Observed", paste("(", -log[10], " p-value)")))) +
        theme_minimal() +
        ylim(c(0, 6)) +
        theme(text = element_text(size = 24), plot.title = element_text(size = 24)) +
        ggtitle("3 Node Mediation No Pleiotropy or LD")

ggsave(
    filename = print(
        file.path(output_mediation_dir, sprintf('%s_QQ_mediation_3node_pval_thresholds.jpeg', Sys.Date()))),
    plot = grouped_mediation_plot,
    units = "in", width = 12, height = 8
    )


## LD r2 threshold
n_pvals_mediation_LD <- nrow(mediation_LD_pvals_df) / (length(pval_thresholds) * length(r2_thresholds) * 2)
lower_bound_mediation_LD <- qbeta(alpha / 2, seq_len(n_pvals_mediation_LD), rev(seq_len(n_pvals_mediation_LD)))
upper_bound_mediation_LD <- qbeta(1 - alpha / 2, seq_len(n_pvals_mediation_LD), rev(seq_len(n_pvals_mediation_LD)))

mediation_LD_pvals_df$r2_threshold_f <- fct_rev(as.factor(mediation_LD_pvals_df$r2_threshold))

grouped_mediation_plot <- mediation_LD_pvals_df %>%
    group_by(pvalue_type, pval_threshold, r2_threshold_f) %>%
    arrange(pvalue) %>%
    mutate(
        pval_threshold = as.factor(pval_threshold),
        observed = -log10(pvalue),
        expected = -log10(ppoints(n_pvals_mediation_LD)),
        lower_bound_log = -log10(lower_bound_mediation_LD),
        upper_bound_log = -log10(upper_bound_mediation_LD)) %>%
    ggplot() +
    facet_grid(
        cols = vars(pvalue_type), rows =  vars(r2_threshold_f),
        labeller = labeller(.rows = function(x) paste0('r2 < ', x))) +
        geom_point(aes(x = expected, y = observed, color = pval_threshold)) +
        geom_ribbon(aes(x = expected, y = observed, ymin = lower_bound_log, ymax = upper_bound_log), alpha = 0.2) +
        geom_abline(slope = 1, intercept = 0, color = "red") +
        xlab(expression(atop("Expected", paste("(", -log[10], " p-value)")))) +
        ylab(expression(atop("Observed", paste("(", -log[10], " p-value)")))) +
        theme_minimal() +
        ylim(c(0, 6)) +
        theme(text = element_text(size = 24), plot.title = element_text(size = 24)) +
        ggtitle("3 Node Mediation No Pleiotropy or LD (Varying r2 threshold)")

ggsave(
    filename = print(
        file.path(output_mediation_LD_dir, sprintf('%s_QQ_mediation_3node_pval_r2_thresholds.jpeg', Sys.Date()))),
    plot = grouped_mediation_plot,
    units = "in", width = 12, height = 10
    )
