
renv::load('/nfs/turbo/sph-jvmorr/NESMR/simulations')
devtools::load_all('/nfs/turbo/sph-jvmorr/NESMR/esmr')
library(ggplot2)
library(dscrutils)

source("/nfs/turbo/sph-jvmorr/NESMR/simulations/R/pvalue_plots.R")

reticulate::use_condaenv('dsc')

# Define the output directory
output_collider_dir <- "/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node_output/figures"
output_dir <- "/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node_mediation_no_LD/figures"
output_mediation_LD_dir <- "/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node_mediation/figures"

dscout <- dscquery(dsc.outdir = "/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node_output",
                   targets    = c("simulate.lrt_pvalue", "simulate.pval_dm"),
                   ignore.missing.files = TRUE)

dscout.mediation.LD <- dscquery(dsc.outdir = "/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node_mediation",
                   targets    = c("simulate.dm_pvalue", "simulate.lrt_pvalue",
                   "simulate.true_model", "simulate.incorrect_model"),
                   ignore.missing.files = TRUE)

dscout.mediation <- dscquery(dsc.outdir = "/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node_mediation_no_LD_test",
                   targets    = c("simulate.dm_pvalue", "simulate.lrt_pvalue",
                   "simulate.true_model", "simulate.incorrect_model"))

# TODO: Now we have p-value thresholds to deal with

# Collider case
sim_pvalues <- dscout$simulate.lrt_pvalue
dm_pvalues <- exp(dscout$simulate.pval_dm)
pvals_df <- data.frame(lrt_pvalues = sim_pvalues, dm_pvalues = dm_pvalues)
pvals_long_df <- rbind(
    data.frame(pvalue = sim_pvalues, pvalue_type = "lrt"),
    data.frame(pvalue = dm_pvalues, pvalue_type = "dm")
)

# Mediation case with LD
pvals_df_LD <- data.frame(
    lrt_pvalues = dscout.mediation.LD$simulate.lrt_pvalue,
    dm_pvalues = dscout.mediation.LD$simulate.dm_pvalue
)

if (with(dscout.mediation.LD, cor(simulate.dm_pvalue, simulate.lrt_pvalue)) == 1) {
    # Note: This is just a temporary fix as I wrote the wrong target
    lrt_pvalue <- unlist(lapply(seq_along(dscout.mediation.LD$simulate.true_model), function(i) {
        true_model <- dscout.mediation.LD$simulate.true_model[[i]]
        incorrect_model <- dscout.mediation.LD$simulate.incorrect_model[[i]]

        true_ll <- logLik.esmr(true_model)
        incorrect_ll <- logLik.esmr(incorrect_model)

        lrt_pvalue <- pchisq(
            - 2 * (true_ll - incorrect_ll),
            df = 1,
            lower.tail = FALSE
        )
    }))

    pvals_df_LD$lrt_pvalues <- lrt_pvalue
}

pval_thresholds <- c(5e-8, 5e-9, 5e-10)

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

n_pvals <- nrow(mediation_pvals_df) / (2 * 2 * 3)
alpha <- 0.05
lower_bound <- qbeta(alpha / 2, seq_len(n_pvals), rev(seq_len(n_pvals)))
upper_bound <- qbeta(1 - alpha / 2, seq_len(n_pvals), rev(seq_len(n_pvals)))

grouped_mediation_plot <- mediation_pvals_df %>%
    group_by(pvalue_type, pval_threshold, beta_var) %>%
    arrange(pvalue) %>%
    mutate(
        pval_threshold = as.factor(pval_threshold),
        observed = -log10(pvalue),
        expected = -log10(ppoints(n_pvals)),
        lower_bound_log = -log10(lower_bound),
        upper_bound_log = -log10(upper_bound)) %>%
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
        filename = 'test_QQ.jpeg',
        plot = grouped_mediation_plot,
        units = "in", width = 12, height = 8
    )

grouped_collider_plot <- pvals_long_df %>%
    group_by(pvalue_type) %>%
    arrange(pvalue) %>%
    mutate(
        observed = -log10(pvalue),
        expected = -log10(ppoints(nrow(pvals_df))),
        lower_bound_log = -log10(lower_bound),
        upper_bound_log = -log10(upper_bound)) %>%
    ggplot() +
    facet_grid(cols = vars(pvalue_type)) +
        geom_point(aes(x = expected, y = observed)) +
        geom_ribbon(aes(x = expected, y = observed, ymin = lower_bound_log, ymax = upper_bound_log), alpha = 0.2) +
        geom_abline(slope = 1, intercept = 0, color = "red") +
        xlab(expression(atop("Expected", paste("(", -log[10], " p-value)")))) +
        ylab(expression(atop("Observed", paste("(", -log[10], " p-value)")))) +
        theme_minimal() +
        ylim(c(0, 6)) +
        theme(text = element_text(size = 24), plot.title = element_text(size = 24)) +
        ggtitle("3 Node Collider with LD")

ggsave(
        filename = 'collider_test_QQ.jpeg',
        plot = grouped_collider_plot,
        units = "in", width = 12, height = 8
    )

# mediation_pvals_df <- data.frame(
#     dm_pvalues = dscout.mediation$simulate.dm_pvalue,
#     lrt_pvalues = dscout.mediation$simulate.lrt_pvalue
#     )

# Just to compare with
# png(paste0(output_collider_dir, '/LRT_qq_plot.png'))
# qqman::qq(pvals_df$lrt_pvalues)
# dev.off()

# png(paste0(output_collider_dir, '/DM_qq_plot.png'))
# qqman::qq(pvals_df$dm_pvalues)
# dev.off()

# png(paste0(output_mediation_LD_dir, '/LRT_qq_plot.png'))
# qqman::qq(pvals_df_LD$lrt_pvalues)
# dev.off()

# png(paste0(output_mediation_LD_dir, '/DM_qq_plot.png'))
# qqman::qq(pvals_df_LD$dm_pvalues)
# dev.off()

# png(paste0(output_dir, '/dm_qq_plot_no_LD.png'))
# qqman::qq(mediation_pvals_df$dm_pvalues)
# dev.off()

# png(paste0(output_dir, '/LRT_qq_plot_no_LD.png'))
# qqman::qq(mediation_pvals_df$lrt_pvalues)
# dev.off()

# Call the function for LRT P-values
create_and_save_plots(
    pvals_df, "lrt_pvalues", "lrt", output_dir = output_collider_dir,
    main = "3 Node Collider")

# Call the function for Delta Method P-values
create_and_save_plots(
    pvals_df, "dm_pvalues", "delta_method", output_dir = output_collider_dir,
    main = "3 Node Collider")

create_and_save_plots(
        mediation_pvals_df,
        pvalue_col = paste0("lrt_pvalue_", pval_thres),
        "pvalue",
        suffix = paste0("medidation_no_pleiotropy_", pval_thres),
        pvalue_type = "pvalue_type",
        group = "pval_threshold",
        output_dir = output_dir,
        main = paste0("3 Node Mediation\nNo Pleiotropy or LD: p < ", pval_thres))

# Mediation p-values
for (pval_thres in pval_thresholds) {
    create_and_save_plots(
        mediation_pvals_df, pvalue_col = paste0("lrt_pvalue_", pval_thres), "lrt", suffix = paste0("medidation_no_pleiotropy_", pval_thres), output_dir,
        main = paste0("3 Node Mediation\nNo Pleiotropy or LD: p < ", pval_thres))

    # Call the function for Delta Method P-values
    create_and_save_plots(
        mediation_pvals_df, pvalue_col = paste0("dm_pvalue_", pval_thres), "delta_method", suffix = paste0("medidation_no_pleiotropy_", pval_thres), output_dir,
        main = paste0("3 Node Mediation\nNo Pleiotropy or LD: p < ", pval_thres))
}


# Mediation p-values with LD
create_and_save_plots(
    pvals_df_LD, "lrt_pvalues", "lrt", suffix = "medidation_with_LD", output_mediation_LD_dir,
    main = "3 Node Mediation with LD")

# Call the function for Delta Method P-values
create_and_save_plots(pvals_df_LD, "dm_pvalues", "delta_method", suffix = "medidation_with_LD", output_mediation_LD_dir,
    main = "3 Node Mediation with LD")

cat('mean(LD LRT p-values < 0.05) =', mean(dscout.mediation.LD$simulate.lrt_pvalue < 0.05, na.rm =TRUE), '\n')
cat('mean(LD Delta method p-values < 0.05) =', mean(dscout.mediation.LD$simulate.dm_pvalue < 0.05, na.rm =TRUE), '\n')

cat('mean(No LD LRT p-values < 0.05) =', mean(dscout.mediation$simulate.lrt_pvalue < 0.05, na.rm =TRUE), '\n')
cat('mean(No LD Delta method p-values < 0.05) =', mean(dscout.mediation$simulate.dm_pvalue < 0.05, na.rm =TRUE), '\n')

