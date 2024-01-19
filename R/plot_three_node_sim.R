
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

dscout.mediation <- dscquery(dsc.outdir = "/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node_mediation_no_LD",
                   targets    = c("simulate.dm_pvalue", "simulate.lrt_pvalue",
                   "simulate.true_model", "simulate.incorrect_model"),
                   ignore.missing.files = TRUE)

# Collider case
sim_pvalues <- dscout$simulate.lrt_pvalue
dm_pvalues <- exp(dscout$simulate.pval_dm)
pvals_df <- data.frame(lrt_pvalues = sim_pvalues, dm_pvalues = dm_pvalues)

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

mediation_pvals_df <- data.frame(
    dm_pvalues = dscout.mediation$simulate.dm_pvalue,
    lrt_pvalues = dscout.mediation$simulate.lrt_pvalue
    )

# Just to compare with
png(paste0(output_collider_dir, '/LRT_qq_plot.png'))
qqman::qq(pvals_df$lrt_pvalues)
dev.off()

png(paste0(output_collider_dir, '/DM_qq_plot.png'))
qqman::qq(pvals_df$dm_pvalues)
dev.off()

png(paste0(output_mediation_LD_dir, '/LRT_qq_plot.png'))
qqman::qq(pvals_df_LD$lrt_pvalues)
dev.off()

png(paste0(output_mediation_LD_dir, '/DM_qq_plot.png'))
qqman::qq(pvals_df_LD$dm_pvalues)
dev.off()

png(paste0(output_dir, '/dm_qq_plot_no_LD.png'))
qqman::qq(mediation_pvals_df$dm_pvalues)
dev.off()

png(paste0(output_dir, '/LRT_qq_plot_no_LD.png'))
qqman::qq(mediation_pvals_df$lrt_pvalues)
dev.off()

# Call the function for LRT P-values
create_and_save_plots(pvals_df, "lrt_pvalues", "lrt", output_dir = output_collider_dir)

# Call the function for Delta Method P-values
# create_and_save_plots(pvals_df, "dm_pvalues", "delta_method", output_dir = output_dir)

# Mediation p-values
create_and_save_plots(mediation_pvals_df, "lrt_pvalues", "lrt", suffix = "medidation_no_pleiotropy", output_dir)

# Call the function for Delta Method P-values
create_and_save_plots(mediation_pvals_df, "dm_pvalues", "delta_method", suffix = "medidation_no_pleiotropy", output_dir)

# Mediation p-values
create_and_save_plots(mediation_pvals_df, "lrt_pvalues", "lrt", suffix = "medidation_no_pleiotropy", output_dir)

# Call the function for Delta Method P-values
create_and_save_plots(mediation_pvals_df, "dm_pvalues", "delta_method", suffix = "medidation_no_pleiotropy", output_dir)


cat('mean(LD LRT p-values < 0.05) =', mean(dscout.mediation.LD$simulate.lrt_pvalue < 0.05, na.rm =TRUE), '\n')
cat('mean(LD Delta method p-values < 0.05) =', mean(dscout.mediation.LD$simulate.dm_pvalue < 0.05, na.rm =TRUE), '\n')

cat('mean(No LD LRT p-values < 0.05) =', mean(dscout.mediation$simulate.lrt_pvalue < 0.05, na.rm =TRUE), '\n')
cat('mean(No LD Delta method p-values < 0.05) =', mean(dscout.mediation$simulate.dm_pvalue < 0.05, na.rm =TRUE), '\n')
