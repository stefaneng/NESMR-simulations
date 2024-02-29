#renv::load('/nfs/turbo/sph-jvmorr/NESMR/simulations')
# devtools::load_all('/nfs/turbo/sph-jvmorr/NESMR/esmr')

library(ggplot2)
library(dscrutils)
library(tidyr)
library(plyr)
library(dplyr)

# reticulate::use_condaenv('dsc')

#sim_path <- '/nfs/turbo/sph-jvmorr/NESMR/simulations/five_node_one_edge_complete'
output_dir <- file.path('results')

threshold <- 0.05

G <- matrix(
  c(0, 0, 0, 0, 0,
    sqrt(0.3), 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, -1 * sqrt(0.1), 0, 0, 0,
    0, -1 * sqrt(0.1), sqrt(0.2), sqrt(0.25), 0),
  nrow = 5,
  byrow = 5
)

h2 <- c(0.5, 0.3, 0.25, 0.4, 0.3)
N <- 40000
J <- 5e5
pi <- 500/J

B_true <- 1 * (abs(G) > 0)

G_edges <- which(lower.tri(G), arr.ind = TRUE)
G_edges_df <- cbind(
    as.data.frame(G_edges),
    beta = G[G_edges]) %>%
    arrange(desc(abs(beta))) %>%
    mutate(
        edge = paste(row, col, sep = '->')
      )

edge_order <- G_edges_df$edge

# Load data from dsc
# dscout <- dscquery(
#     dsc.outdir = sim_path,
#     targets    = c("simulate.results_df", "simulate.model_ll"))
dscout <- readRDS('analysis/five_node_complete_sim_v2.rds')

sim_results_df <- bind_rows(dscout$simulate.results_df, .id = "sim_idx")
long_sim_results_df <- sim_results_df %>%
    pivot_longer(
        cols = ends_with(c('_pvals', '_beta')),
        names_to = c("model", "metric"),
        names_pattern = "(.*)_(log10_pvals|beta)",
        values_to = "value") %>%
    rename(
        beta_true = beta
    ) %>%
    pivot_wider(
        names_from = "metric",
        values_from = "value") %>%
    rename(
        beta_est = beta
    ) %>%
    mutate(
        bias = beta_est - beta_true,
        edge = factor(
            paste(row, col, sep = '->'),
            levels = edge_order)
      )

long_sim_results_df %>%
  filter(model == "5-1" & edge == "5->1")

# write.csv(long_sim_results_df, file.path(output_dir, 'long_sim_results_df.csv'), row.names = FALSE)

bias_boxplot <- ggplot(long_sim_results_df) +
  geom_boxplot(aes(x = edge, y = bias)) +
  geom_hline(yintercept = 0, lty = 2, col = 'red') +
  facet_wrap(~model, ncol = 1) +
  theme_minimal() +
  theme(text = element_text(size = 24), plot.title = element_text(size = 24)) +
  labs(
    title = "Five node simulation complete DAG and one edge bias",
    caption = # add h2, N, J, pi to caption formated as expression. Use an expression or bquote to turn h^2 and pi into mathematically notation
      bquote(
        atop(
            "Sim parameters:" ~ h^2 ~ "=" ~ .(paste0(h2, collapse = ',')) ~ ", N = " ~ .(N) ~ ", J = " ~ .(J) ~ "," ~ pi ~ " = " ~ .(pi),
            "No LD or pleiotropy"
        )
      )
  )

ggsave(
  filename = print(file.path(
    output_dir,
    sprintf("%s_bias_boxplot.jpg", Sys.Date()
    ))
  ),
  plot = bias_boxplot,
  units = "in",
  width = 12, height = 14
)

five_one_log10_pvals <- long_sim_results_df %>%
  filter(model == "5-1" & edge == "5->1") %>%
  pull(log10_pvals)

# hist(10^five_one_log10_pvals)

## P-value histogram plot
pvalue_results <- long_sim_results_df %>%
  filter(correct_edge == FALSE & !is.na(log10_pvals)) %>%
  mutate(pval = 10^log10_pvals,
         model = ifelse(model == 'complete_model', 'Complete DAG', 'One Extra Edge'))

pvalue_hist <- ggplot(pvalue_results) +
    geom_histogram(
      aes(x = pval, after_stat(density)), breaks = seq(0, 1, by = 0.1)) +
  facet_grid(
    rows = vars(model),
    cols = vars(edge),
    scales = 'free_y'
  ) +
  theme_minimal() +
  theme(text = element_text(size = 24), plot.title = element_text(size = 24), axis.text.x = element_text(angle = 90)) +
  xlab("P-values") +
  labs(
    title = "Five node simulation complete DAG and one edge p-value distribution",
    caption = # add h2, N, J, pi to caption formated as expression. Use an expression or bquote to turn h^2 and pi into mathematically notation
      bquote(
        atop(
            "Sim parameters:" ~ h^2 ~ "=" ~ .(paste0(h2, collapse = ',')) ~ ", N = " ~ .(N) ~ ", J = " ~ .(J) ~ "," ~ pi ~ " = " ~ .(pi),
            "No LD or pleiotropy"
        )
      )
  )

ggsave(
  filename = print(file.path(
    output_dir,
    sprintf("%s_pvalue_hist.jpg", Sys.Date()
    ))
  ),
  plot = pvalue_hist,
  units = "in",
  width = 14, height = 10
)

## QQ plot
qqplot_log10 <- pvalue_results %>%
  group_by(model, edge) %>%
  arrange(pval) %>%
  mutate(
    observed = -log10_pvals,
    expected = -log10(ppoints(length(log10_pvals))),
    lower_bound_log = -log10(qbeta(0.025, seq_along(expected), rev(seq_along(expected)))),
    upper_bound_log = -log10(qbeta(0.975, seq_along(expected), rev(seq_along(expected))))
  ) %>%
  ggplot(aes(x = expected, y = observed)) +
  geom_point() +
  geom_ribbon(aes(ymin = lower_bound_log, ymax = upper_bound_log), alpha = 0.2) +
  facet_grid(
    rows = vars(model),
    cols = vars(edge),
    scales = 'free_y'
  ) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  xlab(expression(atop("Expected", paste("(", -log[10], " p-value)")))) +
  ylab(bquote(atop("Observed", paste("(", -log[10], " p-value)")))) +
  theme_minimal() +
  theme(text = element_text(size = 24)) +
    labs(
    title = expression("Five node simulation complete DAG and one edge", -log[10],  "QQ plots"),
    caption = bquote(
        atop(
            "Sim parameters:" ~ h^2 ~ "=" ~ .(paste0(h2, collapse = ',')) ~ ", N = " ~ .(N) ~ ", J = " ~ .(J) ~ "," ~ pi ~ " = " ~ .(pi),
            "No LD or pleiotropy"
        )
      )
  )


ggsave(
  filename = print(file.path(
    output_dir,
    sprintf("%s_log10_qqplot.jpg", Sys.Date()
    ))
  ),
  plot = qqplot_log10,
  units = "in",
  width = 14, height = 10
)


library(kableExtra)
pvalue_results %>%
  group_by(model, edge) %>%
  summarize(
    mean_bias = format(mean(bias), digits = 2, scientific = TRUE)
  ) %>%
  pivot_wider(
    names_from = "model",
    values_from = "mean_bias"
  ) %>%
  # Format the digits as scientific notation
  kable() %>%
  kable_styling()


