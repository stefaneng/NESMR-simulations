library(ggplot2)
library(ggpubr)
library(ggh4x)
# TODO: Point this at cluster results
wib_data <- readRDS('data/2024-06-11_ma_wib_3node_no_errors.rds')
dir <- 'results/ma_wib_3node'

# Set ggplot theme to minimal
theme_set(theme_minimal())

today <- format(Sys.Date(), "%Y-%m-%d")

wib_data <- wib_data %>%
  rename_with(
    ~ gsub("(simulate|fit_esmr)[.]", "", .x),
    starts_with("simulate") | starts_with("fit_esmr")
  ) %>%
  mutate(
    # Edge is null if
    null_edge = (graph_type == "mediation" & edge == "V3_V1") |
      (graph_type == "correlated" & edge == "V2_V1"),
    beta_name = case_when(
      graph_type == "mediation" & edge == "V2_V1" ~ "beta_1",
      graph_type == "mediation" & edge == "V3_V1" ~ "beta_null",
      graph_type == "correlated" & edge == "V3_V1" ~ "beta_1",
      graph_type == "correlated" & edge == "V2_V1" ~ "beta_null",
      edge == "V3_V2" ~ "beta_2"
    )
  ) %>%
  mutate(
    graph_type = gsub("correlated", "confounded", graph_type),
    ma = grepl("ma", model)
  )

grouped_wib_data <- wib_data %>%
  group_by(
    graph_type,
    alpha,
    model,
    N,
    eta,
    edge,
    beta_name
  ) %>%
  summarise(
    mean_bias = mean(bias),
    se_bias = sd(bias, na.rm = TRUE),
    mse = mean(bias^2),
    se_mse = sd(bias^2, na.rm = TRUE),
    reps = n()
  )

# There was an issue in the simulation where "true" variants
# at the specified threshold were not included in the results when alpha = 1
# cursed and true are the same when alpha = 1 so can just replace the model name
replace_true_model <- grouped_wib_data %>%
  filter(alpha == 1 & model == "cursed_model")
replace_true_model$model <- "true_model"

grouped_wib_data <- grouped_wib_data %>%
  filter(!(alpha == 1 & model == "true_model")) %>%
  bind_rows(replace_true_model)

bias_plots <- lapply(unique(grouped_wib_data$N), function(n) {
  grouped_wib_data %>%
    filter(N == n) %>%
    ggplot(., aes(x = -log10(alpha), color = model)) +
    geom_line(aes(y = mean_bias)) +
    geom_point(aes(y = mean_bias)) +
    geom_errorbar(
      aes(ymin = mean_bias - 1.96 * se_bias / sqrt(reps),
          ymax = mean_bias + 1.96 * se_bias / sqrt(reps)),
      width = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_continuous(
      breaks = -log10(c(1, 1e-3, 1e-5, 1e-8, 1e-10)),
      labels = c(expression(10^0), expression(10^-3), expression(10^-5), expression(10^-8), expression(10^-10))
    ) +
    ggh4x::facet_grid2(cols = vars(beta_name), rows = vars(graph_type), scales = "free",
                       labeller = labeller(
                         beta_name = as_labeller(c(beta_1="beta[1]", beta_2="beta[2]", beta_null='beta["null"]'), default = label_parsed)
                       ),
                independent = "y") +
    scale_color_discrete(
      # rename ma_1 to Ma eta=0.5 and ma_2 to Ma eta=1
      labels = c(
        "In-Sample\nCursed",
        "Ma\neta = 0.5",
        "Ma\neta = 1",
        "Out-Sample")
    ) +
    labs(
      x = "Variant inclusion level",
      y = "Bias"
    ) +
    ggtitle(
      sprintf("N = %s", n)
    ) +
    theme(
      strip.text.y = element_text(angle = 0)
    )
})

bias_plots_merged <- ggarrange(
  bias_plots[[1]],
  bias_plots[[2]],
  ncol = 1, nrow = 2,
  common.legend = TRUE,
  legend = "bottom"
)

ggsave(
  filename = print(file.path(dir, paste0(today, "_ma_wib_3node_bias_plot.jpeg"))),
  plot = bias_plots_merged,
  units = "in",
  width = 10, height = 12
)

## TODO: Refactor so that we have one function that plots bias & MSE
mse_plots <- lapply(unique(grouped_wib_data$N), function(n) {
  grouped_wib_data %>%
    filter(N == n) %>%
    ggplot(., aes(x = -log10(alpha), color = model)) +
    geom_line(aes(y = mse)) +
    geom_point(aes(y = mse)) +
    geom_errorbar(
      aes(ymin = mse - 1.96 * se_mse / sqrt(reps),
          ymax = mse + 1.96 * se_mse / sqrt(reps)),
      width = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_continuous(
      breaks = -log10(c(1, 1e-3, 1e-5, 1e-8, 1e-10)),
      labels = c(expression(10^0), expression(10^-3), expression(10^-5), expression(10^-8), expression(10^-10))
    ) +
    ggh4x::facet_grid2(cols = vars(beta_name), rows = vars(graph_type), scales = "free",
                       labeller = labeller(
                         beta_name = as_labeller(c(beta_1="beta[1]", beta_2="beta[2]", beta_null='beta["null"]'), default = label_parsed)
                       ),
                       independent = "y") +
    scale_color_discrete(
      # rename ma_1 to Ma eta=0.5 and ma_2 to Ma eta=1
      labels = c(
        "In-Sample\nCursed",
        "Ma\neta = 0.5",
        "Ma\neta = 1",
        "Out-Sample")
    ) +
    labs(
      x = "Variant inclusion level",
      y = "MSE"
    ) +
    ggtitle(
      sprintf("N = %s", n)
    ) +
    theme(
      strip.text.y = element_text(angle = 0)
    )
})

mse_plots_merged <- ggarrange(
  mse_plots[[1]],
  mse_plots[[2]],
  ncol = 1, nrow = 2,
  common.legend = TRUE,
  legend = "bottom"
)

ggsave(
  filename = print(file.path(dir, paste0(today, "_ma_wib_3node_mse_plot.jpeg"))),
  plot = mse_plots_merged,
  units = "in",
  width = 10, height = 12
)

wib_data %>%
  arrange(
    graph_type,
    model,
    N,
    alpha,
    -log10_pvals
  )

## Calibrated p-values?
pvalue_histograms <- wib_data %>%
  filter(null_edge) %>%
  ggplot(., aes(x = 10^(-log10_pvals))) +
  # x axis labels 0, 0.25, 0.5, 0.75, 1
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), ) +
  labs(
    x = "P-value",
    y = "Count",
  ) +
  facet_grid(
    rows = vars(graph_type, model),
    cols = vars(N, alpha),
    labeller = label_both) +
  geom_histogram(bins = 20, fill = "lightblue") + theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    strip.text.y = element_text(angle = 0))

ggsave(
  filename = print(file.path(dir, paste0(today, "_ma_wib_3node_pvalue_hists.pdf"))),
  plot = pvalue_histograms,
  units = "in",
  width = 14, height = 12
)
