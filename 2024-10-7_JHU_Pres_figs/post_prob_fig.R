library(esmr)
n <- 10
three_dags <- tail(generate_dags_with_permutations(3), n = n)

example_data <- data.frame(
  permn_config = unlist(lapply(three_dags, function(x) paste0(as.numeric(x), collapse = ""))),
  config_axis = unlist(lapply(three_dags, matrix_to_str)),
  percent = c(0.4, 0.2, 0.15, 0.2, 0.07, 0.07, 0.05, 0.025, 0.025, 0.01),
  correct_config = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
)

example_data <- example_data[order(example_data$percent, decreasing = TRUE),]
example_data$config_int <- seq_len(n)

cs <- which(example_data$correct_config)[1]

fig_plot <- ggplot(example_data) +
  # Fill in background area if we have correct configuration
  annotate(
    'rect',
    xmin = example_data$config_int[cs] ,
    xmax = example_data$config_int[cs] + 1,
    ymin = -0.05,
    ymax = 1.05,
    alpha = 0.5,
    fill = "orange"
  ) +
  geom_point(aes(x = config_int + 0.5, y = percent),
             size = 3) +
  scale_x_continuous(
    breaks = example_data$config_int + 0.5,
    labels = example_data$config_axis
  ) +
  ylab("Posterior Probability") +
  xlab('DAG Configuration') +
  coord_cartesian(
    xlim = c(1.2, max(example_data$config_int) + 1),
    ylim = c(0,1), clip = "off"
  ) +
  theme_classic()


today_date <- format(Sys.Date(), "%Y-%m-%d")
ggsave(print(sprintf("~/Projects/NESMR-simulations/results/%s_%s.png", today_date, "post_prob_example_fig")),
       fig_plot, width = 5, height = 2)
