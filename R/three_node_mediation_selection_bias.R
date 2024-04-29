library(GWASBrewer)
library(dplyr)
library(ggplot2)
library(forcats)
library(reshape2)
library(parallel)

devtools::load_all('../esmr/')

# Note: This was not run on cluster but on local machine
#  This illustrates that we have winner's curse

G <- matrix(
  c(0, 0, 0,
    sqrt(0.4), 0, 0,
    0, sqrt(0.2), 0),
  nrow = 3,
  byrow = TRUE
)

B_true <- G
B_true[!B_true == 0] <- 1

# Incorrect configurations
B_inc_mediation <- matrix(
  c(0, 0, 0,
    1, 0, 0,
    1, 1, 0),
  nrow = 3,
  byrow = TRUE
)

res <- parallel::mclapply(1:100, function(x) {
  h2 <- c(0.5, 0.3, 0.25)
  ## simulate summary statistics
  dat <- GWASBrewer::sim_mv(
    G = G,
    N = 40000,
    J = 5e5,
    h2 = h2,
    pi = 500/5e5,
    sporadic_pleiotropy = FALSE,
    est_s = TRUE
  )

  Ztrue <- with(dat, beta_marg/se_beta_hat)
  pval_true <- 2*pnorm(-abs(Ztrue))
  minp <- apply(pval_true, 1, min)

  ix1 <- which(minp < 5e-8)

  true_model <- with(dat,
                          esmr(
                            beta_hat_X = beta_hat[ix1,],
                            se_X = s_estimate[ix1,],
                            pval_thresh = 1,
                            G = diag(3),
                            direct_effect_template = B_true,
                            max_iter = 300))

  incorrect_model <- with(dat,
                               esmr(
                                 beta_hat_X = beta_hat[ix1,],
                                 se_X = s_estimate[ix1,],
                                 pval_thresh = 1,
                                 G = diag(3),
                                 direct_effect_template = B_inc_mediation,
                                 max_iter = 300))

  true_ll <- logLik.esmr(true_model)
  incorrect_ll <- logLik.esmr(incorrect_model)

  lrt_pvalue <- pchisq(
    - 2 * (true_ll - incorrect_ll),
    df = 1,
    lower.tail = FALSE
  )

  dm_pvalue <- exp(incorrect_model$pvals_dm[3, 1])
  setNames(list(true_ll, incorrect_ll, lrt_pvalue, dm_pvalue,
                beta_1 = incorrect_model$direct_effects[2,1],
                beta_2 = incorrect_model$direct_effects[3,2],
                beta_3 = incorrect_model$direct_effects[3,1]),
           c('true_ll', 'incorrect_ll', 'lrt_pvalue', 'dm_pvalue', 'beta_1', 'beta_2', 'beta_3'))
})

res_merged <- melt(res)

# saveRDS(res_merged, 'results/three_node_mediation/three_node_mediation_selection_bias.rds')

true_beta <- data.frame(
  beta = c('beta_1', 'beta_2', 'beta_3'),
  true_value = c(sqrt(0.4), sqrt(0.2), 0)
)

beta_df <- left_join(
  res_merged %>% filter(L2 %in% paste0('beta_', 1:3)) %>% rename(beta = L2),
  true_beta,
  by = 'beta'
) %>%
  mutate(
    bias = value - true_value,
    beta = forcats::fct_recode(
      factor(beta),
      `beta["2,1"] == sqrt(0.4)` = "beta_1",
      `beta["3,2"] == sqrt(0.2)` = "beta_2",
      `beta["3,1"] == 0` = "beta_3"
    )
  )

x_labels <- lapply(levels(beta_df$beta), function(x) parse(text = x))

beta_boxplot <- beta_df %>%
  ggplot(aes(x = beta, y = bias)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty = 2, col = 'red') +
  scale_x_discrete(labels = unlist(x_labels, recursive = F)) +
  xlab("") +
  ylab(expression(hat(beta) - beta)) +
  theme_minimal() +
  theme(text = element_text(size = 24),
        plot.title = element_text(size = 24)) +
  ggtitle("3 Node without in-selection bias")

ggsave(
  filename = print(
    sprintf('%s_beta_in_selection_bias.jpeg', Sys.Date())),
  plot = beta_boxplot,
  units = "in", width = 8, height = 8
)

alpha <- 0.05

pval_data <- res_merged %>%
  filter(L2 %in% c('lrt_pvalue', 'dm_pvalue')) %>%
  mutate(observed = -log10(value)) %>%
  arrange(L2, desc(observed)) %>%
  rename(
    'p_value_method' = L2,
    'iteration' = L1
  )

expected <- -log10(ppoints(nrow(pval_data) / 2))
lower_bound <- qbeta(alpha / 2, seq_along(expected), rev(seq_along(expected)))
upper_bound <- qbeta(1 - alpha / 2, seq_along(expected), rev(seq_along(expected)))
lower_bound_log <- -log10(lower_bound)
upper_bound_log <- -log10(upper_bound)

pval_data$expected <- rep(expected, 2)
pval_data$lower_bound_log <- rep(lower_bound_log, 2)
pval_data$upper_bound_log <- rep(upper_bound_log, 2)

# saveRDS(pval_data, 'results/three_node_mediation/three_node_mediation_qqplot.rds')

pval_qq_plot <- ggplot(pval_data, aes(x = expected, y = observed)) +
  geom_point(aes(group = p_value_method, color = p_value_method)) +
  geom_ribbon(aes(ymin = lower_bound_log, ymax = upper_bound_log), alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  xlab(expression(atop("Expected", paste("(", -log[10], " p-value)")))) +
  ylab(expression(atop("Observed", paste("(", -log[10], " p-value)")))) +
  theme_minimal() +
  scale_color_discrete(breaks = c('dm_pvalue', 'lrt_pvalue'), labels=c("DM","LRT"),
                       name="P-value type") +
  ylim(c(0, 6)) +
  theme(text = element_text(size = 24), plot.title = element_text(size = 20)) +
  ggtitle("3 Node without in-selection bias")

ggsave(
  filename = print(
    sprintf('%s_log10_qqplot_in_selection_bias.jpeg', Sys.Date())),
  plot = pval_qq_plot,
  units = "in", width = 8, height = 8
)
