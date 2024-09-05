project_DAG_bootstrap <- function(
    total_est, total_est_se, reps = 100,
    s = 1.1,
    ...) {
  d <- ncol(total_est)
  non_diag_i <- -seq(1, d^2, by = d + 1)
  replicate(reps, {
    bootstrap_tot_est <- matrix(0, nrow = d, ncol = d)
    bootstrap_tot_est[non_diag_i] <- total_est[non_diag_i] + rnorm(d * (d - 1), mean = 0, sd = total_est_se[non_diag_i])
    project_to_DAG(
      bootstrap_tot_est,
      threshold_to_DAG = TRUE,
      maxit = 2000,
      trace = 5,
      s = s, # TODO: Why do we need s > 1? Always fails when s = 1
      ...
    )
  }, simplify = FALSE)
  # TODO: Standard error for each configuration ?
}

logdet_perm_fas_sim <- function(G) {
  d <- ncol(G)
  stopifnot(nrow(G) == d)
  h2 <- 0.3
  J <- 5000
  N <- 20000
  pi_J <- 0.1
  alpha <- 5e-8
  d <- n
  lambda <- qnorm(1 - alpha / 2)
  dat <- GWASBrewer::sim_mv(
    G = G,
    N = N,
    J = J,
    h2 = h2,
    pi = pi_J,
    sporadic_pleiotropy = TRUE,
    est_s = TRUE
  )

  Ztrue <- with(dat, beta_marg/se_beta_hat)
  pval_true <- 2*pnorm(-abs(Ztrue))
  minp <- apply(pval_true, 1, min)
  ix <- which(minp < alpha)

  # Fit MVMR matrix
  issue_dat <- NULL
  error_msg <- NULL

  rtn_list <- list()

  # Fit MVMR matrix
  mvmr_matrix <- with(dat, nesmr_complete_mvmr(
      beta_hat = beta_hat,
      se_beta_hat = s_estimate,
      pval_select = pval_true
    ))

  nesmr_full <- with(dat, nesmr_complete(
      beta_hat = beta_hat,
      se_beta_hat = s_estimate,
      pval_select = pval_true
    ))

  # TODO: Parametric bootstrap for each

  rtn_list$mvmr_DAG <- project_to_DAG(
    mvmr_matrix$beta,
    threshold_to_DAG = TRUE,
    maxit = 2000
  )

  rtn_list$nesmr_full_DAG <- project_to_DAG(
    nesmr_full$beta,
    threshold_to_DAG = TRUE,
    maxit = 2000
  )

  rtn_list$mvmr_DAG_bootstrap <- project_DAG_bootstrap(
    total_est = mvmr_matrix$beta_hat,
    total_est_se = mvmr_matrix$se_beta_hat
  )

  rtn_list$nesmr_full_bootstrap <- project_DAG_bootstrap(
    total_est = nesmr_full$beta_hat,
    total_est_se = nesmr_full$se_beta_hat
  )

  # TODO: Parameteric bootstrap for each

  rtn_list$all_ordering <- with(dat, nesmr_all_permn(
    beta_hat = beta_hat,
    se_beta_hat = s_estimate,
    B_templates = NULL,
    posterior_probs = TRUE,
    all_DAGs = TRUE
  ))

  return(rtn_list)
}

plot_sim_results <- function(z, G) {
  log_lik_order <- order(z$all_ordering$log_lik, decreasing = TRUE)
  # Sort each list by the log-likelihood
  z$all_ordering <- lapply(z$all_ordering, function(x) {
    x[log_lik_order]
  })

  correct_subset <- unlist(lapply(z$all_ordering$B_template, function(x) {
    all(x[G != 0] == 1)
  }))

  permn_config <- unlist(lapply(z$all_ordering$B_template, function(x) paste0(as.numeric(x), collapse = "")))

  simulation_results <- data.frame(
    correct_subset = correct_subset,
    posterior_probs = z$all_ordering$posterior_probs,
    config = permn_config
  )

  mvmr_bootstraps <- lapply(z$mvmr_DAG_bootstrap, function(x) (x$W != 0) + 0) %>%
    lapply(function(x) paste0(as.numeric(x), collapse = "")) %>%
    unlist() %>%
    table() %>%
    as.data.frame() %>%
    setNames(c('config', 'Freq'))
  mvmr_bootstraps$mvmr_bootstrap_percent <- mvmr_bootstraps$Freq / sum(mvmr_bootstraps$Freq)
  mvmr_bootstraps$Freq <- NULL

  nesmr_full_bootstraps <- lapply(z$nesmr_full_bootstrap, function(x) (x$W != 0) + 0) %>%
    lapply(function(x) paste0(as.numeric(x), collapse = "")) %>%
    unlist() %>%
    table() %>%
    as.data.frame() %>%
    setNames(c('config', 'Freq'))
  nesmr_full_bootstraps$nesmr_bootstrap_percent <- nesmr_full_bootstraps$Freq / sum(nesmr_full_bootstraps$Freq)
  nesmr_full_bootstraps$Freq <- NULL

  simulation_results <- merge(simulation_results, mvmr_bootstraps, by = 'config', all.x = TRUE)
  simulation_results <- merge(simulation_results, nesmr_full_bootstraps, by = 'config', all.x = TRUE)

  simulation_results <- arrange(simulation_results, correct_subset, posterior_probs)
  simulation_results$config_int <- rev(seq_len(nrow(simulation_results)))

  # Wide to long format for simulation results using pivot_longer
  sim_results_long <- pivot_longer(simulation_results,
                                   cols = c('mvmr_bootstrap_percent', 'nesmr_bootstrap_percent', 'posterior_probs'),
                                   names_to = 'method', values_to = 'percent')

  sim_results_long$percent[is.na(sim_results_long$percent)] <- 0

  #simulation_results <- melt(simulation_results, id.vars = c('config', 'correct_subset'))
  p <- ggplot(sim_results_long) +
    geom_point(aes(x = config_int - 0.5, y = percent,
                   color = method), position = position_dodge(0.5)) +
    geom_vline(
      xintercept = simulation_results$config_int - 1,
      color = ifelse(which(rev(!simulation_results$correct_subset))[1] == simulation_results$config_int, "orange", "grey50")
    ) +
    # Remove x-axis labels
    scale_x_continuous(
      breaks = simulation_results$config_int - 0.5,
      labels = simulation_results$config_int) +
    xlab('DAG Configuration') +
    theme_classic()

  return(list(
    data = sim_results_long,
    plot = p
  ))
}
# TODO: Want to plot with d groups:
# One for each number of parameters in the permutations

G3 <- matrix(
  c(0, 0, 0,
    0.1, 0, 0,
    0, 0, 0),
  nrow = 3,
  byrow = TRUE
)

G3_2 <- matrix(
  c(0, 0, 0,
    0.5, 0, 0,
    0.5, 0, 0),
  nrow = 3,
  byrow = TRUE
)

G3_3 <- matrix(
  c(0, 0, 0,
    0.5, 0, 0,
    0, 0.5, 0),
  nrow = 3,
  byrow = TRUE
)

set.seed(13)
G3_sim_results <- logdet_perm_fas_sim(G3)
set.seed(132)
G3_2_sim_results <- logdet_perm_fas_sim(G3_2)
G3_3_sim_results <- logdet_perm_fas_sim(G3_3)


plot_G3_res <- plot_sim_results(G3_sim_results, G3)
plot_G3_2_res <- plot_sim_results(G3_2_sim_results, G3_2)
plot_G3_3_res <- plot_sim_results(G3_3_sim_results, G3_3)

ggsave(print("~/Projects/NESMR-simulations/results/2024-09-05_three_node_one_edge.png"),
       plot_G3_res$plot, width = 12, height = 6)
ggsave(print("~/Projects/NESMR-simulations/results/2024-09-05_three_node_collider.png"),
       plot_G3_2_res$plot, width = 12, height = 6)
ggsave(print("~/Projects/NESMR-simulations/results/2024-09-05_three_node_mediation.png"),
       plot_G3_3_res$plot, width = 12, height = 6)
