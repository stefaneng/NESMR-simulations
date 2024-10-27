library(dplyr)
library(plyr)
library(tidyr)
library(ggplot2)
devtools::load_all("/Users/stefaneng/Projects/esmr")

logdet_perm_fas_sim <- function(G, permute = TRUE, faster_init = FALSE) {
  d <- ncol(G)
  stopifnot(nrow(G) == d)
  h2 <- 0.3
  J <- 5000
  N <- 20000
  pi_J <- 0.1
  alpha <- 5e-8
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

  nesmr_full <- with(dat, nesmr_complete(
      beta_hat = beta_hat,
      se_beta_hat = s_estimate,
      pval_select = pval_true
    ))

  rtn_list$nesmr_full_DAG <- project_to_DAG(
    nesmr_full$beta_hat,
    threshold_to_DAG = TRUE,
    penalty = "L2",
    maxit = 2000
  )

  rtn_list$nesmr_full_bootstrap <- project_to_DAG_bootstrap(
    total_est = nesmr_full$beta_hat,
    total_est_se = nesmr_full$se_beta_hat
  )

  unique_configs <- unique(lapply(rtn_list$nesmr_full_bootstrap, function(x) {
    (x$W != 0) + 0
  }))

  # TODO: For each of the config, do the backselect test
  # Want to find someway to avoid repeating this over again..

  # Generate all permutations
  # Simple way is to just find all the unique configs that are a subset of the unique configs

  # Fit all of the configurations found in the bootstrap
  refit_mods <- lapply(unique_configs, function(B) {
    with(dat, esmr(
      beta_hat_X = beta_hat,
      se_X = s_estimate,
      variant_ix = ix,
      G = diag(d), # required for network problem
      direct_effect_template = B))
  })

  rtn_list$backselect_mods <- nesmr_backselect(
    refit_mods,
    beta_hat = dat$beta_hat,
    se_beta_hat = dat$s_estimate,
    Z_true = Ztrue
    )

  # TODO: Refit all subsets of DAGs that are significant

  all_perms <- generate_dags_with_permutations(d)

  rtn_list$all_ordering <- with(dat, nesmr_all_permn(
    beta_hat = beta_hat,
    se_beta_hat = s_estimate,
    posterior_probs = TRUE,
    all_DAGs = TRUE,
    return_model = TRUE,
    B_templates = all_perms,
    Z_true = Ztrue,
    direct_effect_init = nesmr_full$beta_hat
  ))

  log_lik_order <- order(rtn_list$all_ordering$log_lik, decreasing = TRUE)
  # Sort each list by the log-likelihood
  rtn_list$all_ordering <- lapply(rtn_list$all_ordering, function(x) {
    x[log_lik_order]
  })

  return(rtn_list)
}

process_sim_results <- function(z, G) {
  correct_subset <- unlist(lapply(z$all_ordering$B_template, function(x) {
    all(x[G != 0] == 1)
  }))

  correct_config <- unlist(lapply(z$all_ordering$B_template, function(x) {
    all(x == (G != 0))
  }))

  permn_config <- unlist(lapply(z$all_ordering$B_template, function(x) paste0(as.numeric(x), collapse = "")))
  permn_config_axis <- unlist(lapply(z$all_ordering$B_template, matrix_to_str))

  simulation_results <- data.frame(
    correct_subset = correct_subset,
    correct_config = correct_config,
    posterior_probs = z$all_ordering$posterior_probs,
    posterior_probs_n_edge = z$all_ordering$posterior_probs_n_edge,
    log_lik = z$all_ordering$log_lik,
    AIC = z$all_ordering$aic,
    config = permn_config,
    config_axis = permn_config_axis
  )

  .bootstrap_to_df <- function(x, method = "") {
    tmp_res <- lapply(x, function(x) (x != 0) + 0) %>%
      lapply(function(x) paste0(as.numeric(x), collapse = "")) %>%
      unlist() %>%
      table() %>%
      as.data.frame() %>%
      setNames(c('config', 'Freq'))
    tmp_res[[method]] <- tmp_res$Freq / sum(tmp_res$Freq)
    tmp_res$Freq <- NULL
    tmp_res$config <- as.character(tmp_res$config)
    tmp_res
  }

  nesmr_full_bootstraps <- .bootstrap_to_df(lapply(z$nesmr_full_bootstrap, function(x) x$W), method = "nesmr_project_percent")
  nesmr_full_bootstraps$algo <- "logdet"

  simulation_results <- merge(simulation_results, nesmr_full_bootstraps, by = 'config', all.x = TRUE)

  simulation_results <- arrange(simulation_results, correct_subset, posterior_probs)
  simulation_results$config_int <- rev(seq_len(nrow(simulation_results)))

  simulation_results
}

plot_sim_results <- function(sim_results, G, max_configs = 30) {
  # Remove small results
  sim_results <- sim_results[
    sim_results$config_int <= max_configs,
  ]

  # Wide to long format for simulation results using pivot_longer
  sim_results_long <- pivot_longer(
    sim_results,
    cols = c(
      'nesmr_project_percent',
      'posterior_probs_n_edge',
      'posterior_probs'
    ),
    names_to = 'method', values_to = 'percent') %>%
    arrange(desc(correct_subset), desc(log_lik))

  sim_results_long$algo <- ifelse(grepl('project', sim_results_long$method), 'logdet', 'loglik')

  sim_results_long$percent[is.na(sim_results_long$percent)] <- 0

  cs <- which(sim_results$correct_config)[1]

  min_AIC <- min(sim_results$AIC)

  #sim_results <- melt(sim_results, id.vars = c('config', 'correct_subset'))
  p <- ggplot(sim_results_long) +
    # Fill in background area if we have correct configuration
    annotate(
      'rect',
      xmin = sim_results$config_int[cs] - 1,
      xmax = sim_results$config_int[cs],
      ymin = -0.05,
      ymax = 1.05,
      alpha = 0.5,
      fill = "orange"
    ) +
    geom_point(aes(x = config_int - 0.5, y = percent,
                   color = method,
                   shape = method), position = position_dodge(0.5),
               size = 1.5) +
    geom_vline(
      xintercept = sim_results$config_int - 1,
      color = ifelse(which(rev(!sim_results$correct_subset))[1] == sim_results$config_int, "orange", "grey50")
    ) +
    annotate(
      "text",
      x = sim_results$config_int - 0.5,
      y = 1,
      label = round(min_AIC - sim_results$AIC, 1)
    ) +
    # Remove x-axis labels
    scale_x_continuous(
      breaks = sim_results$config_int - 0.5,
      labels = sim_results$config_axis
      ) +
    annotate(
      "label",
      x = 12.5,
      y = 0.5,
      label = paste0("True graph:\n", paste0(capture.output(print(G)), collapse = "\n"))
    ) +
    xlab('DAG Configuration') +
    coord_cartesian(
      xlim = c(1.2, max(sim_results$config_int) + 1),
      ylim = c(0,1), clip = "off"
    ) +
    theme_classic()

  return(list(
    sim_results_long = sim_results_long,
    plot = p
  ))
}

G3 <- matrix(
  c(0, 0, 0,
    0.05, 0, 0,
    0, 0, 0),
  nrow = 3,
  byrow = TRUE
)

G3_2 <- matrix(
  c(0, 0, 0,
    0.075, 0, 0,
    0.075, 0, 0),
  nrow = 3,
  byrow = TRUE
)

G3_3 <- matrix(
  c(0, 0, 0,
    0.075, 0, 0,
    0, 0.075, 0),
  nrow = 3,
  byrow = TRUE
)


G4 <- matrix(
  c(0, 0, 0, 0,
    0.075, 0, 0, 0,
    0, 0.075, 0, 0,
    0, 0, 0, 0),
  nrow = 4,
  byrow = TRUE
)

sim_params <- list(
  h2 = 0.3,
  J = 5000,
  N = 20000,
  pi_J = 0.1,
  alpha = 5e-8
)

sim_params_str <- paste0(names(sim_params), " = ", sim_params, collapse = ", ")

multi_sim <- function(n, G, faster_init = FALSE) {
  sim_results <- list()
  plot_res <- list()
  seeds <- sample(1:1e6, n)
  for (i in seq_along(seeds)) {
    s <- seeds[[i]]
    set.seed(s)
    sim_results[[i]] <- logdet_perm_fas_sim(G, faster_init = faster_init)
    sim_results[[i]]$sim_results <- process_sim_results(sim_results[[i]], G)
    plot_res[[i]] <- plot_sim_results(sim_results[[i]]$sim_results, G)
  }

  return(list(
    sim_results = sim_results,
    plot_results = plot_res,
    seeds = seeds
  ))
}

plot_n_sim <- function(x, G) {
  all_res <- bind_rows(lapply(x$sim_results, `[[`, "sim_results"), .id = "sim") %>%
    group_by(sim) %>%
    mutate(AIC_diff = min(AIC) - AIC)

  combined_post_prob <- all_res %>%
    ggplot(aes(x = config_int, y = posterior_probs, color = correct_subset, shape = correct_config)) +
    geom_point() +
    facet_wrap(~sim, ncol = 2) +
    # Remove x-axis labels
    scale_x_continuous(
      breaks = unique(all_res$config_int),
      labels = unique(all_res$config_axis)
    ) +
    xlab("DAG Configuration") +
    ylab("Posterior uniform prior probs") +
    labs(
      title = paste0("True graph:\n", paste0(capture.output(print(G)), collapse = "\n")),
      subtitle = sim_params_str
    ) +
    theme(
      title = element_text(size = 12)
    ) +
    theme_minimal()

  combined_AIC <- all_res %>%
    ggplot(aes(x = config_int, y = AIC_diff, color = correct_subset, shape = correct_config)) +
    geom_point() +
    facet_wrap(~sim, ncol = 2) +
    # Remove x-axis labels
    scale_x_continuous(
      breaks = unique(all_res$config_int),
      labels = unique(all_res$config_axis)
    ) +
    coord_cartesian(
      clip = "off"
    ) +
    xlab("DAG Configuration") +
    ylab("MIN(AIC) - AIC") +
    labs(
      title = paste0("True graph:\n", paste0(capture.output(print(G)), collapse = "\n")),
      subtitle = sim_params_str
    ) +
    theme(
      title = element_text(size = 12)
    ) +
    theme_minimal()

  list(
    combined_post_prob = combined_post_prob,
    combined_AIC = combined_AIC
  )
}

set.seed(13)
G3_2_sim_results <- logdet_perm_fas_sim(G3_2)
G3_2_sim_results$sim_results <- process_sim_results(G3_2_sim_results, G3_2)
plot_G3_2_res <- plot_sim_results(G3_2_sim_results$sim_results, G3_2, max_configs = 50)

# TODO: Compare the all permn results to the backselect results
# HOPE is that the backselect results are the highest likelihood models
backselect_ll <- sapply(G3_2_sim_results$backselect_mods, function(x) x$log_lik)
all_permn_log_sum <- log_sum_exp(G3_2_sim_results$all_ordering$log_lik)

# Percentage is high! this is good
sum(exp(backselect_ll - all_permn_log_sum))

stop()
# system.time({
#   G3_sim <- multi_sim(10, G3, faster_init = FALSE)
#   G3_2_sim <- multi_sim(10, G3_2, faster_init = FALSE)
#   G3_3_sim <- multi_sim(10, G3_3, faster_init = FALSE)
# })

# set.seed(13)
# system.time({
#   G3_sim2 <- multi_sim(10, G3, faster_init = T)
#   G3_2_sim2 <- multi_sim(10, G3_2, faster_init = T)
#   G3_3_sim2 <- multi_sim(10, G3_3, faster_init = T)
# })

# stop()
G3_plot_res <- plot_n_sim(G3_sim, G3)
G3_2_plot_res <- plot_n_sim(G3_2_sim, G3_2)
G3_3_plot_res <- plot_n_sim(G3_3_sim, G3_3)

save_plot <- function(plot, suffix) {
  today_date <- format(Sys.Date(), "%Y-%m-%d")
  ggsave(print(sprintf("~/Projects/NESMR-simulations/results/%s_%s.pdf", today_date, suffix)),
         plot, width = 14, height = 14)
}

save_plot(G3_plot_res$combined_post_prob, "G3_10_post_prob")
save_plot(G3_plot_res$combined_AIC, "G3_10_AIC")
save_plot(G3_2_plot_res$combined_post_prob, "G3_2_10_post_prob")
save_plot(G3_2_plot_res$combined_AIC, "G3_2_10_AIC")
save_plot(G3_3_plot_res$combined_post_prob, "G3_3_10_post_prob")
save_plot(G3_3_plot_res$combined_AIC, "G3_3_10_AIC")

set.seed(13)
G4_sim_results <- logdet_perm_fas_sim(G4)
G4_sim_results$sim_results <- process_sim_results(G4_sim_results, G4)
plot_G4_res <- plot_sim_results(G4_sim_results$sim_results, G4, max_configs = 50)

ggsave(print(sprintf("~/Projects/NESMR-simulations/results/%s_four_node_mediation.png", today_date)),
       plot_G4_res$plot, width = 20, height = 8)

refit_mods_log_lik <- unlist(lapply(G4_sim_results$refit_mods, log_py))
