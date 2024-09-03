## Goal of this simulation is to compare the logdet method to the permutation method
# The first test is to have four nodes, but only two nodes are connected
library(esmr)
library(dplyr)


## Two node cases: with heritable confounder U
effect_size_scale <- 1

G <- matrix(
  c(0, 0, 0, 0,
    0.1, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0),
  nrow = 4,
  byrow = TRUE
) * effect_size_scale

traits <- c("Y", "X", "U1", "U2")

rownames(G) <- colnames(G) <- traits

d <- nrow(G)

# Full network without the diagonal
B_full <- matrix(1, nrow = d, ncol = d)
B_lower <- lower.tri(G) + 0
diag(B_full) <- 0

h2 <- 0.3
J <- 5000
N <- 20000
pi_J <- 0.1
alpha <- 5e-8
lambda <- qnorm(1 - alpha / 2)

if (FALSE) {
  res <- replicate(1000, {
    set.seed(NULL)
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

    issue_dat <- NULL
    error_msg <- NULL

    rtn_list <- list()

    # Fit MVMR matrix
    tryCatch({
      rtn_list$mvmr_matrix <- with(dat, nesmr_complete_mvmr(
        beta_hat = beta_hat,
        se_beta_hat = s_estimate,
        pval_select = pval_true
      ))
    }, error = function(e) {
      rtn_list$error_fct <- "n_mvmr"
      rtn_list$issue_dat <- dat
      rtn_list$error_msg <- e
    })

    tryCatch({
      rtn_list$nesmr_full <- with(dat, nesmr_complete(
        beta_hat = beta_hat,
        se_beta_hat = s_estimate,
        pval_select = pval_true
      ))
    }, error = function(e) {
      rtn_list$error_fct <- "nesmr"
      rtn_list$issue_dat <- dat
      rtn_list$error_msg <- e
    })

    return(rtn_list)
  }, simplify = FALSE)

  ## TODO: Histogram of all of the beta estimates
  edge_results <- bind_rows(
    lapply(
      res,
      function(x) {
        mvmr_beta_df <- matrix_to_edgelist(
          x$mvmr_matrix$beta, lower_tri = FALSE, remove_diag = TRUE,
          value = 'beta')

        mvmr_se_df <- matrix_to_edgelist(
          x$mvmr_matrix$se, lower_tri = FALSE, remove_diag = TRUE,
          value = 'se'
          )

        mvmr_df <- merge(
          mvmr_beta_df, mvmr_se_df, by = c("from", "to"))

        complete_beta_df <- matrix_to_edgelist(
          x$nesmr_full$beta, lower_tri = FALSE, remove_diag = TRUE,
          value = 'beta')

        complete_se_df <- matrix_to_edgelist(
          x$nesmr_full$se, lower_tri = FALSE, remove_diag = TRUE,
          value = 'se')

        complete_df <- merge(
          complete_beta_df, complete_se_df, by = c("from", "to"))

        bind_rows(
          cbind(mvmr_df, model = "MVMR-n"),
          cbind(complete_df, model = "NESMR-all")
        )
      }
    ),
    .id = "rep"
  )

  edge_results$pval <- qnorm(.975) * pnorm(-abs(edge_results$beta / edge_results$se))
  edge_results$reject <- edge_results$pval < 0.05
}

saveRDS(edge_results, '~/Projects/NESMR-simulations/data/2024-09-03_mvmr_vs_complete_matrix.rds')

library(ggplot2)
# Create a ggplot histogram facet with rows as "from" and columns from "to"
# The facet should be grouped/colored by "model"
# The histogram should be the beta
ggplot(edge_results, aes(x = beta, fill = model)) +
  geom_histogram(bins = 30) +
  facet_grid(from ~ to) +
  theme_minimal()

edge_results[
  edge_results$from != 2 & edge_results$to != 1,
] %>%
  group_by(from, to, model) %>%
  summarise(
    mean_beta = mean(beta, na.rm = TRUE),
    sd_beta = sd(beta, na.rm = TRUE),
    mean_reject = mean(reject, na.rm = TRUE),
    na_percent = sum(is.na(beta))
  ) %>% View


# Compute difference in bias between MVMR and NESMR for each edge
# Long to wide then subtract the columns
library(tidyr)
edge_results_wide <- pivot_wider(
  edge_results, names_from = model, values_from = c(beta, se), id_cols = c('rep', 'from', 'to'))
head(edge_results)
head(edge_results_wide)

true_edges <- matrix_to_edgelist(
  unname(G), lower_tri = FALSE, remove_diag = TRUE,
  value = 'true_beta')

merge(
  edge_results_wide,
  true_edges,
  by = c('from', 'to')
) %>%
  mutate(
    bias_MVMR_n = `beta_MVMR-n` - true_beta,
    bias_NESMR_all = `beta_NESMR-all` - true_beta,
    null_edge = true_beta == 0
  ) %>%
  group_by(
    from, to, null_edge
  ) %>% summarise(
    mean_bias_MVMR_n = mean(bias_MVMR_n, na.rm = TRUE),
    mean_bias_NESMR_all = mean(bias_NESMR_all, na.rm = TRUE),
    bias_diff = abs(mean_bias_MVMR_n) - abs(mean_bias_NESMR_all),
  ) %>%
  View


## Seems to be slightly underestimating the standard error

# bootstrap_reps <- 100
#
# bootstrap_res <- replicate(bootstrap_reps, {
#   bootstrap_tot_est <- total_est + rnorm(nrow(total_est), 0, total_est_se)
#   dag_proj <- project_to_DAG(
#     bootstrap_tot_est, s = max(bootstrap_tot_est^2) + 0.1, penalty = "L2",
#     zap_to_DAG = FALSE
#   )
#   dag_proj
# }, simplify = FALSE)
#
# do.call('rbind.data.frame', bootstrap_res)
