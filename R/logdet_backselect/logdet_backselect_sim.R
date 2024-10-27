library(esmr)
library(NESMR.Sims)

## Input:
# G
# h2
# J = 5000
# N = 20000
# pi = 0.1 (~ 500 causal variants)
# alpha = 5e-8
# permute = TRUE
# seed ?
#
# Output: four_node_results
#

sim_results <- logdet_perm_fas_sim(
    G = G,
    alpha = alpha,
    include_all_perms = ncol(G) <= 4,
    gwas_params = list(
        h2 = h2,
        J = J,
        N = N,
        pi = 0.1,
        sporadic_pleiotropy = TRUE
        )
    )

summary_sim_results <- logdet_process_sim_results(
    sim_results,
    G = G)