library(dplyr)
library(esmr)
library(igraph)

#threshold1 <- 0.1


minp <- apply(dat$pval_true[dat$ld_list_external, ], 1, min)
nesmr_ix <- dat$ld_list_external[minp < 5e-8]
M <- nrow(dat$direct_trait_effects)
## 1. esmr to generate super graph

mvmr_res_all <- purrr::map_dfr(1:5, function(i){
  cat(i, " ")
  mvmr_minp <- apply(dat$pval_true[dat$ld_list_external,-i], 1, min)
  mvmr_ix  <- which(minp < 5e-8) #which(mvmr_minp < 5e-8)
  myix <- dat$ld_list_external
  mvmr_res <- with(dat,
    esmr(
      beta_hat_X = beta_hat[myix,-i],
      se_X = s_estimate[myix,-i],
      beta_hat_Y = beta_hat[myix,i],
      se_Y = s_estimate[myix,i],
      variant_ix = mvmr_ix,
      max_iter = 300))
  pvals <- 2*pnorm(-abs(mvmr_res$beta$beta_m/mvmr_res$beta$beta_s))
  myres <- data.frame(exposure = (1:5)[-i],
                      outcome = i,
                      pval = pvals,
                      beta = mvmr_res$beta$beta_m,
                      se = mvmr_res$beta$beta_s, 
                      beta_true = dat$direct_trait_effects[-i,i])
  return(myres)
})

mvmr_res_all$padj <- p.adjust(mvmr_res_all$pval, method = "BH")

## Build adj matrix from MVMR results
mvmr_res_all$pval_weights <- pmin(-log10(mvmr_res_all$pval), 20)

wt_matrix <- reshape2::dcast(mvmr_res_all, exposure ~ outcome, value.var = "pval_weights")
wt_matrix <- as.matrix(wt_matrix[,-1])
wt_matrix[is.na(wt_matrix)] <- 0
wt_matrix[wt_matrix < -log10(threshold1)] <- 0

#mvmr_beta_edgelist <- mvmr_res_all %>%
#  filter(padj < threshold1) %>%
#  select(exposure, outcome, pval_weights)


mvmr_res_all <- mutate(mvmr_res_all, 
                       in_discovery = padj <= threshold1)

if(sum(mvmr_res_all$in_discovery) == 0){
    mvmr_res_all <- mutate(mvmr_res_all, 
                           in_adiscovery = NA,
                           disc_beta = NA, 
                           disc_pval = NA)

}else{


#B_disc <- wt_matrix
#B_disc[B_disc != 0] <- 1

#discovery_G <- graph_from_adjacency_matrix(wt_matrix)
wtd_discovery_G <- graph_from_adjacency_matrix(wt_matrix, weighted = TRUE)

# ## 1.5: Check if super graph is acyclic; Otherwise 2.
is_acyclic <- function(g) {
  tryCatch(
    !is.null(topo_sort(g)),
    error = function(e) FALSE)
}
if (! is_acyclic(wtd_discovery_G)) {
  ## 2. Perform feedback arc set with weights on p-values to reduce super graph to "best" acyclic subgraph
  fas <- feedback_arc_set(wtd_discovery_G, algo = "exact")
  print('Removing edges: ')
  print(fas)
  wtd_discovery_G <- wtd_discovery_G - fas
  print('New graph: ')
  print(wtd_discovery_G)
}

gg <- as_edgelist(wtd_discovery_G)
mvmr_res_all <- mutate(mvmr_res_all, 
                       in_adiscovery = paste0(exposure, "-", outcome) %in% paste0(gg[,1], "-", gg[,2]))


## 3. Fit discovery model

## Fit model on discovery set
discovery_adj_mat <- as.matrix(as_adjacency_matrix(wtd_discovery_G))

discovery_model <- with(dat,
                        esmr(
                          beta_hat_X = beta_hat,
                          se_X = s_estimate,
                          variant_ix = nesmr_ix,
                          G = diag(5), # required for network problem
                          direct_effect_template = discovery_adj_mat,
                          max_iter = 300))
disc_pval <- reshape2::melt(discovery_model$pvals_dm) %>%
                rename(exposure = Var1, outcome = Var2) %>%
                mutate(disc_pval = exp(value)) %>%
                select(exposure, outcome, disc_pval) %>%
                filter(exposure != outcome)

disc_beta <- reshape2::melt(discovery_model$direct_effects) %>%
                rename(exposure = Var1, outcome = Var2, disc_beta = value) %>% 
                filter(exposure != outcome)

mvmr_res_all <- full_join(mvmr_res_all, disc_beta) %>% full_join(., disc_pval)
mvmr_res_all$disc_pval[!mvmr_res_all$in_adiscovery] <- NA
mvmr_res_all$disc_beta[!mvmr_res_all$in_adiscovery] <- NA


}
