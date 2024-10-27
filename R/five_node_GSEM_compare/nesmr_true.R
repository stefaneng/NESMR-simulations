library(esmr)
library(dplyr)

B_lower <- lower.tri(G) + 0
minp <- apply(dat$pval_true[dat$ld_list_external, ], 1, min)
ix1 <- dat$ld_list_external[minp < 5e-8]
nesmr_model <- with(
                    dat,
                    esmr(beta_hat_X = beta_hat,
                         se_X = s_estimate,
                         variant_ix = ix1,
                         G = diag(5), # required for network problem
                         direct_effect_template = B_lower,
                         max_iter = 300))
nesmr_result <- reshape2::melt(nesmr_model$direct_effects) %>%
    rename(exposure = Var1, outcome = Var2, nesmr_beta_hat = value) %>%
    filter(exposure > outcome)
nesmr_p <- reshape2::melt(nesmr_model$pvals_dm) %>%
    mutate(nesmr_p = exp(value)) %>%
    rename(exposure = Var1, outcome = Var2) %>%
    select(-value) %>%
    filter(exposure > outcome)
truth <- reshape2::melt(dat$direct_trait_effects) %>%
    rename(exposure = Var1, outcome = Var2, beta_true = value) %>%
    filter(exposure > outcome)
nesmr_result <- nesmr_result %>% full_join(., nesmr_p) %>% full_join(truth, .)
