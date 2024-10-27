library(esmr)
library(dplyr)


minp <- apply(dat$pval_true[dat$ld_list_external, ], 1, min)
ix1 <- dat$ld_list_external[minp < 5e-8]

mvmr_result <- purrr::map_dfr(1:4, function(i){
  cat(i, " ")
  mvmr_minp <- apply(dat$pval_true[dat$ld_list_external,-(1:i), drop = F], 1, min)
  mvmr_ix  <- which(minp < 5e-8) #which(mvmr_minp < 5e-8)
  myix <- dat$ld_list_external
  mvmr_res <- with(dat,
    esmr(
      beta_hat_X = beta_hat[myix,-(1:i), drop = F],
      se_X = s_estimate[myix,-(1:i), drop = F],
      beta_hat_Y = beta_hat[myix,i],
      se_Y = s_estimate[myix,i],
      variant_ix = mvmr_ix,
      max_iter = 300))
  pvals <- 2*pnorm(-abs(mvmr_res$beta$beta_m/mvmr_res$beta$beta_s))
  myres <- data.frame(exposure = (1:5)[-(1:i)],
                      outcome = i,
                      pval = pvals,
                      beta_hat = mvmr_res$beta$beta_m,
                      se = mvmr_res$beta$beta_s, 
                      beta_true = dat$direct_trait_effects[-(1:i),i])
  return(myres)
})

