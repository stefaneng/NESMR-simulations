library(dscrutils)
library(dplyr)
library(tidyr)
#reticulate::use_condaenv('dsc3')

dir <- "/nfs/turbo/sph-jvmorr/NESMR/simulations/2024-04-29_five_node"
dscout1 <- dscquery(dsc.outdir = dir,
                   targets    = c(
                    "nesmr_true.nesmr_result",
                    "simulate",
                    "simulate.x"),
                    #"genomic_sem.gsem_result",
                    #"nesmr_discovery.discovery_res"),
                    module.output.files = c("simulate"),
                   ignore.missing.files = TRUE)

r1 <- purrr::map_dfr(1:length(dscout1$DSC), function(i){
                             res <- dscout1$nesmr_true.nesmr_result[[i]]
                             res$DSC <- dscout1$DSC[i]
                             res$simfile <- dscout1$simulate.output.file[i]
                             res$simulate.x <- dscout1$simulate.x[i]
                             return(res)})
r1 <- r1 %>% rename(beta_hat = nesmr_beta_hat, pval = nesmr_p) %>%
    mutate(Method = "NESMR")

dscout2 <- dscquery(dsc.outdir = dir,
                   targets    = c(
                    "genomic_sem.gsem_result",
                    "simulate",
                    "simulate.x"),
                    #"genomic_sem.gsem_result",
                    #"nesmr_discovery.discovery_res"),
                    module.output.files = c("simulate"),
                   ignore.missing.files = TRUE)

r2 <- purrr::map_dfr(1:length(dscout2$DSC), function(i){
                             res <- dscout2$genomic_sem.gsem_result[[i]]
                             res$DSC <- dscout2$DSC[i]
                             res$simfile <- dscout2$simulate.output.file[i]
                             res$simulate.x <- dscout2$simulate.x[i]
                             return(res)})
r2 <- r2 %>% rename(beta_hat = gsem_beta_hat, pval = gsem_p) %>%
    mutate(Method = "GenomicSEM")

dscout3 <- dscquery(dsc.outdir = dir,
                   targets    = c(
                    "nesmr_discovery.discovery_res",
                    "nesmr_discovery.threshold1",
                    "simulate",
                    "simulate.x"),
                    #"genomic_sem.gsem_result",
                    #"nesmr_discovery.discovery_res"),
                    module.output.files = c("simulate"),
                   ignore.missing.files = TRUE)

thresh <- 0.05
r3 <- purrr::map_dfr(1:length(dscout3$DSC), function(i){
                             res <- dscout3$nesmr_discovery.discovery_res[[i]]
                             pbh <- p.adjust(res$disc_pval, method = "BH")
                             pbf <- p.adjust(res$disc_pval, method = "bonferroni")
                             pu <- res$disc_pval
                             r <- data.frame( DSC = dscout3$DSC[i],
                                             simfile = dscout3$simulate.output.file[i],
                                             simulate.x = dscout3$simulate.x[i],
                                             threshold = dscout3$nesmr_discovery.threshold1[i]) %>%
                                  mutate(ntrue_pbf = sum(pbf < thresh & !is.na(pbf) & res$beta_true != 0),
                                             nfalse_pbf = sum(pbf < thresh& !is.na(pbf) & res$beta_true == 0),
                                         nfalse_pu = sum(pu < thresh & !is.na(pu) & res$beta_true == 0),
                                         ntrue_pu = sum(pu < thresh & !is.na(pu) & res$beta_true != 0),
                                             ntrue_pbh = sum(pbh < thresh& !is.na(pbh) & res$beta_true != 0),
                                             nfalse_pbh = sum(pbh <thresh & !is.na(pbh) & res$beta_true == 0),
                                             real_fdr = nfalse_pbh/ntrue_pbh)
                             if(r$ntrue_pbh == 0) r$real_fdr = 0

                             return(r)})

r4 <- purrr::map_dfr(1:length(dscout3$DSC), function(i){
                             res <- dscout3$nesmr_discovery.discovery_res[[i]] %>%
                                 select(exposure, outcome, beta, pval, beta_true) %>%
                                 filter(exposure > outcome) %>%
                                 rename(beta_hat = beta) %>% mutate(Method = "MVMR_all")
                             res$DSC <- dscout3$DSC[i]
                             res$simfile <- dscout3$simulate.output.file[i]
                             res$simulate.x <- dscout3$simulate.x[i]
                             res$threshold <- dscout3$nesmr_discovery.threshold1[i]
                             return(res)}) %>% filter(threshold == 1) %>% select(-threshold)

dscout4 <- dscquery(dsc.outdir = dir,
                   targets    = c(
                    "mvmr_true.mvmr_result",
                    "simulate",
                    "simulate.x"),
                    #"genomic_sem.gsem_result",
                    #"nesmr_discovery.discovery_res"),
                    module.output.files = c("simulate"),
                   ignore.missing.files = TRUE)

r5 <- purrr::map_dfr(1:length(dscout4$DSC), function(i){
                             res <- dscout4$mvmr_true.mvmr_result[[i]]
                             res$DSC <- dscout4$DSC[i]
                             res$simfile <- dscout4$simulate.output.file[i]
                             res$simulate.x <- dscout4$simulate.x[i]
                             return(res)}) %>% select(-se)
r5 <- r5 %>%
    mutate(Method = "MVMR")

est_results <- bind_rows(r1, r2) %>% bind_rows(r4) %>% bind_rows(r5)

saveRDS(est_results, file = paste0(dir, "/est_results.RDS"))
saveRDS(r3, file = paste0(dir, "/disc_results.RDS"))
