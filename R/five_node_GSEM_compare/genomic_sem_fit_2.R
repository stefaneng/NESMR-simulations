library(lavaan)
library(GenomicSEM)
library(dplyr)
library(stringr)

# TODO: Generalize this
trait_names <- paste0("V", 1:5)
sem_model <- "
    V1 ~ V2 + V3 + V4 + V5
    V2 ~ V3 + V4 + V5
    V3 ~ V4 + V5
    V4 ~ V5
  "

#sem_mod_ML <- sem(
#  model = sem_model,
#  sample.cov = ldsc_res$S,
#  sample.nobs = 200,
#  estimator = "ML"
#)
#
#W <- solve(ldsc_res$V)
#
#sem_mod_DWLS <- sem(
#  model = sem_model,
#  sample.cov = ldsc_res$S,
#  WLS.V = W,
#  sample.nobs = 2,
#  estimator = "DWLS",
#  std.lv = FALSE,
#  ordered = TRUE,
#  optim.dx.tol = +Inf
#)
#
genomic_sem_fit <- usermodel(
  covstruc = ldsc_res,
  model = sem_model
)


truth <- reshape2::melt(dat$direct_trait_effects) %>%
        rename(exposure = Var1, outcome = Var2, beta_true = value) %>%
        filter(exposure > outcome)

gsem_result <- genomic_sem_fit$results  %>%
        mutate(exposure = as.numeric(str_replace(rhs, "V", "")), 
               outcome = as.numeric(str_replace(lhs, "V", ""))) %>%
        rename(gsem_beta_hat = Unstand_Est, gsem_p = p_value) %>%
        select(exposure, outcome, gsem_beta_hat, gsem_p) %>%
        filter(exposure != outcome) %>%
        full_join(truth, .) 
