library(lavaan)
library(GenomicSEM)

trait_names <- paste0('V', 1:ncol(G))
sem_model <- '
    V1 ~ V2 + V3 + V4 + V5
    V2 ~ V3 + V4 + V5
    V3 ~ V4 + V5
    V4 ~ V5
  '

sem_mod_ML <- sem(
  model = sem_model,
  sample.cov = ldsc_res$S,
  sample.nobs = 200,
  estimator = "ML"
)

W <- solve(ldsc_res$V)

sem_mod_DWLS <- sem(
  model = sem_model,
  sample.cov = ldsc_res$S,
  WLS.V = W,
  sample.nobs = 2,
  estimator = "DWLS",
  std.lv = FALSE,
  ordered = TRUE,
  optim.dx.tol = +Inf
)

genomic_sem_fit <- usermodel(
  covstruc = ldsc_res
)
