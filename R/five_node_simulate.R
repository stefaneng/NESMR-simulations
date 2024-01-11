renv::load('/nfs/turbo/sph-jvmorr/NESMR/simulations')
devtools::load_all('/nfs/turbo/sph-jvmorr/NESMR/esmr')
library(GWASBrewer)
library(igraph)
library(dplyr)

threshold <- 0.05

is_acyclic <- function(g) {
  tryCatch(
    !is.null(topo_sort(g)),
    error = function(e) FALSE)
}

G <- matrix(
  c(0, 0, 0, 0, 0,
    sqrt(0.3), 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, -1 * sqrt(0.1), 0, 0, 0,
    0, -1 * sqrt(0.1), sqrt(0.2), sqrt(0.25), 0),
  nrow = 5,
  byrow = 5
)

h2 <- c(0.5, 0.3, 0.25, 0.4, 0.3)
## simulate summary statistics
data(ld_mat_list)
data(AF)
dat <- sim_mv(
  G = G,
  N = 40000,
  J = 5e5,
  h2 = h2,
  pi = 500/5e5,
  R_LD = ld_mat_list,
  af = AF,
  est_s = TRUE
)

Z <- with(dat, beta_hat/s_estimate)
dat$pval <- 2*pnorm(-abs(Z))
minp <- apply(dat$pval, 1, min)

# ld pruning
dat$ld_list_minp <- sim_ld_prune(dat, R_LD = ld_mat_list, pvalue = minp)