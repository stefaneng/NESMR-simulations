renv::load('/nfs/turbo/sph-jvmorr/NESMR/simulations')
devtools::load_all('/nfs/turbo/sph-jvmorr/NESMR/esmr')
library(GWASBrewer)

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
data(ld_mat_list)
data(AF)
## simulate summary statistics
# Parameters from dsc
# N <- 40000
# J <- 5e5
#pi <- 500 / J
dat <- sim_mv(
  G = G,
  N = N,
  J = J,
  h2 = h2,
  pi = pi,
  R_LD = GWASBrewer::ld_mat_list,
  af = AF,
  sporadic_pleiotropy = FALSE,
  est_s = TRUE
)
