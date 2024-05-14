#renv::load('/nfs/turbo/sph-jvmorr/NESMR/simulations')

G <- matrix(
  c(0, 0, 0,
    sqrt(0.4), 0, 0,
    0, sqrt(0.2), 0),
  nrow = 3,
  byrow = TRUE
)

G <- scale_effects * G

B_true <- G
B_true[!B_true == 0] <- 1
B_lower <- lower.tri(B_true) + 0

dat <- GWASBrewer::sim_mv(
  G = G,
  N = n,
  J = J,
  h2 = h2,
  pi = pi_J,
  sporadic_pleiotropy = FALSE,
  est_s = TRUE
)