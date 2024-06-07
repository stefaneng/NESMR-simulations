library(GWASBrewer)

if (graph_type == "mediation") {
  G <- matrix(
    c(0, 0, 0,
      sqrt(0.4), 0, 0,
      0, sqrt(0.2), 0),
    nrow = 3,
    byrow = TRUE
  )
} else if (graph_type == "correlated") {
  G <- matrix(
    c(0, 0, 0,
      0, 0, 0,
      sqrt(0.4), sqrt(0.2), 0),
    nrow = 3,
    byrow = TRUE
  )
} else if (graph_type == "collider") {
  G <- matrix(
    c(0, 0, 0,
      sqrt(0.2), 0, 0,
      sqrt(0.4), 0, 0),
    nrow = 3,
    byrow = TRUE
  )
} else {
  stop("Invalid graph type")
}

G <- G * effect_scale_factor

#h2 <- 0.3
#J <- 5000
# N <- 30000, 2e5
#N <- 30000
#pi_J <- 0.1
# 1, 1e-3, 1e-5, 1e-8, 1e-10
#alpha <- 5e-8
#lambda <- qnorm(1 - alpha / 2)

dat <- GWASBrewer::sim_mv(
  G = G,
  N = N,
  J = J,
  h2 = h2,
  pi = pi_J,
  sporadic_pleiotropy = TRUE,
  est_s = TRUE
  )

true_beta <- extract_lower_triangular(G, "G_")
