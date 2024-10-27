DSC:
  run: simulate
  exec_path: R
  output: 2024-06-04_uvmr_ma

simulate: renv.R + uvmr_ma_compare/uvmr_ma_sim.R + R(
  res <- uvmr_rand_wc_sim(
    n = n,
    J = J,
    true_beta = true_beta,
    pi_J = pi_J,
    h2 = h2))
  n: 10000, 20000, 40000
  J: 500000
  pi_J: 0.001
  h2: 0.25, 0.5
  true_beta: 0.1
  $res: res
