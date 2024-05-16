DSC:
  define:
    analyze: winners_curse_nesmr
  run: simulate * analyze
  exec_path: R
  output: 2024-05-14_winners_curse

simulate: renv.R + winners_curse/sim_ma_winner_sim.R
  n: 40000
  J: 500000
  pi_J: 0.001
  scale_effects: 0.2, 0.5, 1
  h2: c(0.5, 0.3, 0.25)
  $N: n
  $G: G
  $dat: dat
  $seed: DSC_SEED

winners_curse_nesmr: renv.R + winners_curse/fit_ma_nesmr.R
  eta: 0.25, 0.5, 0.75
  alpha: 1e-7, 5e-8, 1e-9
  dat: $dat
  N: $N
  G: $G
  $cursed_model: cursed_model
  $cursed_ix: cursed_ix
  $ma_ix: ma_ix
  $ma_SNP_model: ma_SNP_model
  $seed: DSC_SEED