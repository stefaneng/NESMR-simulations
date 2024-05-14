#nohup dsc --replicate 100  --host config.yml -c 4 xxxxx.dsc > xxxxx.out &

DSC:
  define:
    analyze: winners_curse_nesmr
  run: simulate * analyze
  exec_path: R
  output: 2024-05-13_winners_curse

simulate: renv.R + winners_curse/sim_ma_winner_sim.R
  n: 40000
  J: 500000
  pi_J: 0.001
  scale_effects: 0.1, 0.5, 1
  h2: c(0.5, 0.3, 0.25)
  $N: n
  $G: G
  $dat: dat
  $seed: DSC_SEED

winners_curse_nesmr: renv.R + winners_curse/fit_ma_nesmr.R
  eta: 0.1, 0.25, 0.5, 0.7, 0.9
  alpha: 1e-6, 1e-7, 5e-8
  dat: $dat
  N: $N
  G: $G
  $cursed_model: cursed_model
  $cursed_ix: cursed_ix
  $ma_ix: ma_ix
  $ma_SNP_model: ma_SNP_model
  $seed: DSC_SEED