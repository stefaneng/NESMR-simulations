
library(ggplot2)
library(dscrutils)
library(tidyr)
library(plyr)
library(dplyr)

reticulate::use_condaenv('dsc')

G <- matrix(
  c(0, 0, 0, 0, 0,
    sqrt(0.3), 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, -1*sqrt(0.1), 0, 0, 0,
    0, -1*sqrt(0.1), sqrt(0.2), sqrt(0.25), 0),
  nrow = 5,
  byrow = 5
)

# Load your data
dscout <- dscquery(dsc.outdir = "/nfs/turbo/sph-jvmorr/NESMR/simulations/five_node_discovery",
                   targets    = c(
                    "discovery_algo.loglik_results", "discovery_algo.results_df", "discovery_algo.DSC_TIME",
                    "discovery_algo.backward_select_edges", "discovery_algo.discovery_model",
                    "discovery_algo.last_backward_mod", "discovery_algo.backward_select_adj_mat",
                    "discovery_algo.backward_results", "discovery_algo.backward_select_pvals",
                    "discovery_algo.backward_select_edges", "discovery_algo.backward_mod_results"
                    ))
