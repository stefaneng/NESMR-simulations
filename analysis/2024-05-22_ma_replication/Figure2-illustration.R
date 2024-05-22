library(ggpubr)
library(ggplot2)
library(grid)
bias.carve = NULL
sd.carve = NULL

bias.select = NULL
sd.select = NULL

one_sim <- function(eta, pvalue, lambda, gamma, se.gamma.hat, N = 1e6, bins = 100) {
  xlimits <- c(-1, 1) * (gamma + 4 * se.gamma.hat)
  if (abs(gamma) > 0.02) {
    xlimits <- gamma + c(-1, 1) * 3 * se.gamma.hat
  }
  Z = rnorm(Nrep, mean = 0, sd = eta)

  gamma.hat = rnorm(Nrep, mean = gamma, sd = se.gamma.hat)

  ############ Carving
  ### Selection variable
  S = abs(gamma.hat/se.gamma.hat + Z) - lambda

  ### Number of selection
  mean(S > 0)

  ### No selection mean
  cat('No selection mean =', mean(gamma.hat) - gamma, '\n')

  ### Post-selection mean
  bias.select = (mean(gamma.hat[S >0]) - gamma)/sd((gamma.hat[S >0]))
  sd.select = sd((gamma.hat[S >0]))
  ###
  selected.set = S > 0
  selected.set.org = (abs(gamma.hat/se.gamma.hat) - lambda > 0)

  trun.center2 = (lambda - gamma.hat/se.gamma.hat)/eta
  trun.center1 = (- lambda - gamma.hat/se.gamma.hat)/eta
  gamma.carve =se.gamma.hat*( gamma.hat/se.gamma.hat - (1/eta) * ( - dnorm(trun.center1) + dnorm(trun.center2))/
                                (1 - pnorm(trun.center2) +  pnorm(trun.center1) ))

  trun.center2.plugin = (lambda - gamma.hat/se.gamma.hat)
  trun.center1.plugin = (- lambda - gamma.hat/se.gamma.hat)
  gamma.plugin = gamma.hat - se.gamma.hat*( - dnorm(trun.center1) + dnorm(trun.center2))/(1 - pnorm(trun.center2) +  pnorm(trun.center1) )



  ### Standized bias after carving correction
  bias.carve = (mean(gamma.carve[selected.set])- gamma)/sd(gamma.carve[selected.set])
  sd.carve = sd(gamma.carve[selected.set])
  # hist( (gamma.carve[selected.set] - gamma)/sd(gamma.carve[selected.set]))

  data.gamma1 = data.frame(
    Method = factor(c(rep("Rao-Blackwellization", sum(selected.set)), rep("No correction", sum(selected.set.org)))),
    Estimate = c(gamma.carve[selected.set],
                 gamma.hat[selected.set.org]))

  # Create a bquote that include eta and gamma
  annot_grob <- bquote(eta == .(eta) ~ gamma == .(round(gamma, 4)))

  grob <- grobTree(textGrob(annot_grob, x=0.1,  y=0.95, hjust=0,
                            gp=gpar(fontsize=28, fontface="bold")))

  h2 = gghistogram(
    data.gamma1, x = "Estimate", y = "..density..",
    add = "mean",
    color = "Method", palette = c("#00AFBB", "#E7B800"),
    fill = "Method",
    add_density = TRUE, bins = bins,
    ylab = "Density",
    xlim = xlimits,
    #xlab = expression(paste("selected ",hat(beta))),
    #title = expression(paste("(B). Case when ",beta[x[j]]/sigma[x[j]], " = " , lambda))
  ) + geom_vline(xintercept = gamma, linetype="dotted",
                 color = "red", size=1.5)+
    theme(axis.title = element_text(face="bold", size = rel(1.4)),
          axis.text = element_text( size = rel(1.5)),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0, unit='pt')),
          legend.position = "bottom",
          legend.title = element_text(face = "bold", margin = margin(r = 10, unit = 'pt'), size = rel(1.4)),
          legend.text = element_text(margin = margin(r = 10, unit = 'pt'), size = rel(1.4) ) ) +
    annotation_custom(grob)

  list(
    plot = h2,
    plot_data = data.gamma1
  )
}

Nrep = 1e7
eta = c(0.25, 0.5, 1, 1.25)
pvalue = 5*10^(-5)
lambda = abs(qnorm( 1-pvalue/2 ))
se.gamma.hat = sqrt(1/1e5)

gamma = c(
  lambda/10 * se.gamma.hat,
  lambda * se.gamma.hat,
  4*lambda * se.gamma.hat
  )

# params <- expand.grid(gamma = gamma, eta = eta)

# plot_results <- apply(params, 1, function(x) {
#   one_sim(eta = x['eta'], pvalue, lambda, gamma = x['gamma'], se.gamma.hat, N = Nrep)$plot
# })

separate_eta_plots <- lapply(eta, function(x) {
  p <- lapply(gamma, function(y) {
    one_sim(eta = x, pvalue, lambda, gamma = y, se.gamma.hat, N = Nrep)$plot
  })

  call_args <- c(
    p, common.legend = TRUE, legend = "bottom", ncol = 3, nrow = 1
  )

  do.call('ggarrange', call_args)
})

call_args <- c(
  separate_eta_plots, common.legend = TRUE, legend = "bottom", ncol = 1, nrow = length(eta)
)
merged_plots <- do.call('ggarrange', call_args)

ggsave(
  "RB-illustration_eta.pdf",
  merged_plots,
  width = 20,
  height = 30,
)
