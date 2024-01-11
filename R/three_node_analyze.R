renv::load('/nfs/turbo/sph-jvmorr/NESMR/simulations')
devtools::load_all('/nfs/turbo/sph-jvmorr/NESMR/esmr')

true_ll <- logLik(true_model)
full_ll <- logLik(full_model)

lrt_pvalue <- pchisq(
    - 2 * (true_ll - full_ll),
    df = sum(B_alt1) - sum(B_true),
    lower.tail = FALSE
)

true_fbar <- true_model$f$fbar
full_fbar <- full_model$f$fbar