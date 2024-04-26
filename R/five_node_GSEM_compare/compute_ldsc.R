# Calculate LD scores and LDSC
dat <- GWASBrewer:::calc_ld_scores(dat, R_LD = ld_mat_list)

# Calculate ldsc matrices from Genomic SEM function
ldsc_res <- ldsc_GWASBrewer(
  dat, N
)
