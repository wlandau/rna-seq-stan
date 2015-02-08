# laplace_mvn

Each data frame (RDS file) in this directory corresponds to the simulated dataset with the given "size" and "rep" numbers. (File names are of the form, "pprob_lap_mvn<size>-<rep>"). Each row corresponds to a gene. The variables (columns) are:

#### geneid

ID of each gene (row): for example, "AC148152.3_FG001" 

#### phph

Probability of high parent heterosis for each gene. Probabilities were not conditioned on estimates of mean effects. 

#### plph

Same as phph, but for low parent heterosis.


# post_probs_horseshoe.rds

# post_probs_laplace.rds

# post_probs.rds

# rev_probs_laplace.rds

# rev_probs_mu_cov.rds

# revised_pvals.rds

# stan_probs_corr.rds

# stan_probs.rds
