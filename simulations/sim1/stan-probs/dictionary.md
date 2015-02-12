# laplace_mvn

Each data frame (RDS file) in this directory corresponds to the simulated dataset with the given "size" and "rep" numbers. (File names are of the form, "pprob_lap_mvn<size>-<rep>"). Each row corresponds to a gene. The variables (columns) are:

#### geneid

ID of each gene (row): for example, "AC148152.3_FG001" 

#### phph

Probability of high parent heterosis for each gene. Probabilities were not conditioned on estimates of mean effects. 

#### plph

Same as phph, but for low parent heterosis.


# post_probs_horseshoe.rds

A single data frame (RDS). Rows correspond to a gene in a particular simulation. Only involved one sample.size (16) setting.

#### sample.size

Sample size setting for the simulation. 16. 

#### sim

ID for simulation within a sample.size. int type, one of 1:10. 

#### sim_gene

ID for gene within a simulation. one of 1:25000.

#### geneid

ID of each gene; for example, "AC148152.3_FG001"

#### phph

Posterior probability of HPH

#### plph

Posterior probability of LPH

# post_probs_laplace.rds
A single data frame (RDS). Rows correspond to a gene in a particular simulation.

#### sample.size

Sample size setting for the simulation. num type, one of 4, 8, 16. 

#### sim

ID for simulation within a sample.size. int type, one of 1:10. 

#### sim_gene

ID for gene within a simulation. one of 1:25000.

#### geneid

ID of each gene; for example, "AC148152.3_FG001"

#### phph

Posterior probability of HPH

#### plph

Posterior probability of LPH


# rev_probs_laplace.rds
A single data frame (RDS). Rows correspond to a gene in a particular simulation. This file was result of multivariate normal prior on mu1, mu2, mu3.

#### sample.size

Sample size setting for the simulation. num type, one of 4, 8, 16. 

#### sim

ID for simulation within a sample.size. int type, one of 1:10. 

#### geneid

ID of each gene; for example, "AC148152.3_FG001"

#### old_prob

max(plph,phph)

#### new_prob

if \hat{mu_3} \in (\mu_1,\mu_2), this is zero,
otherwise it is max(plph,phph)

# rev_probs_mu_cov.rds
A single data frame (RDS). Rows correspond to a gene in a particular simulation. This file was result of multivariate normal prior on mu1, mu2, mu3.

#### sample.size

Sample size setting for the simulation. num type, one of 4, 8, 16. 

#### sim

ID for simulation within a sample.size. int type, one of 1:10. 

#### geneid

ID of each gene; for example, "AC148152.3_FG001"

#### old_prob

max(plph,phph)

#### new_prob

if \hat{mu_3} \in (\mu_1,\mu_2), this is zero,
otherwise it is max(plph,phph)

# revised_pvals.rds

A single data frame (RDS). Rows correspond to a gene in a particular simulation. This file was result of independent normal priors on mu1, mu2, mu3.

#### sample.size

Sample size setting for the simulation. num type, one of 4, 8, 16. 

#### sim

ID for simulation within a sample.size. int type, one of 1:10. 

#### geneid

ID of each gene; for example, "AC148152.3_FG001"

#### V1

if \hat{mu_3} \in (\mu_1,\mu_2), this is zero,
otherwise it is max(plph,phph)


# stan_probs_corr.rds
A single data frame (RDS). Rows correspond to a gene in a particular simulation.

#### sample.size

Sample size setting for the simulation. num type, one of 4, 8, 16. 

#### sim

ID for simulation within a sample.size. int type, one of 1:10. 

#### sim_gene

ID for gene within a simulation. one of 1:25000.

#### geneid

ID of each gene; for example, "AC148152.3_FG001"

#### phph

Posterior probability of HPH

#### plph

Posterior probability of LPH


# stan_probs

A single data frame (RDS). Rows correspond to a gene in a particular simulation.

#### sample.size

Sample size setting for the simulation. num type, one of 4, 8, 16. 

#### sim

ID for simulation within a sample.size. int type, one of 1:10. 

#### sim_gene

ID for gene within a simulation. one of 1:25000.

#### geneid

ID of each gene; for example, "AC148152.3_FG001"

#### phph

Posterior probability of HPH

#### plph

Posterior probability of LPH
