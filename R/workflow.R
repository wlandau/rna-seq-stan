# Run this file in its containing directory with 
# nohup nice -19 R CMD BATCH workflow.R &

source("setup.R")

#for(size in work.parms("size"))
#  simulator(size = size)

unpack.stan()
#ranks(mtds = "ShrinkBayes", reps = reps.on.server, ncpus = ncpus)
ranks(mtds = "ShrinkBayes", reps = 2, ncpus = 1)
rocs(mtd = c("stan_corr", "stan", "stan_laplace"))
aucs()
fdr()
makePlots()
makePlots(exclude="stan_laplace")