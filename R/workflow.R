# Run this file in its containing directory with 
# nohup nice -19 R CMD BATCH workflow.R &

source("setup.R")

#for(size in work.parms("size"))
#  simulator(size = size)

#unpack.stan()
ranks(mtds = c(    
#    "edgeR", 
#    "baySeq", 
#    "ShrinkBayes",
#    "ShrinkBayesMu",
#   "fullyBayes",
     "DESeq2"
), 
reps = reps.on.server, ncpus = ncpus)
#ranks(mtds = "ShrinkBayesMu", reps = 1, ncpus = 1)
#rocs(mtd = c("stan_corr", "stan", "stan_laplace"))
#aucs()
#fdr()
#makePlots()
#makePlots(exclude="stan_laplace")