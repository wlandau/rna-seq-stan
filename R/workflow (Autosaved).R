# Run this file in its containing directory with 
# nohup nice -19 R CMD BATCH workflow.R &

source("setup.R")

#simulator(which.datasets)
pvals(which.datasets, ncpus)
rocs()
aucs()