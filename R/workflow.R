# Run this file in its containing directory with 
# nohup nice -19 R CMD BATCH workflow.R &

source("setup.R")

#for(sample.size in sample.sizes)
#  simulator(sample.size = sample.size)

#pvals(which.datasets, ncpus, sample.sizes)
#unpack.stan()
#rocs()
#aucs(which.datasets, sample.sizes)

makePlots()