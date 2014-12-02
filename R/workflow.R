# Run this file in its containing directory with 
# nohup nice -19 R CMD BATCH workflow.R &

source("setup.R")

#for(size in work.parms("size"))
#  simulator(size = size)

#unpack.stan()
ranks(mtds = "ShrinkBayes", reps = reps.on.server)
#rocs()
#aucs()
#fdr()
#makePlots()