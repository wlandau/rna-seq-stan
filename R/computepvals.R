compute.pvals = function(which.datasets = 1:100, ncpus = 2){
  group = readRDS("../data/group.rds")
  for(i in which.datasets){
    logfile("Calculating p-values for dataset", i)
    cts = readRDS(paste("../simulations/sim", i, ".rds", sep=""))
    pv = pvals1dataset(cts, group, ncpus = ncpus)
    saveRDS(pv, paste("../pvals/pvals", i, ".rds", sep=""))
  }
}