compute.pvals = function(which.datasets = 1:100, ncpus = 2){
  group = readRDS("../data/group.rds")
  for(i in which.datasets){
    print(paste("Calculating p-values for dataset", i))
    cts = readRDS(paste("../simulations/sim", i, ".rds", sep=""))
    pv = pvals1dataset(cts, group, ncpus = ncpus)
    saveRDS("../pvals/pvals", i, ".rds", sep="")
  }
}