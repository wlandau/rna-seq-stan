fdr = function(pkgs = c("edgeR", "baySeq", "ShrinkBayes", "stan"), sample.sizes = c(4, 8, 16), reps = 1:10){

  dan = NULL
  jarad = NULL

  for(pkg in pkgs)
    for(sample.size in sample.sizes)
      for(rep in reps){
        print(paste(pkg, sample.size, rep))

        d = data.frame(
          cutoff = readRDS(paste("../pvals/", pkg, rep, "-", sample.size, ".rds", sep="")),
          FDR = 1 - abs(readRDS(paste("../simulations/truth", rep, "-", sample.size, ".rds", sep=""))))
        d = d[order(d$cutoff),]

        fdr = as.data.frame(apply(d, 2, function(x){cumsum(x)/1:length(x)}))
        fdr$FDRminusCutoff = fdr$FDR - fdr$cutoff
        fdr$pkg = pkg
        fdr$sample.size = sample.size
        fdr$rep = rep
        saveRDS(fdr, paste("../fdr/dan-", pkg, rep, "-", sample.size, ".rds", sep=""))
      
        fdr = as.data.frame(t(sapply(1:(dim(d)[1] - 100), function(x){
          apply(d[x:(x + 100),], 2, mean)
        })))

        fdr$FDRminusCutoff = fdr$FDR - fdr$cutoff
        fdr$pkg = pkg
        fdr$sample.size = sample.size
        fdr$rep = rep
        saveRDS(fdr, paste("../fdr/jarad-", pkg, rep, "-", sample.size, ".rds", sep=""))
      }

  dan = NULL
  jarad = NULL

  for(pkg in pkgs)
    for(sample.size in sample.sizes)
      for(rep in reps){
        print(paste(pkg, sample.size, rep))
        dan = rbind(dan, readRDS(paste("../fdr/dan-", pkg, rep, "-", sample.size, ".rds", sep="")))
        jarad = rbind(jarad, readRDS(paste("../fdr/jarad-", pkg, rep, "-", sample.size, ".rds", sep="")))
      }

  saveRDS(as.data.frame(dan), "../fdr/dan.rds")
  saveRDS(as.data.frame(jarad), "../fdr/jarad.rds")
}