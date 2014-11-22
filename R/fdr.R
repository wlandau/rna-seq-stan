prepFDRdf = function(fdr, upper = 0.15, npts = 100, pkg, sample.size, rep){
  apx = approx(x = fdr$cutoff, fdr$FDR, xout = (1:npts/npts)*0.15)
  apxfdr = data.frame(cutoff = apx$x, FDR = apx$y)

  apxfdr$FDRminusCutoff = apxfdr$FDR - apxfdr$cutoff
  apxfdr$pkg = pkg
  apxfdr$sample.size = sample.size
  apxfdr$rep = rep
  apxfdr = apxfdr[apxfdr$cutoff <= upper,]

  return(apxfdr)
}

fdr = function(pkgs = c("edgeR", "baySeq", "ShrinkBayes", "stan_corr", "stan", "stan_horseshoe"), sample.sizes = c(4, 8, 16), reps = 1:10, jarad = T, dan = T){

  for(pkg in pkgs)
    for(sample.size in sample.sizes)
      for(rep in reps){
        print(paste(pkg, sample.size, rep))

        d = data.frame(
          cutoff = readRDS(paste("../pvals/", pkg, rep, "-", sample.size, ".rds", sep="")),
          FDR = 1 - abs(readRDS(paste("../simulations/truth", rep, "-", sample.size, ".rds", sep=""))))
        d = d[order(d$cutoff),]

        if(dan){
          fdr = prepFDRdf(as.data.frame(apply(d, 2, function(x){cumsum(x)/1:length(x)})), 
                                     pkg = pkg, sample.size = sample.size, rep = rep)
          saveRDS(fdr, paste("../fdr/dan-", pkg, rep, "-", sample.size, ".rds", sep=""))
        }

        if(jarad){
          fdr = prepFDRdf(as.data.frame(t(sapply(1:(dim(d)[1] - 100), function(x){
            apply(d[x:(x + 100),], 2, mean)
          }))), pkg = pkg, sample.size = sample.size, rep = rep)

          saveRDS(fdr, paste("../fdr/jarad-", pkg, rep, "-", sample.size, ".rds", sep=""))
        }
      }

  dan.df = NULL
  jarad.df = NULL

  for(pkg in pkgs)
    for(sample.size in sample.sizes)
      for(rep in reps){
        print(paste(pkg, sample.size, rep))
        if(dan)
          dan.df = rbind(dan.df, readRDS(paste("../fdr/dan-", pkg, rep, "-", sample.size, ".rds", sep="")))
        if(jarad)
          jarad.df = rbind(jarad.df, readRDS(paste("../fdr/jarad-", pkg, rep, "-", sample.size, ".rds", sep="")))
      }

  if(dan)
    saveRDS(as.data.frame(dan.df), "../fdr/dan.rds")
  if(jarad)
    saveRDS(as.data.frame(jarad.df), "../fdr/jarad.rds")
}