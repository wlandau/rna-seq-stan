collectShrinkBayesMu = function(){
d = NULL
for(size in c(4, 8, 16)) for(rep in 1:10){
  file = paste("../priors/ShrinkBayesMu-", size, "-", rep, ".rds", sep="")
  if(file.exists(file)){
    x = readRDS(file)
    d = rbind(d, c(
      size = size,
      rep = rep,
      mu.hybrid = x$pmlist$mufixed,
      mu.parent1 = x$pmlist$mumu.parent1,
      mu.parent2 = x$pmlist$mumu.parent2,
      mu.dispersion = x$pmlist$mudisp,
      prec.hybrid = x$pmlist$precfixed,
      prec.parent1 = x$pmlist$precmu.parent1,
      prec.parent2 = x$pmlist$precmu.parent2,
      prec.dispersion = x$pmlist$precdisp
    ))
  }
}
  saveRDS(as.data.frame(d), "../priors/priors-ShrinkBayesMu.rds")
}