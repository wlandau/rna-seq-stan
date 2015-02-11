collectShrinkBayesPhiAlpDel = function(){
d = NULL
for(size in c(4, 8, 16)) for(rep in 1:10){
  file = paste("../priors/ShrinkBayes-", size, "-", rep, ".rds", sep="")
  if(file.exists(file)){
    x = readRDS(file)
    d = rbind(d, c(
      size = size,
      rep = rep,
      mean.phi = x$pmlist$mufixed,
      mean.alpha = x$pmlist$mualp,
      mean.delta = x$pmlist$mudel,
      mean.dispersion = x$pmlist$mudisp,
      prec.phi = x$pmlist$precfixed,
      prec.alpha = x$pmlist$precalp,
      prec.delta = x$pmlist$precdel,
      prec.dispersion = x$pmlist$precdisp
    ))
  }
}
  saveRDS(as.data.frame(d), "../priors/priors-ShrinkBayes-phi-alpha-delta.rds")
}

collectShrinkBayesMu = function(){
d = NULL
for(size in c(4, 8, 16)) for(rep in 1:10){
  file = paste("../priors/ShrinkBayesMu-", size, "-", rep, ".rds", sep="")
  if(file.exists(file)){
    x = readRDS(file)
    d = rbind(d, c(
      size = size,
      rep = rep,
      mu.phi = x$pmlist$mufixed,
      mu.alpha = x$pmlist$mumu.alpha,
      mu.delta = x$pmlist$mumu.delta,
      mu.dispersion = x$pmlist$mudisp,
      prec.phi = x$pmlist$precfixed,
      prec.alpha = x$pmlist$precmu.alpha,
      prec.delta = x$pmlist$precmu.delta,
      prec.dispersion = x$pmlist$precdisp
    ))
  }
}
  saveRDS(as.data.frame(d), "../priors/priors-ShrinkBayesMu.rds")
}