work.parms = function(item){
  mtd = c("edgeR", "ShrinkBayes", "baySeq", "stan_corr", "stan", "stan_laplace")
  switch(item, 
             "mtd" = ordered(mtd, levels = mtd),
             "size" = c(4, 8, 16),
             "rep" = 1:10,
             "keep" = 25000)
}

file.name = function(mtd, size, rep){
  paste(mtd, "-", size, "-", rep, ".rds", sep="")
}

loopify = function(fun, mtd = work.parms("mtd"), size = work.parms("size"), rep = work.parms("rep")){
  .grid = expand.grid(rep, size, mtd)
  colnames(.grid) = c("rep", "size", "mtd")
  as.data.frame(t(Vectorize(fun, c("mtd", "size", "rep"))(.grid$mtd, .grid$size, .grid$rep)))
}

unpack.stan = function(mtds = c("stan_corr", "stan", "stan_laplace"), check.names = F){
  for(mtd in mtds){
    file = switch(mtd,
                        "stan_corr" = "../ranks/stan_probs_corr.rds",
                        "stan" = "../ranks/stan_probs.rds",
                        "stan_laplace" = "../ranks/post_probs_laplace.rds")

    all.probs = readRDS(file)

    ddply(all.probs, .variables = .(sample.size, sim), .fun = function(x){
      probs = 1 - (x$phph + x$plph)
      names(probs) = x$geneid
      saveRDS(probs, paste("../ranks/", file.name(mtd, x$sample.size[1], x$sim[1]), sep=""))
    })
  }
}