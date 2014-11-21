unpack.stan = function(check.names = F){
  for(file in c("../pvals/stan_probs.rds", "../pvals/stan_probs_corr.rds")){

    sp = readRDS(file)

     if(file == "../pvals/stan_probs.rds")
       mtd = "stan"
     else if(file == "../pvals/stan_probs_corr.rds")
       mtd = "stan_corr"

    ddply(sp, .variables = .(sample.size, sim), .fun = function(x){
      pv = 1 - (x$phph + x$plph)
      names(pv) = x$geneid
      saveRDS(pv, paste("../pvals/", mtd, x$sim[1], "-", x$sample.size[1], ".rds", sep=""))
    })

    # check if names agree
    if(check.names){
      for(sample.size in c(4, 8, 16)){
        for(sim in 1:10){
          pv = readRDS(paste("../pvals/", mtd, sim, "-", sample.size, ".rds", sep=""))
          tr = readRDS(paste("../simulations/", "truth", sim, "-", sample.size, ".rds", sep=""))
          ct = readRDS(paste("../simulations/", "truth", sim, "-", sample.size, ".rds", sep=""))
          print(all(names(pv) == names(tr)))
          print(all(names(pv) == names(ct)))
          print(all(names(tr) == names(ct)))
        }
      }
    }
  }
}