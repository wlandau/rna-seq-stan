work.parms = function(item){
  mtd = c("edgeR", "baySeq", "stan_corr", "ShrinkBayes", "stan", "stan_laplace", "stan_laplace_mvn", "ShrinkBayesMu", "pprob_lap_phimix", "rev_probs_laplace", "rev_probs_mu_cov", "revised_pvals")
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

unpack.stan = function(mtds){
  for(mtd in intersect(mtds, c("stan_corr", "stan", "stan_laplace"))){
    file = switch(mtd,
                        "stan_corr" = "../ranks/stan_probs_corr.rds",
                        "stan" = "../ranks/stan_probs.rds",
                        "stan_laplace" = "../ranks/post_probs_laplace.rds")

    all.probs = readRDS(file)

    ddply(all.probs, .variables = .(sample.size, sim), .fun = function(x){

#      probs = 1 - pmax(x$phph, x$plph)
      probs = 1 - (x$phph + x$plph)
     
      names(probs) = x$geneid
      saveRDS(probs, paste("../ranks/", file.name(mtd, x$sample.size[1], x$sim[1]), sep=""))
    })
  }

  for(mtd in intersect(mtds, c("rev_probs_laplace", "rev_probs_mu_cov", "revised_pvals"))){
    file = switch(mtd,
                        "rev_probs_laplace" = "../other_Eric_posterior_probs/rev_probs_laplace.rds", 
                        "rev_probs_mu_cov" = "../other_Eric_posterior_probs/rev_probs_mu_cov.rds", 
                        "revised_pvals" = "../other_Eric_posterior_probs/revised_pvals.rds")

    all.probs = readRDS(file)
    if(mtd == "revised_pvals")
      all.probs$new_prob = all.probs$V1

    ddply(all.probs, .variables = .(sample.size, sim), .fun = function(x){

#      probs = 1 - pmax(x$phph, x$plph)
#      probs = 1 - (x$phph + x$plph)
    probs = 1- x$new_prob
     
      names(probs) = x$geneid
      saveRDS(probs, paste("../ranks/", file.name(mtd, x$sample.size[1], x$sim[1]), sep=""))
    })
  }


  if(mtd == "stan_laplace_mvn"){
    file.begin = "../laplace_mvn_posterior_probs/pprob_lap_mvn"
    for(size in work.parms("size")){
      for(rep in work.parms("rep")){
        x = readRDS(paste(file.begin, size, "-", rep, sep=""))
#       probs = 1 - pmax(x$phph, x$plph)
        probs = 1 - (x$phph + x$plph)
        saveRDS(probs, paste("../ranks/", file.name("stan_laplace_mvn", size, rep), sep=""))
      }
    }
  }

  if(mtd == "pprob_lap_phimix"){
    file.begin = "../laplace_mvn_posterior_probs/pprob_lap_phimix"
    for(size in work.parms("size")){
      for(rep in work.parms("rep")){
        x = readRDS(paste(file.begin, size, "-", rep, sep=""))
#       probs = 1 - pmax(x$phph, x$plph)
        probs = 1 - (x$phph + x$plph)
        saveRDS(probs, paste("../ranks/", file.name("pprob_lap_phimix", size, rep), sep=""))
      }
    }
  }
}