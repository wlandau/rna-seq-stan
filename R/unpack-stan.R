laplace_mvn1 = function(){
  file.begin = paste(work.parms("path"), "stan-probs/laplace_mvn/pprob_lap_mvn", sep="")
  for(size in work.parms("size")){
    for(rep in work.parms("rep")){
      x = readRDS(paste(file.begin, size, "-", rep, sep=""))
#     probs = 1 - pmax(x$phph, x$plph)
      probs = 1 - (x$phph + x$plph)
      saveRDS(probs, paste(work.parms("path"), "ranks/", file.name("stan_laplace_mvn", size, rep), sep=""))
    }
  }
}

# only has sample.size = 16
post_probs_horseshoe.rds1 = function(){
  d = readRDS(paste(work.parms("path"), "stan-probs/post_probs_horseshoe.rds", sep=""))
  ddply(d, .variables = .(sample.size, sim), .fun = function(x){
#    probs = 1 - pmax(x$phph, x$plph)
    probs = 1 - (x$phph + x$plph)   
    names(probs) = x$geneid
    saveRDS(probs, paste(work.parms("path"), "ranks/", file.name("stan_horseshoe", x$sample.size[1], x$sim[1]), sep=""))
  })  
}

post_probs_laplace.rds1 = function(){
  d = readRDS(paste(work.parms("path"), "stan-probs/post_probs_laplace.rds", sep=""))
  ddply(d, .variables = .(sample.size, sim), .fun = function(x){
#    probs = 1 - pmax(x$phph, x$plph)
    probs = 1 - (x$phph + x$plph)   
    names(probs) = x$geneid
    saveRDS(probs, paste(work.parms("path"), "ranks/", file.name("stan_laplace", x$sample.size[1], x$sim[1]), sep=""))
  })  
}

rev_probs_laplace.rds1 = function(){
  d = readRDS(paste(work.parms("path"), "stan-probs/rev_probs_laplace.rds", sep=""))
  ddply(d, .variables = .(sample.size, sim), .fun = function(x){
#    probs = 1 - pmax(x$phph, x$plph)
#    probs = 1 - (x$phph + x$plph)   
    probs = 1 - x$new_prob   
    names(probs) = x$geneid
    saveRDS(probs, paste(work.parms("path"), "ranks/", file.name("rev_laplace", x$sample.size[1], x$sim[1]), sep=""))
  })  
}

rev_probs_mu_cov.rds1 = function(){
  d = readRDS(paste(work.parms("path"), "stan-probs/rev_probs_mu_cov.rds", sep=""))
  ddply(d, .variables = .(sample.size, sim), .fun = function(x){
#    probs = 1 - pmax(x$phph, x$plph)
#    probs = 1 - (x$phph + x$plph)   
    probs = 1 - x$new_prob
    names(probs) = x$geneid
    saveRDS(probs, paste(work.parms("path"), "ranks/", file.name("rev_probs_mu_cov", x$sample.size[1], x$sim[1]), sep=""))
  })  
}

revised_pvals.rds1 = function(){
  d = readRDS(paste(work.parms("path"), "stan-probs/revised_pvals.rds", sep=""))
  ddply(d, .variables = .(sample.size, sim), .fun = function(x){
#    probs = 1 - pmax(x$phph, x$plph)
#    probs = 1 - (x$phph + x$plph)   
    probs = 1 - x$V1
    names(probs) = x$geneid
    saveRDS(probs, paste(work.parms("path"), "ranks/", file.name("revised_pvals", x$sample.size[1], x$sim[1]), sep=""))
  })  
}

stan_probs_corr.rds1 = function(){
  d = readRDS(paste(work.parms("path"), "stan-probs/stan_probs_corr.rds", sep=""))
  ddply(d, .variables = .(sample.size, sim), .fun = function(x){
#    probs = 1 - pmax(x$phph, x$plph)
    probs = 1 - (x$phph + x$plph)   
    names(probs) = x$geneid
    saveRDS(probs, paste(work.parms("path"), "ranks/", file.name("stan_corr", x$sample.size[1], x$sim[1]), sep=""))
  })  
}

stan_probs.rds1 = function(){
  d = readRDS(paste(work.parms("path"), "stan-probs/stan_probs_corr.rds", sep=""))
  ddply(d, .variables = .(sample.size, sim), .fun = function(x){
#    probs = 1 - pmax(x$phph, x$plph)
    probs = 1 - (x$phph + x$plph)   
    names(probs) = x$geneid
    saveRDS(probs, paste(work.parms("path"), "ranks/", file.name("stan", x$sample.size[1], x$sim[1]), sep=""))
  })  
}

unpack.stan.sim1 = function(mtds){
  laplace_mvn1()                         # stan_laplace_mvn
#  post_probs_horseshoe.rds()  # stan_horseshoe
  post_probs_laplace.rds1()       # stan_laplace
  rev_probs_laplace.rds1()         # rev_laplace
  rev_probs_mu_cov.rds1()        # rev_probs_mu_cov
  revised_pvals.rds1()                # revised_pvals
  stan_probs_corr.rds1()            # stan_corr
  stan_probs.rds1()                    # stan
}


post_laplace2 = function(){
  file.begin = paste(work.parms("path"), "stan-probs/laplace_normal/post_laplace", sep="")
  for(size in work.parms("size")){
    for(rep in work.parms("rep")){
      x = readRDS(paste(file.begin, size, "-", rep, ".rds", sep=""))
#     probs = 1 - pmax(x$phph, x$plph)
      probs = 1 - (x$phph + x$plph)
      saveRDS(probs, paste(work.parms("path"), "ranks/", file.name("stan_laplace", size, rep), sep=""))
    }
  }
}

post_phimix2 = function(){
  file.begin = paste(work.parms("path"), "stan-probs/laplace_normal/post_phimix", sep="")
  for(size in work.parms("size")){
    for(rep in work.parms("rep")){
      x = readRDS(paste(file.begin, size, "-", rep, ".rds", sep=""))
#     probs = 1 - pmax(x$phph, x$plph)
      probs = 1 - (x$phph + x$plph)
      saveRDS(probs, paste(work.parms("path"), "ranks/", file.name("stan_laplace_phimix", size, rep), sep=""))
    }
  }
}

unpack.stan.sim2 = function(mtds){
  post_laplace2()                         # stan_laplace
  post_phimix2()                         # stan_laplace_phimix
}