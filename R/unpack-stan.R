laplace_mvn = function(){
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
post_probs_horseshoe.rds = function(){
  d = readRDS(paste(work.parms("path"), "stan-probs/post_probs_horseshoe.rds", sep=""))
  ddply(d, .variables = .(sample.size, sim), .fun = function(x){
#    probs = 1 - pmax(x$phph, x$plph)
    probs = 1 - (x$phph + x$plph)   
    names(probs) = x$geneid
    saveRDS(probs, paste(work.parms("path"), "ranks/", file.name("stan_horseshoe", x$sample.size[1], x$sim[1]), sep=""))
  })  
}

post_probs_laplace.rds = function(){
  d = readRDS(paste(work.parms("path"), "stan-probs/post_probs_laplace.rds", sep=""))
  ddply(d, .variables = .(sample.size, sim), .fun = function(x){
#    probs = 1 - pmax(x$phph, x$plph)
    probs = 1 - (x$phph + x$plph)   
    names(probs) = x$geneid
    saveRDS(probs, paste(work.parms("path"), "ranks/", file.name("stan_laplace", x$sample.size[1], x$sim[1]), sep=""))
  })  
}

rev_probs_laplace.rds = function(){
  d = readRDS(paste(work.parms("path"), "stan-probs/rev_probs_laplace.rds", sep=""))
  ddply(d, .variables = .(sample.size, sim), .fun = function(x){
#    probs = 1 - pmax(x$phph, x$plph)
#    probs = 1 - (x$phph + x$plph)   
    probs = 1 - x$new_prob   
    names(probs) = x$geneid
    saveRDS(probs, paste(work.parms("path"), "ranks/", file.name("rev_laplace", x$sample.size[1], x$sim[1]), sep=""))
  })  
}

rev_probs_mu_cov.rds = function(){
  d = readRDS(paste(work.parms("path"), "stan-probs/rev_probs_mu_cov.rds", sep=""))
  ddply(d, .variables = .(sample.size, sim), .fun = function(x){
#    probs = 1 - pmax(x$phph, x$plph)
#    probs = 1 - (x$phph + x$plph)   
    probs = 1 - x$new_prob
    names(probs) = x$geneid
    saveRDS(probs, paste(work.parms("path"), "ranks/", file.name("rev_probs_mu_cov", x$sample.size[1], x$sim[1]), sep=""))
  })  
}

revised_pvals.rds = function(){
  d = readRDS(paste(work.parms("path"), "stan-probs/revised_pvals.rds", sep=""))
  ddply(d, .variables = .(sample.size, sim), .fun = function(x){
#    probs = 1 - pmax(x$phph, x$plph)
#    probs = 1 - (x$phph + x$plph)   
    probs = 1 - x$V1
    names(probs) = x$geneid
    saveRDS(probs, paste(work.parms("path"), "ranks/", file.name("revised_pvals", x$sample.size[1], x$sim[1]), sep=""))
  })  
}

stan_probs_corr.rds = function(){
  d = readRDS(paste(work.parms("path"), "stan-probs/stan_probs_corr.rds", sep=""))
  ddply(d, .variables = .(sample.size, sim), .fun = function(x){
#    probs = 1 - pmax(x$phph, x$plph)
    probs = 1 - (x$phph + x$plph)   
    names(probs) = x$geneid
    saveRDS(probs, paste(work.parms("path"), "ranks/", file.name("stan_corr", x$sample.size[1], x$sim[1]), sep=""))
  })  
}

stan_probs.rds = function(){
  d = readRDS(paste(work.parms("path"), "stan-probs/stan_probs_corr.rds", sep=""))
  ddply(d, .variables = .(sample.size, sim), .fun = function(x){
#    probs = 1 - pmax(x$phph, x$plph)
    probs = 1 - (x$phph + x$plph)   
    names(probs) = x$geneid
    saveRDS(probs, paste(work.parms("path"), "ranks/", file.name("stan", x$sample.size[1], x$sim[1]), sep=""))
  })  
}

unpack.stan = function(mtds){
  laplace_mvn()                         # stan_laplace_mvn
#  post_probs_horseshoe.rds()  # stan_horseshoe
  post_probs_laplace.rds()       # stan_laplace
  rev_probs_laplace.rds()         # rev_laplace
  rev_probs_mu_cov.rds()        # rev_probs_mu_cov
  revised_pvals.rds()                # revised_pvals
  stan_probs_corr.rds()            # stan_corr
  stan_probs.rds()                    # stan
}