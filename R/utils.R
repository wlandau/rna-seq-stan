work.parms = function(item){
  mtd = c(
    "edgeR", 
    "baySeq", 
    "ShrinkBayes",
    "ShrinkBayesMu", 

    "stan_laplace_mvn",     # ../stan-probs/laplace_mvn/pprob_lap_mvn...
#    "stan_horseshoe",        # ../stan-probs/post_probs_horseshoe.rds, only has sample.size = 16
    "stan_laplace",             # ../stan-probs/post_probs_laplace.rds 
    "rev_laplace",              # ../stan-probs/rev_probs_laplace.rds
    "rev_probs_mu_cov",  # ../stan-probs/rev_probs_mu_cov.rds
    "revised_pvals",           # ../stan-probs/revised_pvals.rds
    "stan_corr",                  # ../stan-probs/stan_probs_corr.rds
    "stan"                           # ../stan-probs/stan_probs_corr.rds
  )

  simno = 2

  switch(item, 
             "sim" = simno,
             "path" = paste("../simulations/sim", simno, "/", sep=""),
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