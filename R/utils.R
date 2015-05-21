work.parms = function(item){

  simno = 1

  mtd = c(
    "edgeR", 
    "baySeq", 
    "ShrinkBayes",
    "ShrinkBayesMu", 
    "DESeq2",
    "fullyBayes"
  )

  if(simno == 1){
    mtd = c(mtd,
      "stan_laplace_mvn",     # ../simulations/sim1/stan-probs/laplace_mvn/pprob_lap_mvn...
 #     "stan_horseshoe",        # ../simulations/sim1/stan-probs/post_probs_horseshoe.rds, only has sample.size = 16
      "stan_laplace",             # ../simulations/sim1/stan-probs/post_probs_laplace.rds 
      "rev_laplace",               # ../simulations/sim1/stan-probs/rev_probs_laplace.rds
      "rev_probs_mu_cov",   # ../simulations/sim1/stan-probs/rev_probs_mu_cov.rds
      "revised_pvals",           # ../simulations/sim1/stan-probs/revised_pvals.rds
      "stan_corr",                  # ../simulations/sim1/stan-probs/stan_probs_corr.rds
      "stan"                           # ../simulations/sim1/stan-probs/stan_probs_corr.rds
    )
  } else if (simno == 2){
    mtd = c(mtd,
      "stan_laplace",                 # ../simulations/sim2/stan-probs/laplace_normal/posts_laplace...
      "stan_laplace_phimix"      # ../simulations/sim2/stan-probs/laplace_normal/posts_phimix...
    )
  }

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