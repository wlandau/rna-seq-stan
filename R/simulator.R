simulator = function(rep = 1:10, size = 4, keep = work.parms("keep")){
  
  set.seed(10292014+size)
  seeds = sample.int(1e8, 50)
  group = as.factor(rep(c("parent1", "parent2", "hybrid"), each = size))
  saveRDS(group, "../data/group.rds")
  
  initial = readRDS("../data/initial.rds") #empirical estimates from real data obtained with edgeR
  G = ncol(initial$B)
  
  truth = with(initial,apply(B,MARGIN=2,FUN = 
                               function(x) {as.numeric(x[3] > pmax(x[2],x[1])) - 
                                              as.numeric(x[3] < pmin(x[1],x[2]))}))
  
  id = names(truth)
  
  Vectorize(function(size, rep){
    print(paste("Simulation", size, rep))
    set.seed(seeds[rep])
    print("  Calling ldply")
    df <-
      ldply(1:G,function(x){
        mu = exp(rep(initial$B[,x],each=size) + rep(initial$c, each = size / 4))
        y=rnbinom(3*size,mu=mu,size=exp(-initial$lpsi[x]))
        t(y)
      })
    
    print("  Trimming genes")  
    trimmed = trim_genes(counts=df,group=rep(1:3,each=size),geneid=id)
    
    print("  Cleaning up")
    counts.t = trimmed[[1]]
    geneid.t = trimmed[[2]]
    rownames(counts.t) = geneid.t
    keep = sort(sample(geneid.t, keep))
    
    counts.keep = counts.t[keep,]
    truth.keep = truth[keep]
    
    saveRDS(counts.keep, file=sprintf("../simulations/sim-%i-%i.rds", size, rep))
    saveRDS(truth.keep, file=sprintf("../simulations/truth-%i-%i.rds", size, rep))
  }, "rep")(size, rep)
}