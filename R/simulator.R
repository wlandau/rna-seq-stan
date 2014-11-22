simulator = function(which.datasets = 1:10, sample.size = 4, n.keep = work.parms("keep")){
  
  set.seed(10292014+sample.size)
  seeds = sample.int(1e8, 50)
  group = as.factor(rep(c("parent1", "parent2", "hybrid"), each = sample.size))
  saveRDS(group, "../data/group.rds")
  
  initial = readRDS("../data/initial.rds") #empirical estimates from real data obtained with edgeR
  G = ncol(initial$B)
  
  truth = with(initial,apply(B,MARGIN=2,FUN = 
                               function(x) {as.numeric(x[3] > pmax(x[2],x[1])) - 
                                              as.numeric(x[3] < pmin(x[1],x[2]))}))
  
  # UNCOMMENT BELOW TO FORCE MORE EXTREME HETEROSIS
  # initial$B = apply(initial$B, 2, function(x){
  #    if(as.numeric(x[3] > pmax(x[2],x[1])))
  #       x[3] = x[3] + 1
  #    else if(as.numeric(x[3] < pmin(x[1],x[2])))
  #       x[3] = x[3] - 1 
  #   return(x)
  #  })
  
  id = names(truth)
  
  for (i in which.datasets){
    print(paste("Simulating dataset", i))
    set.seed(seeds[i])
    print("  Calling ldply")
    df <-
      ldply(1:G,function(x){
        mu = exp(rep(initial$B[,x],each=sample.size) + rep(initial$c, each = sample.size / 4))
        y=rnbinom(3*sample.size,mu=mu,size=exp(-initial$lpsi[x]))
        t(y)
      })
    
    print("  Trimming genes")  
    trimmed = trim_genes(counts=df,group=rep(1:3,each=sample.size),geneid=id)
    
    print("  Cleaning up")
    counts.t = trimmed[[1]]
    geneid.t = trimmed[[2]]
    rownames(counts.t) = geneid.t
    keep = sort(sample(geneid.t, n.keep))
    
    counts.keep = counts.t[keep,]
    truth.keep = truth[keep]
    
    saveRDS(counts.keep, file=sprintf("../simulations/sim%i-%i.rds",i,sample.size))
    saveRDS(truth.keep, file=sprintf("../simulations/truth%i-%i.rds",i,sample.size))
  }
}