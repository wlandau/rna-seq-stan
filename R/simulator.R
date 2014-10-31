simulate.datasets = function(which.datasets = 100){

  set.seed(10292014)

  initial = readRDS("../data/initial.rds") #empirical estimates from real data obtained with edgeR
  G = ncol(initial$B)

  truth = with(initial,apply(B,MARGIN=2,FUN = 
                     function(x) {as.numeric(x[3] > pmax(x[2],x[1])) - 
                                        as.numeric(x[3] < pmin(x[1],x[2]))}))
  id = names(truth)

  for (i in which.datasets){
    print(paste("Simulating dataset", i))

    print("  Calling ldply")
    df <-
      ldply(1:G,function(x){
        mu = exp(rep(initial$B[,x],each=4)+initial$c)
        y=rnbinom(12,mu=mu,size=exp(-initial$lpsi[x]))
        t(y)
      })

    print("  Trimming genes")  
    trimmed = trim_genes(counts=df,group=rep(1:3,each=4),geneid=id)

    print("  Cleaning up")
    counts.t = trimmed[[1]]
    geneid.t = trimmed[[2]]
    rownames(counts.t) = geneid.t
    keep = sort(sample(geneid.t,25000))

    saveRDS(counts.t[keep,], file=sprintf("../simulations/sim%i.rds",i))
    saveRDS(truth[keep], file=sprintf("../simulations/truth%i.rds",i))
  }
}