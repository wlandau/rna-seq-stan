pr = Vectorize(function(cutoff, ranks, n){
  sum(pvals <= cutoff)/n
}, "cutoff")

roc = function(ranks, truth){
  truth = truth[order(ranks)]
  ranks = sort(ranks)

  ranks.false = ranks[truth == 0]
  ranks.true = ranks[truth != 0]

  fpr = pr(ranks, ranks.false, length(ranks.false))
  tpr = pr(ranks, ranks.true, length(ranks.true))

  data.frame(fpr = fpr, tpr = tpr)
}

rocs = Vectorize(function(mtd, size, rep){
  print(paste("ROC", mtd, s, i))
  ranks = readRDS(paste("../pvals/", pkg, i, "-", s, ".rds", sep=""))
  truth = readRDS(paste("../simulations/truth", i, "-", s, ".rds", sep=""))^2
  saveRDS(roc(ranks, truth), paste("../roc/", pkg, i, "-", s, ".rds", sep=""))
}, c("mtd", "size", "rep"))

auc = function(.roc, upper = 1e-1){
  m = max(which(.roc$fpr <= upper))
  u = 2:m
  l = 1:(m - 1)
  sum(.roc$tpr[l] * (.roc$fpr[u] - .roc$fpr[l]))
}

aucs = function(which.datasets = 1:10, which.sizes = c(4, 8, 16), pkgs = c("edgeR", "baySeq", "ShrinkBayes", "stan_corr", "stan", "stan_horseshoe")){
  ret = NULL
  for(pkg in pkgs){
    for(s in which.sizes){
      for(i in which.datasets){
        print(paste(s, i))
        logfile(pkg, "AUC dataset", i, s)
        .roc = readRDS(paste("../roc/", pkg, i, "-", s, ".rds", sep=""))
        .auc = auc(.roc)
        ret = rbind(ret, c(.auc, pkg, i, s))
      }   
    }
  }

  ret = data.frame(ret)
  colnames(ret) = c("AUC", "Method", "Replicate", "SampleSize")
  ret$AUC = as.numeric(as.vector(ret$AUC))
  
  saveRDS(ret, "../auc/auc.rds")
}