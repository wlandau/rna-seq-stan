pr = Vectorize(function(cutoff, pvals, n){
  sum(pvals <= cutoff)/n
}, "cutoff")

roc = function(pv, truth){
  truth = truth[order(pv)]
  pv = sort(pv)

  pvals.false = pv[truth == 0]
  pvals.true = pv[truth != 0]

  fpr = pr(pv, pvals.false, length(pvals.false))
  tpr = pr(pv, pvals.true, length(pvals.true))

  data.frame(fpr = fpr, tpr = tpr)
}

rocs = function(which.datasets = 1:10, which.sizes = c(4, 8, 16), pkgs = c("edgeR", "baySeq", "ShrinkBayes")){
  for(pkg in pkgs){
    for(s in which.sizes){
      for(i in which.datasets){
        print(paste(s, i))      
        logfile(pkg, "ROC dataset", i, s)
        pvals = readRDS(paste("../pvals/", pkg, i, "-", s, ".rds", sep=""))
        truth = readRDS(paste("../simulations/truth", i, "-", s, ".rds", sep=""))^2
        saveRDS(roc(pvals, truth), paste("../roc/", pkg, i, "-", s, ".rds", sep=""))
      }
    }
  }
}

auc = function(.roc, upper = 1e-1){
  m = max(which(.roc$fpr <= upper))
  u = 2:m
  l = 1:(m - 1)
  sum(.roc$tpr[l] * (.roc$fpr[u] - .roc$fpr[l]))
}

aucs = function(which.datasets = 1:10, which.sizes = c(4, 8, 16), pkgs = c("edgeR", "baySeq", "ShrinkBayes")){
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