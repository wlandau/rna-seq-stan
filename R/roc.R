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

rocs = function(which.datasets = 1:100){
#  for(pkg in c("edgeR")){
  for(pkg in c("edgeR", "baySeq", "ShrinkBayes")){
    for(i in which.datasets){
      logfile(pkg, "ROC dataset", i)
      pvals = readRDS(paste("../pvals/", pkg, i, ".rds", sep=""))
      truth = readRDS(paste("../simulations/truth", i, ".rds", sep=""))^2
      saveRDS(roc(pvals, truth), paste("../roc/", pkg, i, ".rds", sep=""))
    }
  }
}

auc = function(.roc, upper = 1e-1){
  m = max(which(.roc$fpr <= upper))
  u = 2:m
  l = 1:(m - 1)
  sum(.roc$tpr[l] * (.roc$fpr[u] - .roc$fpr[l]))
}

aucs = function(which.datasets = 1:100){
  auc.list = list()
#  for(pkg in c("edgeR")){ 
  for(pkg in c("edgeR", "baySeq", "ShrinkBayes")){ 
   auc.list[[pkg]] = c()
    for(i in which.datasets){
      logfile(pkg, "AUC dataset", i)
      .roc = readRDS(paste("../roc/", pkg, i, ".rds", sep=""))
      auc.list[[pkg]] = c(auc.list[[pkg]], auc(.roc))
    }   
  }
  
  ret = data.frame(
#    auc = c(auc.list[["edgeR"]]),
    auc = c(auc.list[["edgeR"]], auc.list[["baySeq"]], auc.list[["ShrinkBayes"]]),
#    pkg = rep(c("edgeR"), each = length(which.datasets))
    pkg = rep(c("edgeR", "baySeq", "ShrinkBayes"), each = length(which.datasets))
  )
  
  saveRDS(ret, "../auc/auc.rds")
}