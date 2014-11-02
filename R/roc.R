pr = Vectorize(function(cutoff, pvals, n){
  sum(pvals <= cutoff)/n
}, "cutoff")

roc = function(pvals, truth, upper = 1e-1, npoints = 1e3){
  truth = truth[order(pvals)]
  pvals = sort(pvals)
  cutoffs = 1:npoints/npoints

  pvals.false = pvals[truth == 0]
  pvals.true = pvals[truth == 1]

  fpr = pr(cutoffs, pvals.false, length(pvals.false))
  tpr = pr(cutoffs, pvals.true, length(pvals.true))

  stepf <- stepfun(x = fpr, y = c(0, tpr))
  fpr = cutoffs*upper
  tpr = stepf(fpr)

  data.frame(fpr = fpr, tpr = tpr)
}

rocs = function(which.datasets = 1:100, upper = 1e-1){
  for(pkg in c("edgeR")){
#  for(pkg in c("edgeR", "baySeq", "ShrinkBayes")){
    for(i in which.datasets){
      pvals = readRDS(paste("../pvals/", pkg, i, ".rds", sep=""))
      truth = readRDS(paste("../simulations/truth", i, ".rds", sep=""))^2
      saveRDS(roc(pvals, truth, upper = upper), paste("../roc/", pkg, i, ".rds", sep=""))
    }
  }
}

auc = function(.roc){
  trapz(x = .roc$fpr, y = .roc$tpr)
}

aucs = function(which.datasets = 1:100){
  auc.list = list()
  for(pkg in c("edgeR")){ 
#  for(pkg in c("edgeR", "baySeq", "ShrinkBayes")){ 
   ret[[pkg]] = c()
    for(i in which.datasets){
      .roc = readRDS(paste("../roc/", pkg, i, ".rds", sep=""))
      auc.list[[pkg]] = c(auc.list[[pkg]], auc(.roc))
    }   
  }
  
  ret = data.frame(
    auc = c(auc.list[["edgeR"]]),
#    auc = c(auc.list[["edgeR"]], auc.list[["baySeq"]], auc.list[["ShrinkBayes"]]),
    pkg = rep(c("edgeR"), each = length(which.datasets))
#    pkg = rep(c("edgeR", "baySeq", "ShrinkBayes"), each = length(which.datasets))
  )
  
  saveRDS(ret, "../auc/auc.rds")
}