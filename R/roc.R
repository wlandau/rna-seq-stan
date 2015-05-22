fpr.point = Vectorize(function(cutoff, ranks, n){
  sum(ranks <= cutoff)/n
}, "cutoff")

roc = function(ranks, truth){
  ranks[is.na(ranks)] = 1

  truth = truth[order(ranks)]
  ranks = sort(ranks)


  ranks.false = ranks[truth == 0]
  ranks.true = ranks[truth != 0]

  fpr = fpr.point(ranks, ranks.false, length(ranks.false))
  tpr = fpr.point(ranks, ranks.true, length(ranks.true))

  data.frame(fpr = fpr, tpr = tpr)
}

rocs = function(mtd = work.parms("mtd"), size = work.parms("size"), rep = work.parms("rep")){
  irrelevance = loopify(function(mtd, size, rep){
    print(paste("ROC", mtd, size, rep))
    name = file.name(mtd, size, rep)
    ranks = readRDS(paste(work.parms("path"), "ranks/", name, sep=""))
    truth = readRDS(paste(work.parms("path"), "truth/", file.name("truth", size, rep), sep=""))^2
    saveRDS(roc(ranks, truth), paste(work.parms("path"), "roc/", name, sep=""))
  }, mtd, size, rep)
}

auc = function(.roc, upper = 1e-1){
  m = max(which(.roc$fpr <= upper))
  u = 2:m
  l = 1:(m - 1)
  sum(.roc$tpr[l] * (.roc$fpr[u] - .roc$fpr[l]))
}

aucs = function(mtd = work.parms("mtd"), size = work.parms("size"), rep = work.parms("rep")){
  ret = loopify(function(mtd, size, rep){
    print(paste("AUC", mtd, size, rep))
    .auc = auc(readRDS(paste(work.parms("path"), "roc/", file.name(mtd, size, rep), sep="")))
    c(.auc, mtd, size, rep)
  }, mtd, size, rep)

  colnames(ret) = c("auc", "mtd", "size", "rep")
  ret$auc = as.numeric(as.vector(ret$auc))
  ret$mtd = ordered(ret$mtd, labels = mtd)
  ret$size.short = ordered(ret$size, levels = size)
  ret$size = ordered(ret$size.short, labels = paste(size, "samples / group"))

  smry = ddply(ret, .variables = .(mtd, size), .fun = function(x){
    c(mean(x$auc), mean(x$auc) + c(-1, 1) * 1.96  * sd(x$auc)/sqrt(length(x$auc)))
  })
  smry = smry[rep(1:dim(smry)[1], each = 10),]

  ret = cbind(ret, mean = smry$V1, lower = smry$V2, upper = smry$V3)
  saveRDS(ret, paste(work.parms("path"), "auc/auc.rds", sep=""))
}