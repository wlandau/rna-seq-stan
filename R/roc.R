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

  as.data.frame(fdr = fdr, tpr = tpr)
}

auc = function(.roc){
  trapz(x = .roc$fpr, y = .roc$tpr)
}