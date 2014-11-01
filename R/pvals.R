edgeR.fit = function(counts, group){
  dge <- DGEList(counts = counts, group = group)
  design <- model.matrix(~ -1 + group)
  dge <- calcNormFactors(dge)
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)

  glmFit(dge, design)
}

edgeR.ord = function(fit){
  apply(fit$coef, 1, function(x){
    if(x[2] < min(x[c(1, 3)]))
      return("hybrid < min")
    else if(x[2] > max(x[c(1, 3)]))
      return("hybrid > max")
    else
      return("min <= hybrid <= max")
  })
}

pvals1dataset = function(pkg, counts, group, ncpus = 2){
 
  t = proc.time()

  if(pkg == "edgeR"){
    fit = edgeR.fit(counts, group)
    ord = edgeR.ord(fit)

    p1 = glmLRT(fit, contrast = c(-1, 1, 0))$table$PValue
    p3 = glmLRT(fit, contrast = c(0, 1, -1))$table$PValue

    ret = apply(cbind(ord, p1, p3), 1, function(x){
      if(x[1] == "min <= hybrid <= max")
        return(1)
      else if(x[1] == "hybrid < min")
        return(min(as.numeric(x[2:3]))/2)
      else if(x[1] == "hybrid > max")
        return(max(as.numeric(x[2:3]))/2)
    })
  } else if(pkg == "baySeq"){

    ord = edgeR.ord(edgeR.fit(counts, group))

    models = list(
      m111 = rep(1, length(group)),
      m113 = c(rep(1, sum(group != 3)), rep(3, sum(group == 3))),
      m122 = c(rep(1, sum(group == 1)), rep(2, sum(group != 1))),
      m121 = c(rep(1, sum(group == 1)), rep(2, sum(group == 2)), rep(1, sum(group == 3))),
      m123 = group
    )

    cl = makeCluster(ncpus, "SOCK")
    cd = new("countData", data = counts, replicates = group, groups= models)
    libsizes(cd) = getLibsizes(cd)
    cd = getPriors.NB(cd, cl = cl)
    cd = getLikelihoods.NB(cd, cl = cl)
    stopCluster(cl)

    post = apply(exp(cd@posteriors)[,c("m121", "m123")], 1, sum)
    post[ord == "min <= hybrid <= max"] = 0
    ret = 1 - post
  } else if(pkg == "ShrinkBayes"){

    phi = rep(1, length(group))
    alp = (group == 3) - (group == 1)
    del = as.integer(group == 2)

    size <- calcNormFactors(counts)
    libsize <- apply(counts, 2, sum)
    libsize <- libsize/exp(mean(log(libsize)))
    logsize <- log(size) + log(libsize)
    logsize <- logsize - mean(logsize)
    lgsz = offset(logsize)

    form = y ~ 0 + phi + alp + del + lgsz
    lc = inla.make.lincombs(alp = c(1,1), del = c(1, -1), phi = c(0, 0))

    priors = ShrinkSeq(form = form, dat = counts, shrinkfixed = "phi", shrinkaddfixed = list("alp", "del"), fams = "nb", shrinkdisp = T, addfixedmeanzero = F, ncpus = ncpus, ntag = ceiling(dim(counts)[1]/2)) # high ntag

    # finalprior = T to get rid of mixture prior
    # ncpus doesn't work in FitAllShrink
    fit = FitAllShrink(forms = form, dat = counts, shrinksimul = priors, fams="nb", finalprior = T, lincomb = lc) 

    postsLc1 = BFUpdatePosterior(fit, priors, shrinklc = "lc1", ncpus = ncpus) # LPH
    postsLc2 = BFUpdatePosterior(fit, priors, shrinklc = "lc2", ncpus = ncpus) # HPH

    post1 = SummaryWrap(postsLc1, thr = 0, direction = "greater", ncpus = ncpus) # P(alp + del > 0) = P(del > -alp)
    post2 = SummaryWrap(postsLc2, thr = 0, direction = "lesser", ncpus = ncpus) # P(alp - del < 0) = P(del > alp)  
    post3 = SummaryWrap(postsLc2, thr = 0, direction = "greater", ncpus = ncpus) # P(alp - del > 0) = P(del < alp)
    post4 = SummaryWrap(postsLc1, thr = 0, direction = "lesser", ncpus = ncpus) # P(alp + del < 0) = P(del < -alp)

    postsAlp = BFUpdatePosterior(fit, priors, shrinkpara = "alp", ncpus = ncpus)
    postmeansAlp = SummaryWrap(postsAlp, summary="postmean", ncpus = ncpus)

    postsDel = BFUpdatePosterior(fit, priors, shrinkpara = "alp", ncpus = ncpus)
    postmeansDel = SummaryWrap(postsDel, summary="postmean", ncpus = ncpus)

    df = cbind(postmeansAlp, postmeansDel, post1, post2, post3, post4)
    colnames(df) = c("alp", "del", paste("post", 1:4, sep=""))
 
    ret = apply(df, 1, function(x){
      if(abs(x["del"]) < abs(x["alp"]))
        return(0)
      else if(x["del"] > -x["alp"])
        return(x["post1"])
      else if(x["del"] > x["alp"])
        return(x["post2"])
      else if(x["del"] < x["alp"])
        return(x["post3"])
      else if(x["del"] < -x["alp"])
        return(x["post4"])    
    })
  }

  s = proc.time()
  logfile(paste(s - t, collapse = T))
  return(ret)
}

pvals = function(which.datasets = 1:100, ncpus = 2){
  group = readRDS("../data/group.rds")
  
  for(pkg in c("edgeR", "baySeq", "ShrinkBayes")){
    for(i in which.datasets){
      logfile(pkg, "pvals data", i)
      cts = readRDS(paste("../simulations/sim", i, ".rds", sep=""))
      saveRDS(pvals1dataset(pkg, cts, group, ncpus = ncpus), paste("../pvals/", pkg, i, ".rds", sep=""))
    }
    logfile("Done with", pkg, ".\n")
  }
}