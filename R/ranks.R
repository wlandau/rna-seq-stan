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
    if(x["grouphybrid"] < min(x[c("groupparent1", "groupparent2")]))
      return("hybrid < min")
    else if(x["grouphybrid"] > max(x[c("groupparent1", "groupparent2")]))
      return("hybrid > max")
    else
      return("min <= hybrid <= max")
  })
}

ranks1dataset = function(mtd, counts, group, ncpus = 2){
 
  t = proc.time()

  if(mtd == "edgeR"){
    fit = edgeR.fit(counts, group)

    p1 = glmLRT(fit, contrast = (unique(group) == "hybrid") - (unique(group) == "parent1"))$table$PValue
    p2 = glmLRT(fit, contrast = (unique(group) == "hybrid") - (unique(group) == "parent2"))$table$PValue

    ret = apply(cbind(fit$coef, p1, p2), 1, function(x){
      if(x["grouphybrid"] < x["groupparent1"] && x["groupparent1"] <= x["groupparent2"]){
        return(x["p1"]/2)
      } else if(x["grouphybrid"] < x["groupparent2"] && x["groupparent2"] <= x["groupparent1"]){
        return(x["p2"]/2)
      } else if(x["grouphybrid"] > x["groupparent1"] && x["groupparent1"] >= x["groupparent2"]){
        return(x["p1"]/2)
      } else if(x["grouphybrid"] > x["groupparent2"] && x["groupparent2"] >= x["groupparent1"]){
        return(x["p2"]/2)
      } else {
        return(1)
      }
    })
  } else if(mtd == "baySeq"){

    ord = edgeR.ord(edgeR.fit(counts, group))

    p1.diff = as.vector(group)
    p1.diff[p1.diff == "hybrid"] = p1.diff[p1.diff == "parent2"] = "others"
    p1.diff = as.factor(p1.diff)

    p2.diff = as.vector(group)
    p2.diff[p2.diff == "hybrid"] = p2.diff[p2.diff == "parent1"] = "others"
    p2.diff = as.factor(p2.diff)

    hy.diff = as.vector(group)
    hy.diff[hy.diff == "parent1"] = hy.diff[hy.diff == "parent2"] = "parents"
    hy.diff = as.factor(hy.diff)

    models = list(
      all.same = rep("same", length(group)),
      p1.diff = p1.diff,
      p2.diff = p2.diff,
      hy.diff = hy.diff,
      all.diff = group
    )

    cl = makeCluster(ncpus, "SOCK")
    cd = new("countData", data = counts, replicates = group, groups= models)
    libsizes(cd) = getLibsizes(cd)
    cd = getPriors.NB(cd, cl = cl)
    cd = getLikelihoods.NB(cd, cl = cl)
    stopCluster(cl)

    post = apply(exp(cd@posteriors)[, c("hy.diff", "all.diff")], 1, sum)
    post[ord == "min <= hybrid <= max"] = 0
    ret = 1 - post
  } else if(mtd == "ShrinkBayes"){

    phi = rep(1, length(group))
    alp = (group == "parent2") - (group == "parent1")
    del = as.integer(group == "hybrid")

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

    postsDel = BFUpdatePosterior(fit, priors, shrinkpara = "del", ncpus = ncpus)
    postmeansDel = SummaryWrap(postsDel, summary="postmean", ncpus = ncpus)

    df = cbind(postmeansAlp, postmeansDel, post1, post2, post3, post4)
    colnames(df) = c("alp", "del", paste("post", 1:4, sep=""))
 
    ret = apply(df, 1, function(x){
      if(abs(x["del"]) <= abs(x["alp"]))
        return(0)
#    else if(x["del"] > -x["alp"]) # before 12/2/14
      else if(x["del"] > -x["alp"] && -x["alp"] >= 0) # 12/2/14: changed in preparation for next round of simulations
        return(x["post1"])
#   else if(x["del"] > x["alp"] && x["alp"] >= 0) # before 12/2/14
     else if(x["del"] > x["alp"] && x["alp"] >= 0) # 12/2/14: changed in preparation for next round of simulations
        return(x["post2"])
#    else if(x["del"] < x["alp"]) # before 12/2/14
      else if(x["del"] < x["alp"] && x["alp"] <= 0) # 12/2/14: changed in preparation for next round of simulations
        return(x["post3"])
#    else if(x["del"] < -x["alp"]) # before 12/2/14
      else if(x["del"] < -x["alp"] && -x["alp"] <= 0) # 12/2/14: changed in preparation for next round of simulations
        return(x["post4"])    
    })
  }

  s = proc.time()
  logfile(paste(s - t, collapse = " "))
  return(ret)
}

ranks = function(mtds = c("edgeR", "baySeq", "ShrinkBayes"), sizes = c(4, 8, 16), reps = 1:10, ncpus = 2){
  
  for(mtd in mtds){
    for(size in sizes){
      group = as.factor(rep(c("parent1", "parent2", "hybrid"), each = size))

      for(rep in reps){
        logfile(mtd, "ranks", size, rep)
        cts = readRDS(paste("../simulations/", file.name("sim", size, rep), sep=""))
        saveRDS(ranks1dataset(mtd, cts, group, ncpus = ncpus), paste("../ranks/", file.name(mtd, size, rep), sep=""))
      }
    }
    logfile("Done with", mtd, ".\n")
  }
}