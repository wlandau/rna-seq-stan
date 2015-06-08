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

ranks1dataset = function(mtd, size, rep, counts, group, ncpus = 2){
 
# counts = counts[1:100,]

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

    nsize <- calcNormFactors(counts)
    libsize <- apply(counts, 2, sum)
    libsize <- libsize/exp(mean(log(libsize)))
    logsize <- log(nsize) + log(libsize)
    logsize <- logsize - mean(logsize)
    logsize <<- logsize

    form = y ~ 0 + phi + alp + del + offset(logsize)
    lc = inla.make.lincombs(alp = c(1,1), del = c(1, -1), phi = c(0, 0))

    priors = ShrinkSeq(form = form, dat = counts, shrinkfixed = "phi", shrinkaddfixed = list("alp", "del"), fams = "nb", shrinkdisp = T, addfixedmeanzero = F, ncpus = ncpus, ntag = ceiling(dim(counts)[1]/2)) # high ntag

   saveRDS(priors, paste(work.parms("path"), "priors/", "ShrinkBayes-", size, "-", rep, ".rds", sep=""))

    # finalprior = T to get rid of mixture prior
    # ncpus doesn't work in FitAllShrink
    fit = FitAllShrink(forms = form, dat = counts, shrinksimul = priors, fams="nb", finalprior = T, lincomb = lc) 

#    postsLc1 = BFUpdatePosterior(fit, priors, shrinklc = "lc1", ncpus = ncpus) # LPH
#    postsLc2 = BFUpdatePosterior(fit, priors, shrinklc = "lc2", ncpus = ncpus) # HPH

    postsLc1 = andyspost(fit, priors, shrinklc = "lc1", ncpus = ncpus) # LPH
    postsLc2 = andyspost(fit, priors, shrinklc = "lc2", ncpus = ncpus) # HPH

    post1 = SummaryWrap(postsLc1, thr = 0, direction = "greater", ncpus = ncpus) # P(alp + del > 0) = P(del > -alp)
    post2 = SummaryWrap(postsLc2, thr = 0, direction = "lesser", ncpus = ncpus) # P(alp - del < 0) = P(del > alp)  
    post3 = SummaryWrap(postsLc2, thr = 0, direction = "greater", ncpus = ncpus) # P(alp - del > 0) = P(del < alp)
    post4 = SummaryWrap(postsLc1, thr = 0, direction = "lesser", ncpus = ncpus) # P(alp + del < 0) = P(del < -alp)

#    postsAlp = BFUpdatePosterior(fit, priors, shrinkpara = "alp", ncpus = ncpus)
    postsAlp = andyspost(fit, priors, shrinkpara = "alp", ncpus = ncpus)

    postmeansAlp = SummaryWrap(postsAlp, summary="postmean", ncpus = ncpus)

#    postsDel = BFUpdatePosterior(fit, priors, shrinkpara = "del", ncpus = ncpus)
    postsDel = andyspost(fit, priors, shrinkpara = "del", ncpus = ncpus)

    postmeansDel = SummaryWrap(postsDel, summary="postmean", ncpus = ncpus)

    df = cbind(postmeansAlp, postmeansDel, post1, post2, post3, post4)
    colnames(df) = c("alp", "del", paste("post", 1:4, sep=""))
 
    ret = 1 - apply(df, 1, function(x){
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
  } else if(mtd == "ShrinkBayesMu") {

#####################################
### CELL MEANS PARAMETERIZATION ###
#####################################

print(group)
print(head(counts))
str(str(counts))

    mu.parent1 = as.integer(group == "parent1")
    mu.parent2 = as.integer(group == "parent2")
    mu.hybrid = as.integer(group == "hybrid")

str(mu.parent1)
str(mu.parent2)
str(mu.hybrid)

    nsize <- calcNormFactors(counts)
    libsize <- apply(counts, 2, sum)
    libsize <- libsize/exp(mean(log(libsize)))
    logsize <- log(nsize) + log(libsize)
    logsize <- logsize - mean(logsize)
    logsize <<- logsize

    form = y ~ 0 + mu.parent1 + mu.parent2 + mu.hybrid + offset(logsize)
    lc = inla.make.lincombs(mu.parent1 = c(1,0), mu.parent2 = c(0, 1), mu.hybrid = c(-1, -1))

print(form)
print(lc)

print("priors")

    priors = ShrinkSeq(form = form, dat = counts, shrinkfixed = "mu.hybrid", shrinkaddfixed = list("mu.parent1", "mu.parent2"), fams = "nb", shrinkdisp = T, addfixedmeanzero = F, ncpus = ncpus, ntag = ceiling(dim(counts)[1]/2)) # high ntag

   saveRDS(priors, paste(work.parms("path"),"priors/", "ShrinkBayesMu-", size, "-", rep, ".rds", sep=""))

print("fitallshrink")

    # finalprior = T to get rid of mixture prior
    # ncpus doesn't work in FitAllShrink
    fit = FitAllShrink(forms = form, dat = counts, shrinksimul = priors, fams="nb", finalprior = T, lincomb = lc) 

print("bfupdateposterior")

#    postsLc1 = BFUpdatePosterior(fit, priors, shrinklc = "lc1", ncpus = ncpus) # LPH
#    postsLc2 = BFUpdatePosterior(fit, priors, shrinklc = "lc2", ncpus = ncpus) # HPH

    postsLc1 = andyspost(fit, priors, shrinklc = "lc1", ncpus = ncpus) # LPH
    postsLc2 = andyspost(fit, priors, shrinklc = "lc2", ncpus = ncpus) # HPH

print("summarywrap")

    p1.g.h = SummaryWrap(postsLc1, thr = 0, direction = "greater", ncpus = ncpus) # P(mu.parent1 - mu.hybrid > 0) = P(mu.parent1 > mu.hybrid)
    h.g.p2 = SummaryWrap(postsLc2, thr = 0, direction = "lesser", ncpus = ncpus) # P(mu.parent2 - mu.hybrid < 0) = P(mu.hybrid > mu.parent2)  
    p2.g.h = SummaryWrap(postsLc2, thr = 0, direction = "greater", ncpus = ncpus) # P(mu.parent2 - mu.hybrid > 0) = P(mu.parent2 > mu.hybrid)
    h.g.p1 = SummaryWrap(postsLc1, thr = 0, direction = "lesser", ncpus = ncpus) # P(mu.parent1 - mu.hybrid  < 0) = P(mu.hybrid > mu.parent1)

print("bfupdateposterior")

#  posts.mu.parent1 = BFUpdatePosterior(fit, priors, shrinkpara = "mu.parent1", ncpus = ncpus)
  posts.mu.parent1 = andyspost(fit, priors, shrinkpara = "mu.parent1", ncpus = ncpus)
  postmeans.mu.parent1 = SummaryWrap(posts.mu.parent1, summary="postmean", ncpus = ncpus)

#  posts.mu.parent2 = BFUpdatePosterior(fit, priors, shrinkpara = "mu.parent2", ncpus = ncpus)
  posts.mu.parent2 = andyspost(fit, priors, shrinkpara = "mu.parent2", ncpus = ncpus)
  postmeans.mu.parent2 = SummaryWrap(posts.mu.parent2, summary="postmean", ncpus = ncpus)
 
#  posts.mu.hybrid = BFUpdatePosterior(fit, priors, shrinkpara = "mu.hybrid", ncpus = ncpus)
  posts.mu.hybrid = andyspost(fit, priors, shrinkpara = "mu.hybrid", ncpus = ncpus)
  postmeans.mu.hybrid = SummaryWrap(posts.mu.hybrid, summary="postmean", ncpus = ncpus)

df = cbind(postmeans.mu.parent1, postmeans.mu.parent2, postmeans.mu.hybrid, p1.g.h, p2.g.h, h.g.p1, h.g.p2)
colnames(df) = c("p1", "p2", "h", "p1.g.h", "p2.g.h", "h.g.p1", "h.g.p2")
 
    ret = 1 - apply(df, 1, function(x){
      if((min(x["p1"], x["p2"]) <= x["h"]) && (x["h"] <= max(x["p1"], x["p2"])))
        return(0)
      else if(x["h"] > x["p1"] && (x["p1"] >= x["p2"]))
        return(x["h.g.p1"])
      else if(x["h"] > x["p2"] && (x["p2"] >= x["p1"]))
        return(x["h.g.p2"])
      else if(x["h"] < x["p1"] && (x["p1"] <= x["p2"]))
        return(x["p1.g.h"])
      else if(x["h"] < x["p2"] && (x["p2"] <= x["p1"]))
        return(x["p2.g.h"])
    })

  } else if(mtd == "DESeq2"){

    colData = DataFrame(Treatment = factor(rep(c("parent1", "parent2", "hybrid"), each = size)))

    se = SummarizedExperiment(assays = SimpleList(counts = counts), colData = colData)

    dds <- DESeqDataSet(se = se, design = ~ Treatment - 1)
    dds <- DESeq(dds, betaPrior = F)

    cf = coef(dds)

    p1 = cf[, "Treatmentparent1"]
    p2 = cf[, "Treatmentparent2"]
    h = cf[, "Treatmenthybrid"]

    p1.g.h = results(dds, altHypothesis = "greater", 
      contrast = c(Treatmenthybrid = -1, Treatmentparent1 = 1, Treatmentparent2 = 0))$padj

    h.g.p2 = results(dds, altHypothesis = "less", 
      contrast = c(Treatmenthybrid = -1, Treatmentparent1 = 0, Treatmentparent2 = 1))$padj

    p2.g.h = results(dds, altHypothesis = "greater", 
      contrast = c(Treatmenthybrid = -1, Treatmentparent1 = 0, Treatmentparent2 = 1))$padj

    h.g.p1 = results(dds, altHypothesis = "less", 
      contrast = c(Treatmenthybrid = -1, Treatmentparent1 = 1, Treatmentparent2 = 0))$padj

    df = cbind(p1, p2, h, p1.g.h, p2.g.h, h.g.p1, h.g.p2)

    ret = apply(df, 1, function(x){
      if((min(x["p1"], x["p2"]) <= x["h"]) && (x["h"] <= max(x["p1"], x["p2"])))
        return(0)
      else if(x["h"] > x["p1"] && (x["p1"] >= x["p2"]))
        return(x["h.g.p1"])
      else if(x["h"] > x["p2"] && (x["p2"] >= x["p1"]))
        return(x["h.g.p2"])
      else if(x["h"] < x["p1"] && (x["p1"] <= x["p2"]))
        return(x["p1.g.h"])
      else if(x["h"] < x["p2"] && (x["p2"] <= x["p1"]))
        return(x["p2.g.h"])
      else
        return(1)
    })

  } else if(mtd == "fullyBayes"){
    g1 = group
    group = (g1 == "parent1") + 2*(g1 == "hybrid") + 3*(g1 == "parent2")
    chain = Chain(counts, group, Configs(iterations = 1e4, burnin = 1e5, returns = character(0)))
    o = heterosis(chain)
    ret = 1 - (o@hph + o@lph)
  }

  s = proc.time()
  logfile(paste(s - t, collapse = " "))
  names(ret) = rownames(counts)
  return(ret)
}

ranks = function(mtds = c("edgeR", "baySeq", "ShrinkBayes"), sizes = c(4, 8, 16), reps = 1:10, ncpus = 2){
  
  for(mtd in mtds){
    for(size in sizes){
      group = as.factor(rep(c("parent1", "parent2", "hybrid"), each = size))

      for(rep in reps){
        logfile(mtd, "ranks", size, rep)

        print(paste(size, rep))

        cts = readRDS(paste(work.parms("path"), "datasets/", file.name("sim", size, rep), sep=""))
        saveRDS(ranks1dataset(mtd, size, rep, cts, group, ncpus = ncpus), paste(work.parms("path"),"ranks/", file.name(mtd, size, rep), sep=""))
      }
    }
    logfile("Done with", mtd, ".\n")
  }
}