trim.whitespace = function(x){
  gsub("^\\s+|\\s+$", "", x)
}

figdir = function(mtd = work.parms("mtd")){
  if(length(mtd) == length(work.parms("mtd")) && all(sort(mtd) == sort(work.parms("mtd")))){
    dr = paste(work.parms("path"), "fig-all/", sep="")
  } else {
    m = paste(paste(mtd, collapse="-"))
    if(length(mtd))
      m = paste("-", m, sep="")
    dr = paste(work.parms("path"), "fig", m, "/", sep="")
  }
  if(!file.exists(dr))
    dir.create(dr)
  return(dr)
}

relevel = function(d){
  if(!is.null(d$mtd)){
    lvl = levels(d$mtd)
    lvl[lvl == "stan"] = "RStan (ind)"
    lvl[lvl == "stan_corr"] = "RStan (cov)"
    lvl[lvl == "stan_laplace"] = "RStan (Laplace)"
    lvl[lvl == "stan_laplace_mvn"] = "RStan (Laplace MVN)"
    levels(d$mtd) = lvl   
  }
  return(d)
}

exampleROCdf = function(upper = 1e-1, mtd = work.parms("mtd")){
  l = loopify(function(mtd, size, rep){
    .roc = readRDS(paste(work.parms("path"), "roc/", file.name(mtd, size, 1), sep=""))
    .roc = .roc[.roc$fpr <= upper,]
    .roc$mtd = mtd
    .roc$size = size
    .roc
  }, mtd = mtd, rep = 1)

  d = data.frame(
    fpr = do.call("c", l$fpr),
    tpr = do.call("c", l$tpr),
    mtd = ordered(do.call("c", lapply(l$mtd, as.character)), levels = rev(work.parms("mtd"))),
    size = ordered(do.call("c", l$size), labels = paste(work.parms("size"), "samples / group"))
  )

  return(relevel(d))
}

plotExampleROC = function(upper = c(1e-1, 1), mtd = work.parms("mtd")){ 
  Vectorize(function(upper){
    pl = ggplot(exampleROCdf(upper, mtd), aes(x = fpr, y = tpr)) +
           facet_grid(~ size) + 
           geom_abline() +  
           geom_line(aes(color = mtd)) + 
           labs(color = "Method") +
           xlab("\nFalse positive rate (FPR)") + 
           ylab("True positive rate (TPR)\n") + 
           theme(axis.text = element_text(family = "Helvetica", color = 'black'),
                      axis.text.x = element_text(angle = -80, hjust = 0),
                      panel.background = element_rect(fill='white'),
                      panel.border = element_rect(color="black", fill = NA),
                      panel.grid.major = element_line(color="lightgray"),
                      text = element_text(family = "Helvetica", colour= "black"))

    dr = figdir(mtd)

    ggsave(paste(dr, "exampleROC", gsub("\\.", "_", as.character(upper)),
                 ".pdf", sep=""), pl, width = 8, height = 5, dpi = 1600)
  }, "upper")(upper)
}

plotAUCfacet = function(facet.by.size = c(T, F), mtd = work.parms("mtd")){
  Vectorize(function(facet.by.size){
    df = readRDS(paste(work.parms("path"), "auc/auc.rds", sep=""))
    df = df[df$mtd %in% mtd, ]
    df$mtd = ordered(df$mtd, work.parms("mtd"))

   dr = figdir(mtd)
   df = relevel(df)

    if(facet.by.size){
      pl = ggplot(df, aes(x = mtd, y = auc)) +
             facet_grid(~ size, scales = "free_x") + 
             theme(axis.text.x = element_text(angle = -80, hjust = 0)) + 
             xlab("\nMethod")
    } else {
      pl = ggplot(df, aes(x = size.short, y = auc)) +
             facet_grid(~ mtd, scales = "free_x") +
             xlab("\nSample Size per Treatment Group")
    }
     pl = pl +
           geom_line(aes(group = rep), alpha = .5) +  
           geom_crossbar(aes(ymin = lower, y = mean, ymax = upper), fatten = 1, width = .5, color = "blue") + 
           ylab("Area under Receiver Operating Characteristic\n(ROC) curve (false positive rate < 0.1)\n") + 
           theme(axis.text = element_text(family = "Helvetica", colour = 'black'),
                      legend.position="none",
                      panel.background = element_rect(fill='white'),
                      panel.border = element_rect(color="black", fill = NA),
                      panel.grid.major = element_line(color="lightgray"),
                      text = element_text(family = "Helvetica", colour= "black"))
    
    ggsave(paste(dr, "auc-facet-", facet.by.size, ".pdf", sep=""), pl, width = 8, height = 5, dpi = 1600)
  }, "facet.by.size")(facet.by.size)
}

plotAUCcolor = function(jitter = c(T, F), mtd = work.parms("mtd")){
  Vectorize(function(jitter){
    df = readRDS(paste(work.parms("path"), "auc/auc.rds", sep=""))
    df = df[df$mtd %in% mtd, ]
    df$mtd = ordered(df$mtd, rev(mtd))

   dr = figdir(mtd)
   df = relevel(df)

    pl = ggplot(df, aes(x = size.short, y = auc, color = mtd), na.rm=TRUE) +
           geom_crossbar(aes(ymin = lower, y = mean, ymax = upper), fatten = 1, width = .25) + 
           labs(color = "Method") +
           xlab("\nSample Size per Treatment Group") +
           ylab("Area under Receiver Operating Characteristic\n(ROC) curve (false positive rate < 0.1)\n") + 
           theme(axis.text = element_text(family = "Helvetica", colour = 'black'),
                      panel.background = element_rect(fill='white'),
                      panel.border = element_rect(color="black", fill = NA),
                      panel.grid.major = element_line(color="lightgray"),
                      text = element_text(family = "Helvetica", colour= "black"))
    
    if(jitter)
      pl = pl + geom_jitter()
    else
      pl = pl + geom_point()    
    
    ggsave(paste(dr, "auc-color-", jitter, ".pdf", sep=""), pl, width = 8, height = 5, dpi = 1600)
  }, "jitter")(jitter)
}

plotFDRfacet = function(facet.direction = c(T, F), who = c("dan", "jarad"), y.axis = c("subtract", "leave"),
reverse.x = c(" ", " no "), mtd = work.parms("mtd")){
  grid = expand.grid(who, facet.direction, y.axis, reverse.x)
  
  Vectorize(function(who, facet.direction, y.axis, reverse.x){
    d = readRDS(paste(work.parms("path"), "fdr/",who, ".rds", sep=""))
    d = d[d$mtd %in% mtd, ]

   dr = figdir(mtd)
   d = relevel(d)
  
    if(who == "dan") reverse.x = " no "

    x.lab = ifelse(who == "dan", "est. Bayesian FDR", paste("mean post. probability of", reverse.x, "heterosis", sep=""))
    y.lab = ifelse(y.axis == "subtract", paste("FDP", "-", x.lab), "FDP")

    if(y.axis == "leave" && reverse.x == " no ")
      pl = ggplot(d, aes(x = cutoff, y = fdr)) 
    else if(y.axis == "subtract" && reverse.x == " no ")
      pl = ggplot(d, aes(x = cutoff, y = fdr.minus.cutoff))  
    else if(y.axis == "leave" && reverse.x == " ")
      pl = ggplot(d, aes(x = one.minus.cutoff, y = fdr)) 
    else if(y.axis == "subtract" && reverse.x == " ")
      pl = ggplot(d, aes(x = one.minus.cutoff, y = fdr.minus.cutoff))

    pl = pl + geom_abline(slope = as.integer(y.axis == "leave"), 
                                 intercept = as.integer(reverse.x == " ") * as.integer(y.axis == "leave"), 
                                 color = "blue") + 
           geom_line(aes(group = rep), alpha = 0.5) + 
           xlab(paste("\n", x.lab, sep="")) + 
           ylab(paste(y.lab, "\n", sep="")) + 
           theme(axis.text.x = element_text(angle = -80, hjust = 0))

    if(facet.direction) pl = pl + facet_grid(mtd ~ size) else pl = pl + facet_grid(size ~ mtd)

     if(reverse.x == " ")
    pl = pl + scale_x_reverse()

     ggsave(paste(dr, "fdr-facet-",who, "-", facet.direction, "-", y.axis, "-", trim.whitespace(reverse.x), ".pdf", sep=""), pl, width = 8, height = 8, dpi = 1600)
  })(grid[[1]], grid[[2]], grid[[3]], grid[[4]])
}

plotFDRindiv = function(who = c("dan", "jarad"), mtd = work.parms("mtd"), y.axis = c("subtract", "leave"), 
reverse.x = c(" ", " no ")){

  grid = expand.grid(who, mtd, y.axis, reverse.x)
   dr = figdir(mtd)  

  Vectorize(function(who, mtd, y.axis, reverse.x){
    d = readRDS(paste(work.parms("path"), "fdr/",who, ".rds", sep=""))
    d = d[d$mtd == mtd,]
    d = relevel(d)

    if(who == "dan") reverse.x = " no "

    x.lab = ifelse(who == "dan", "est. Bayesian FDR", paste("mean post. probability of", reverse.x, "heterosis", sep=""))
    y.lab = ifelse(y.axis == "subtract", paste("FDP", "-", x.lab), "FDP")

    if(y.axis == "leave" && reverse.x == " no ")
      pl = ggplot(d, aes(x = cutoff, y = fdr)) 
    else if(y.axis == "subtract" && reverse.x == " no ")
      pl = ggplot(d, aes(x = cutoff, y = fdr.minus.cutoff))  
    else if(y.axis == "leave" && reverse.x == " ")
      pl = ggplot(d, aes(x = one.minus.cutoff, y = fdr)) 
    else if(y.axis == "subtract" && reverse.x == " ")
      pl = ggplot(d, aes(x = one.minus.cutoff, y = fdr.minus.cutoff))

    pl = pl + facet_grid(~ size) +
           geom_abline(slope = as.integer(y.axis == "leave"), 
                                 intercept = as.integer(reverse.x == " ") * as.integer(y.axis == "leave"), 
                                 color = "blue") + 
           geom_line(aes(group = rep), alpha = 0.5) + 
           xlab(paste("\n", x.lab, sep="")) + 
           ylab(paste(y.lab, "\n", sep="")) + 
           theme(axis.text.x = element_text(angle = -80, hjust = 0))

    if(reverse.x == " ")
      pl = pl + scale_x_reverse()

     ggsave(paste(dr, "fdr-indiv-",who, "-", mtd, "-", y.axis, "-", trim.whitespace(reverse.x), ".pdf", sep=""), pl, width = 8, height = 5, dpi = 1600)
  })(grid[[1]], grid[[2]], grid[[3]], grid[[4]])
}

makePlots = function(mtd = work.parms("mtd")){
  plotExampleROC(mtd = mtd)
#  plotAUCcolor(mtd = mtd)
  plotAUCfacet(mtd = mtd)
  plotFDRfacet(mtd = mtd)
#  plotFDRindiv(mtd = mtd)
}