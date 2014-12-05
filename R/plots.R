exampleROCdf = function(upper = 1e-1){
  l = loopify(function(mtd, size, rep){
    .roc = readRDS(paste("../roc/", file.name(mtd, size, 1), sep=""))
    .roc = .roc[.roc$fpr <= upper,]
    .roc$mtd = mtd
    .roc$size = size
    .roc
  }, rep = 1)

  data.frame(
    fpr = do.call("c", l$fpr),
    tpr = do.call("c", l$tpr),
    mtd = ordered(do.call("c", lapply(l$mtd, as.character)), levels = rev(work.parms("mtd"))),
    size = ordered(do.call("c", l$size), labels = paste(work.parms("size"), "samples / group"))
  )
}

plotExampleROC = function(upper = c(1e-1, 1)){ 
  Vectorize(function(upper){
    pl = ggplot(exampleROCdf(upper), aes(x = fpr, y = tpr)) +
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

    ggsave(paste("../fig/exampleROC", upper, ".pdf", sep=""), pl, width = 8, height = 5, dpi = 1600)
  }, "upper")(upper)
}

plotAUCfacet = function(facet.by.size = c(T, F)){
  Vectorize(function(facet.by.size){
    df = readRDS("../auc/auc.rds")

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
    
    ggsave(paste("../fig/auc-facet-", facet.by.size, ".pdf", sep=""), pl, width = 8, height = 5, dpi = 1600)
  }, "facet.by.size")(facet.by.size)
}

plotAUCcolor = function(jitter = c(T, F)){
  Vectorize(function(jitter){
    df = readRDS("../auc/auc.rds")
    df$mtd = ordered(df$mtd, rev(work.parms("mtd")))
    
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
    
    ggsave(paste("../fig/auc-color-", jitter, ".pdf", sep=""), pl, width = 8, height = 5, dpi = 1600)
  }, "jitter")(jitter)
}

plotFDRfacet = function(facet.direction = c(T, F), who = c("dan", "jarad"), y.axis = c("subtract", "leave"),
reverse.x = c(" ", " no ")){
  grid = expand.grid(who, facet.direction, y.axis, reverse.x)
  
  Vectorize(function(who, facet.direction, y.axis, reverse.x){
    d = readRDS(paste("../fdr/",who, ".rds", sep=""))
  
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

     ggsave(paste("../fig/fdr-facet-",who, "-", facet.direction, "-", y.axis, "-", trim(reverse.x), ".pdf", sep=""), pl, width = 8, height = 5, dpi = 1600)
  })(grid[[1]], grid[[2]], grid[[3]], grid[[4]])
}

plotFDRindiv = function(who = c("dan", "jarad"), mtd = work.parms("mtd"), y.axis = c("subtract", "leave"), 
reverse.x = c(" ", " no ")){
  grid = expand.grid(who, mtd, y.axis, reverse.x)
  
  Vectorize(function(who, mtd, y.axis, reverse.x){
    d = readRDS(paste("../fdr/",who, ".rds", sep=""))
    d = d[d$mtd == mtd,]

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

     ggsave(paste("../fig/fdr-indiv-",who, "-", mtd, "-", y.axis, "-", trim(reverse.x), ".pdf", sep=""), pl, width = 8, height = 5, dpi = 1600)
  })(grid[[1]], grid[[2]], grid[[3]], grid[[4]])
}

makePlots = function(){
  plotExampleROC()
  plotAUCcolor()
  plotAUCfacet()
  plotFDRfacet()
  plotFDRindiv()
}