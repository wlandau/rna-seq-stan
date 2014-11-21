plotExampleROC = function(upper = 1e-1, pkgs = c("edgeR", "baySeq", "ShrinkBayes", "stan_corr", "stan"), sample.sizes = c(4, 8, 16)){

  df = NULL
  for(pkg in pkgs){
    for(s in sample.sizes){
      r = readRDS(paste("../roc/", pkg, "1-", s, ".rds", sep=""))
      r = r[r$fpr <= upper,]
      df = rbind(df, cbind(r$fpr, r$tpr, s, pkg))
    }
  }

  df = data.frame(df)
  colnames(df) = c("FPR", "TPR", "SampleSize", "Method")

  df$FPR = as.numeric(as.vector(df$FPR))
  df$TPR = as.numeric(as.vector(df$TPR))
  df$Method = factor(df$Method, levels = c("edgeR", "ShrinkBayes", "baySeq", "stan_corr", "stan"))
  df$SampleSize = factor(df$SampleSize, levels = c(4, 8, 16))
  df$SampleSize = revalue(df$SampleSize, c("4" = "4 samples per group", "8" = "8 samples per group", "16" = "16 samples per group"))

  pl = ggplot(df, aes(x = FPR, y = TPR), na.rm=TRUE) + #+ geom_line() +
         geom_abline() +  
         geom_line(aes(color = Method)) + 
         xlab("\nFalse positive rate (FPR)") + 
         ylab("True positive rate (TPR)\n") +
         facet_grid(~ SampleSize) + 
         theme(axis.text = element_text(family = "Helvetica", color = 'black'),
                    axis.text.x = element_text(angle = -45, hjust = 0),
                    panel.background = element_rect(fill='white'),
                    panel.border = element_rect(color="black", fill = NA),
                    panel.grid.major = element_line(color="lightgray"),
                    text = element_text(family = "Helvetica", colour= "black"))

  ggsave(paste("../fig/exampleROC", upper, ".pdf", sep=""), pl, width = 8, height = 5, dpi = 1600)

#  pl = pl + theme(legend.position="none")
#  ggsave(paste("../fig/exampleROC-nolegend", upper, ".pdf", sep=""), pl, width = 8, height = 5, dpi = 1600)
}

plotAUCfacet = function(file = "../auc/auc.rds", facet.by = "SampleSize"){
  a = readRDS(file)

  a$AUC = as.numeric(as.vector(a$AUC))
  a$Method = factor(a$Method, levels = c("edgeR", "ShrinkBayes", "baySeq", "stan_corr", "stan"))
  a$SampleSize = factor(a$SampleSize, levels = c(4, 8, 16))

  df = ddply(a, .variables = .(Method, SampleSize), .fun = function(x){
    c(mean(x$AUC), mean(x$AUC) + c(-1, 1) * 1.96  * sd(x$AUC)/sqrt(length(x$AUC)))
  })

  colnames(df) = c("Method", "SampleSize", "Mean", "Lower", "Upper")

  df = rbind(df[df$Method == "edgeR",], df[df$Method == "baySeq",], df[df$Method == "ShrinkBayes",], 
                  df[df$Method == "stan_corr",], df[df$Method == "stan",])
  df =  df[rep(1:dim(df)[1], each = 10),]

  a = cbind(a, Mean = df$Mean, Lower = df$Lower, Upper = df$Upper)

  if(facet.by == "SampleSize"){
    a$SampleSize = revalue(a$SampleSize, c("4" = "4 samples per group", "8" = "8 samples per group", "16" = "16 samples per group"))
    pl = ggplot(a, aes(x = Method, y = AUC), na.rm=TRUE) +
           facet_grid(~ SampleSize, scales = "free_x") + 
           theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
           xlab("\nMethod")
  } else {
    pl = ggplot(a, aes(x = SampleSize, y = AUC), na.rm=TRUE) +
           facet_grid(~ Method, scales = "free_x") +
           xlab("\nSample Size per Treatment Group")
  }
   pl = pl +
         geom_line(aes(group = Replicate), alpha = .5) +  
 #        geom_point(aes(y = Mean, group = Replicate), color = "blue", fill = "white", pch = 3, cex = 6) + 
 #        geom_errorbar(aes(ymin = Lower, ymax = Upper), color = "blue", position = "identity", width = .5) +
         geom_crossbar(aes(ymin = Lower, y = Mean, ymax = Upper), fatten = 1, width = .5, color = "blue") + 
         ylab("Area under Receiver Operating Characteristic\n(ROC) curve (false positive rate < 0.1)\n") + 
         theme(axis.text = element_text(family = "Helvetica", colour = 'black'),
                    legend.position="none",
                    panel.background = element_rect(fill='white'),
                    panel.border = element_rect(color="black", fill = NA),
                    panel.grid.major = element_line(color="lightgray"),
                    text = element_text(family = "Helvetica", colour= "black"))
    
    ggsave(paste("../fig/auc-", facet.by, ".pdf", sep=""), pl, width = 8, height = 5, dpi = 1600)
}

plotAUC = function(file = "../auc/auc.rds", jitter = T){
  a = readRDS(file)

  a$AUC = as.numeric(as.vector(a$AUC))
  a$Method = factor(a$Method, levels = c("edgeR", "ShrinkBayes", "baySeq", "stan_corr", "stan"))
  a$SampleSize = factor(a$SampleSize, levels = c(4, 8, 16))

  df = ddply(a, .variables = .(Method, SampleSize), .fun = function(x){
    c(mean(x$AUC), mean(x$AUC) + c(-1, 1) * 1.96  * sd(x$AUC)/sqrt(length(x$AUC)))
  })

  colnames(df) = c("Method", "SampleSize", "Mean", "Lower", "Upper")

  df = rbind(df[df$Method == "edgeR",], df[df$Method == "baySeq",], df[df$Method == "ShrinkBayes",], 
         df[df$Method == "stan_corr",], df[df$Method == "stan",])
  df =  df[rep(1:dim(df)[1], each = 10),]

  a = cbind(a, Mean = df$Mean, Lower = df$Lower, Upper = df$Upper)


  pl = ggplot(a, aes(x = SampleSize, y = AUC, color = Method), na.rm=TRUE) +
#         geom_point(aes(y = Mean), pch = 3, cex = 6) + 
 #        geom_errorbar(aes(ymin = Lower, ymax = Upper), width = .25) +
         geom_crossbar(aes(ymin = Lower, y = Mean, ymax = Upper), fatten = 1, width = .25) + 
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

    ggsave(paste("../fig/auc", jitter, ".pdf", sep=""), pl, width = 8, height = 5, dpi = 1600)
}

plotFDR = function(type = "dan", facet.direction = T){
  d = readRDS(paste("../fdr/", type, ".rds", sep=""))
  for(i in 1:dim(d)[2])
    d[[i]][is.na(d[[i]])] = 0
 

  pl = ggplot(d, aes(x = cutoff, y = FDRminusCutoff)) + 
         geom_abline(slope = 0, intercept = 0, color = "blue") + 
         geom_line(aes(group = rep)) + 
         xlab("\nCutoff (p-value or posterior probability)") + 
         ylab("FDR - cutoff\n") + 
         facet_grid(pkg ~ sample.size)

  if(facet.direction)
    pl = pl + facet_grid(pkg ~ sample.size)
  else
    pl = pl + facet_grid(sample.size ~ pkg)

   ggsave(paste("../fig/fdr-", type, "-", facet.direction, ".pdf", sep=""), pl, width = 8, height = 5, dpi = 1600)
}

makePlots = function(){
  plotExampleROC(1e-1)
  plotExampleROC(1)
  plotAUC(jitter = F)
  plotAUC(jitter = T) 
  plotAUCfacet(facet.by = "SampleSize")
  plotAUCfacet(facet.by = "Method")
  plotFDR("dan", T)
  plotFDR("dan", F)
  plotFDR("jarad", T)
  plotFDR("jarad", F)
}