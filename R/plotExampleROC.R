plotExampleROC = function(upper = 1e-1, pkgs = c("edgeR", "baySeq", "ShrinkBayes"), sample.sizes = c(4, 8, 16)){

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
  df$Method = factor(df$Method, levels = c("edgeR", "ShrinkBayes", "baySeq"))
  df$SampleSize = factor(df$SampleSize, levels = c(4, 8, 16))
  df$SampleSize = revalue(df$SampleSize, c("4" = "4 samples per group", "8" = "8 samples per group", "16" = "16 samples per group"))

  pl = ggplot(df, aes(x = FPR, y = TPR), na.rm=TRUE) + #+ geom_line() +
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

  ggsave("../fig/exampleROC.pdf", pl, width = 8, height = 5, dpi = 1600)

  pl = pl + theme(legend.position="none")

  ggsave("../fig/exampleROC-nolegend.pdf", pl, width = 8, height = 5, dpi = 1600)
}