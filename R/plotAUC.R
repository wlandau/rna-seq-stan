plotAUC = function(file = "../auc/auc.rds"){
  a = readRDS(file)

  a$AUC = as.numeric(as.vector(a$AUC))
  a$Method = factor(a$Method, levels = c("edgeR", "ShrinkBayes", "baySeq"))
  a$SampleSize = factor(a$SampleSize, levels = c(4, 8, 16))
  a$SampleSize = revalue(a$SampleSize, c("4" = "4 samples per group", "8" = "8 samples per group", "16" = "16 samples per group"))

  df = ddply(a, .variables = .(Method, SampleSize), .fun = function(x){
    c(mean(x$AUC), quantile(x$AUC, c(0.025, 0.975)))
  })

  colnames(df) = c("Method", "SampleSize", "Mean", "Lower", "Upper")

  df = rbind(df[df$Method == "edgeR",], df[df$Method == "baySeq",], df[df$Method == "ShrinkBayes",])
  df =  df[rep(1:dim(df)[1], each = 10),]

  a = cbind(a, Mean = df$Mean, Lower = df$Lower, Upper = df$Upper)

  pl = ggplot(a, aes(x = Method, y = AUC), na.rm=TRUE) + 
         geom_line(aes(group = Replicate), alpha = .5) +  
         geom_point(aes(y = Mean, group = Replicate), color = "blue", cex = 2) + 
         geom_errorbar(aes(ymin = Lower, ymax = Upper), color = "blue") +
         facet_grid(~ SampleSize, scales = "free_x") + 
         xlab("\nMethod") +
         ylab("Area under Receiver Operating Characteristic (ROC) curve \n(false positive rate < 0.1)\n") + 
         theme(axis.text = element_text(family = "Helvetica", colour = 'black'),
                    axis.text.x = element_text(angle = -45, hjust = 0),
                    legend.position="none",
                    panel.background = element_rect(fill='white'),
                    panel.border = element_rect(color="black", fill = NA),
                    panel.grid.major = element_line(color="lightgray"),
                    text = element_text(family = "Helvetica", colour= "black"))
    
    ggsave("../fig/auc.pdf", pl, width = 8, height = 5, dpi = 1600)
}