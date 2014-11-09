plotExampleROC = function(file = "../roc/edgeR1-16.rds", upper = 1e-1){
  r = readRDS(file)
  r = r[r$fpr <= upper,]

  pl = ggplot(r, aes(x = fpr, y = tpr), na.rm=TRUE) + 
    geom_line() +
    xlab("\nFalse positive rate (FPR)") + 
    ylab("True positive rate (TPR)\n") +
    theme(text = element_text(family = "Helvetica", size = 12, colour= "black", face="plain"),
          legend.title= element_text(family = "Helvetica", size=12, face="plain"),
          axis.text = element_text(family = "Helvetica", colour = 'black'),
          axis.line=element_line())

  pl = pl + theme(panel.background = element_rect(fill='white'),
                  panel.grid.major = element_line(color="lightgray"),
                  panel.border = element_rect(color="black", fill = NA))

  ggsave("../fig/exampleROC.pdf", pl, width = 5, height = 5, dpi = 1600)
}