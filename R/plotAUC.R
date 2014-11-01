plotAUC = function(file = "../auc/auc.rds"){
  a = readRDS(file)

  pl = ggplot(a, aes(x = pkg, y = auc), na.rm=TRUE) + 
           geom_boxplot(aes(x = pkg, y = auc), fill = "gray") + 
           xlab("\nR package") + 
           ylab("Area under Receiver Operating Characteristic (ROC) curve \n(false positive rate < 0.1)\n") +
      theme(axis.text.x = element_text(family = "Helvetica"),
            text = element_text(family = "Helvetica", colour= "black"),
            axis.text = element_text(family = "Helvetica", colour = 'black'),
            axis.line=element_line(),
            legend.position="none")
    
    pl = pl + theme(panel.background = element_rect(fill='white'),
                    panel.grid.major = element_line(color="lightgray"),
                    panel.border = element_rect(color="black", fill = NA)) +
              ylab("Area under ROC curve (false positive rate < 0.1)\n")
    
    ggsave("../fig/auc.pdf", pl, width = 5, height = 5, dpi = 1600)
}