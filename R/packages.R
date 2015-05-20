is.installed = function(pkg){
  is.element(pkg, installed.packages()[,1])
} 

safeLoad = Vectorize(function(pkg){
  if(is.installed(pkg))
    library(pkg, character.only = T)
  else
    print(paste(pkg, "package not found."))
}, "pkg")

safeLoad(c(
  "ggplot2",
  "plyr",
  "pracma",
  "R.utils",
  "edgeR",
  "baySeq",
  "ShrinkBayes",
  "DESeq2",
  "heterosis"
))