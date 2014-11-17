for(file in list.files())
  if(file != "setup.R" && file != "workflow.R" && tolower(substr(x = file, start = nchar(file) - 1, stop = nchar(file))) == ".r")
    source(file)

for(dir in c("../auc", "../data/", "../fig",  "../log", "../pvals", "../roc", "../simulations"))
  if(!file.exists(dir))
    dir.create(dir)

host = System$getHostname()
sample.sizes = c(4, 8, 16)

if(host == "impact2.stat.iastate.edu"){ 
  ncpus = 4
  which.datasets = 1:2
} else if(host == "impact3.stat.iastate.edu") {
  ncpus = 6
  which.datasets = 3:6
} else if(host == "impact4.stat.iastate.edu") {
  ncpus = 6
  which.datasets = 7:10
} else {
  ncpus = 2
  which.datasets = 1:2
}