for(file in list.files())
  if(file != "setup.R" && file != "workflow.R" && tolower(substr(x = file, start = nchar(file) - 1, stop = nchar(file))) == ".r")
    source(file)

for(dir in c("../auc", "../data/", "../fig",  "../log", "../pvals", "../roc", "../simulations"))
  if(!file.exists(dir))
    dir.create(dir)

host = System$getHostname()

if(host == "impact1.stat.iastate.edu"){
  ncpus = 4
  which.datasets = 1:33
} else if(host == "impact2.stat.iastate.edu"){ 
  ncpus = 4
  which.datasets = 1:25
} else if(host == "impact3.stat.iastate.edu") {
  ncpus = 6
  which.datasets = 26:62
} else if(host == "impact4.stat.iastate.edu") {
  ncpus = 6
  which.datasets = 63:100
} else {
  ncpus = 2
  which.datasets = 1:2
}