for(file in list.files())
  if(file != "setup.R" && file != "workflow.R" && tolower(substr(x = file, start = nchar(file) - 1, stop = nchar(file))) == ".r")
    source(file)

for(dir in c("../auc", "../data/", "../fdr", "../fig",  "../log", "../priors", "../pvals", "../roc", "../simulations"))
  if(!file.exists(dir))
    dir.create(dir)

host = System$getHostname()

if(host == "impact1.stat.iastate.edu"){
  ncups = 12
  reps.on.server = 1:2
} else if(host == "impact2.stat.iastate.edu"){ 
  ncpus = 6
  reps.on.server = 1:4
} else if(host == "impact3.stat.iastate.edu"){
  ncpus = 2
  reps.on.server = 5:6
} else if(host == "impact4.stat.iastate.edu"){
  ncpus = 6
  reps.on.server = 7:10
} else {
  ncpus = 2
  reps.on.server = 1:2
}