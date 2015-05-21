for(file in list.files())
  if(file != "setup.R" && file != "workflow.R" && tolower(substr(x = file, start = nchar(file) - 1, stop = nchar(file))) == ".r")
    source(file)

for(dir in paste(work.parms("path"), c("auc", "fdr", "log", "priors", "ranks", "roc"), sep=""))
  if(!file.exists(dir))
    dir.create(dir)

host = System$getHostname()

if(host == "impact1.stat.iastate.edu"){
  ncups = 12
  reps.on.server = 1:3
} else if(host == "impact2.stat.iastate.edu"){ 
  ncpus = 6
  reps.on.server = 1:4
} else if(host == "impact3.stat.iastate.edu"){
  ncpus = 6
  reps.on.server = 1:10
} else if(host == "impact4.stat.iastate.edu"){
  ncpus = 6
  reps.on.server = 8:10
} else {
  ncpus = 2
  reps.on.server = 1:10
}