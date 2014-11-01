# Run this file in its containing directory with 
# nohup nice -19 R CMD BATCH workflow.R &

for(file in list.files())
  if(file != "workflow.R" && tolower(substr(x = file, start = nchar(file) - 1, stop = nchar(file))) == ".r")
    source(file)

for(dir in c("../auc",  "../log", "../pvals", "../roc", "../simulations"))
  if(!file.exists(dir))
    dir.create(dir)

host = System$getHostname()

if(host == "impact1.stat.iastate.edu"){
  ncpus = 4
  which.datasets = 1:25
} else if(host == "impact2.stat.iastate.edu"){ 
  ncpus = 5
  which.datasets = 26:65  
} else if(host == "impact3.stat.iastate.edu") {
  ncpus = 2
  which.datasets = 66:75
} else if(host == "impact4.stat.iastate.edu") {
  ncpus = 5
  which.datasets = 76:100
}

#s = proc.time()
#simulate.datasets(which.datasets)
#t = proc.time()
#logfile("Time spent simulating datasets:", t - s)

s = proc.time()
compute.pvals(which.datasets, ncpus)
t = proc.time()
logfile("Time spent computing pvals:", t - s)