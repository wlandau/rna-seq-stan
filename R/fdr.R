dan.fdr = function(d){
  apply(d, 2, function(x){cumsum(x)/1:length(x)})
}

jarad.fdr = function(d){
  t(sapply(1:(dim(d)[1] - 100), function(n){apply(d[n:(n + 100),], 2, mean)}))
}

save.fdr.1rep = function(.fdr, who, mtd, size, rep, upper = 0.15, npts = 100){
  .fdr = as.data.frame(.fdr)
  colnames(.fdr) = c("cutoff", "fdr")
  apx = approx(x = .fdr$cutoff, .fdr$fdr, xout = (1:npts/npts)*0.15)
  .fdr = data.frame(cutoff = apx$x, fdr = apx$y)

  .fdr$fdr.minus.cutoff = .fdr$fdr - .fdr$cutoff
  .fdr$one.minus.cutoff = 1 - .fdr$cutoff
  .fdr$mtd = mtd
  .fdr$size.short = ordered(size, work.parms("size"))
  .fdr$rep = ordered(rep, work.parms("rep"))

  .fdr = .fdr[.fdr$cutoff <= upper,]

  saveRDS(.fdr, paste(work.parms("path"), "fdr/", who, "-", file.name(mtd, size, rep), sep=""))
}

collect.fdr = function(who = "dan"){
  .fdr.list = loopify(function(mtd, size, rep){
    n = paste(work.parms("path"), "fdr/", who, "-", file.name(mtd, size, rep), sep="")
    if(file.exists(n))
      readRDS(n)
  }, work.parms("mtd"), work.parms("size"), work.parms("rep"))

  .fdr = list()
  for(n in names(.fdr.list))
    .fdr[[n]] = do.call("c", lapply(.fdr.list[[n]], as.vector))
  .fdr = as.data.frame(.fdr)

  .fdr$size = ordered(paste(as.character(.fdr$size.short), "samples / group"), labels = paste(work.parms("size"), "samples / group"))
  .fdr$mtd = ordered(.fdr$mtd, levels = rev(work.parms("mtd")))

  .fdr = .fdr[is.finite(.fdr$fdr),]
  saveRDS(.fdr, paste(work.parms("path"), "fdr/", who, ".rds", sep=""))
}

fdr = function(people = c("dan", "jarad"), 
                      mtd = work.parms("mtd"), size = work.parms("size"), rep = work.parms("rep")){

  for(who in people){
    loopify(function(mtd, size, rep){
      print(paste("FDR", who, mtd, size, rep))
      
      d = data.frame(
        ranks = readRDS(paste(work.parms("path"), "ranks/", file.name(mtd, size, rep), sep="")),
        truth = 1 - abs(readRDS(paste(work.parms("path"), "truth/", file.name("truth", size, rep), sep=""))))
      
      d = d[order(d$ranks),]
     
      if(who == "dan")
        .fdr = dan.fdr(d)
      else if(who == "jarad")
        .fdr = jarad.fdr(d)

      save.fdr.1rep(.fdr, who, mtd, size, rep)
    }, mtd, size, rep)
    
    collect.fdr(who)
  }
}