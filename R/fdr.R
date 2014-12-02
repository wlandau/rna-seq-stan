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
  .fdr$mtd = mtd
  .fdr$size = size
  .fdr$rep = rep
  .fdr = .fdr[.fdr$cutoff <= upper,]

  saveRDS(.fdr, paste("../fdr/", who, "-", file.name(mtd, size, rep), sep=""))
}

fdr = function(people = c("dan", "jarad"), 
                      mtd = work.parms("mtd"), size = work.parms("size"), rep = work.parms("rep")){

  for(who in people){
    loopify(function(mtd, size, rep){
      print(paste("FDR", who, mtd, size, rep))
      
      d = data.frame(
        ranks = readRDS(paste("../ranks/", file.name(mtd, size, rep), sep="")),
        truth = 1 - abs(readRDS(paste("../simulations/", file.name("truth", size, rep), sep=""))))
      
      d = d[order(d$ranks),]
     
      if(who == "dan")
        .fdr = dan.fdr(d)
      else if(who == "jarad")
        .fdr = jarad.fdr(d)

      save.fdr.1rep(.fdr, who, mtd, size, rep)
    }, mtd, size, rep)

    .fdr = loopify(function(mtd, size, rep){
      readRDS(paste("../fdr/", who, "-", file.name(mtd, size, rep), sep=""))
    }, mtd, size, rep)

    ret = data.frame(
      cutoff = do.call("c", .fdr$cutoff),
      one.minus.cutoff = 1 - do.call("c", .fdr$cutoff),
      fdr = do.call("c", .fdr$fdr),
      fdr.minus.cutoff = do.call("c", .fdr$fdr.minus.cutoff),
      mtd = ordered(do.call("c", lapply(.fdr$mtd, as.character)), levels = rev(work.parms("mtd"))),
      size.short = ordered(do.call("c", .fdr$size), work.parms("size")),
      size = ordered(do.call("c", .fdr$size), labels = paste(work.parms("size"), "samples / group")),
      rep = do.call("c", .fdr$rep)
    )

    ret = ret[is.finite(ret$fdr),]
    saveRDS(ret, paste("../fdr/", who, ".rds", sep=""))
  }
}