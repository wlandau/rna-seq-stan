logfile = function(...){
  fileConn<-file(paste(work.parms("path"), "log/log-", System$getHostname(), ".txt", sep=""), open="a")
  writeLines(paste(...), fileConn)
  close(fileConn)
}