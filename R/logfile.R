logfile = function(...){
  fileConn<-file(paste("../log/log-", System$getHostname(), ".txt", sep=""), open="a")
  writeLines(paste(...), fileConn)
  close(fileConn)
}