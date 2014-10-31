trim_genes <- function(counts, group, geneid, z.allow.tr = 3, mean.thr = exp(1)){
  #counts is a num_gene x num_sample matrix
  #group is a num_sample length factor mapping columns to groups
  #geneid is a num_gene length vector rows to genes, i.e. 1:G or descriptive names
  #z.allow.tr is maximum zeros allowable for any treatment(gene)
  #mean.thr is a minimum average count for a gene
  G = length(geneid)
  if(G != nrow(counts)) stop("length(geneid) must match nrow(counts)")
  if(length(group) != ncol(counts)) stop("length(group) must match ncol(counts)")
  key_for_id = cbind(geneid,1:G)
  counts = as.matrix(counts)
  
  #Which genes have a treatment with excessive zeros?
  flag_zeros <- NULL
  for(i in group){
    flags  = apply(counts, MARGIN=1,FUN= function(x){
        
        sum(x[which(group==i)]==0) > z.allow.tr       #indicator for group i (for all genes) 
        })
        
    flag_zeros = c(flag_zeros, which(flags))
  }                  
    
  #Which genes have too low average expression?
  flags2 = apply(counts, MARGIN=1, FUN = function(x) mean(x) < mean.thr)
  flag_low = which(flags2)
  
  allflags = sort(unique(c(flag_zeros,flag_low)))
  
  
  return(list(counts[-allflags,],geneid[-allflags]))
}