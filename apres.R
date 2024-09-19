# Codes of affinity propagation algorithm

# incrementally add cluster points to find the clustering with maximum
# silhouette, version for AP algorithm
# apres is the clustering result of AP algorithm
# ddatax is the distance matrix
apres.silhouette <- function(apres,ddatax){
    
  len_exemplars <- sapply(apres@clusters,length)
  exemplars <- apres@exemplars
  exemplars <- exemplars[order(len_exemplars,decreasing = TRUE)]
  maxsi <- 0
  si.ind <- 0
  for(i in 2:length(exemplars)){
    temp_exemplars <- exemplars[1:i]
    d2 <- ddatax[temp_exemplars,]
    result <- apply(d2,2,which.min)
    if(length(unique(result)) < 2) datasi <- 0
    else{
      datasi <- silhouette(result,ddatax)
      datasi <- mean(datasi[,3])  
    }
    if(datasi > maxsi){
      maxsi <- datasi
      si.ind <- i
    } 
  }
  
  temp_exemplars <- exemplars[1:si.ind]
  d2 <- ddatax[temp_exemplars,]
  result <- apply(d2,2,which.min)
  
  return(list(si=maxsi,si.ind=si.ind,clustering=result))
}