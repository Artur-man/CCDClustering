# codes for Fuzzy C-means clustering

# the clustering function for the fuzzy cmeans clustering
# library flexclust is needed
# centers: found centers of cmeans algorithm
# datax: data
# cls: actual classes 
cmeans.cluster<- function(centers,datax){
  
  # find closest points to centers, check their distances
  centtodata <- as.matrix(dist2(centers,datax))
  result <- apply(centtodata,2,which.min)
  
  return(result)
}

# the clustering function for the fuzzy cmeans clustering
# library flexclust is needed
# centers: found centers of cmeans algorithm
# datax: data
# cls: actual classes 
cmeans.valid<- function(centers,datax,cls){
  
  # find closest points to centers, check their distances
  centtodata <- as.matrix(dist2(centers,datax))
  ind.neigh <- apply(centtodata,1,which.min)
  cl <- cls[ind.neigh]
  
  result <- apply(centtodata,2,function(x){
    ind <- cl[which.min(x)]
    return(ind)
  })
  
  return(result)
}