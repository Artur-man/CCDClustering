# Codes of RK-CCD clustering 

# ccd clustering that find the dominating set with greedy algorithm
# the boxes are given with a K function estimation
# old name = ccd3.clustering
# datax, the data set
# low.num, lowest number of a box cardinality 
# r.seq, the number of radii used for each window, proximity region
# simul, provided simulations, if null, compute
rccd.clustering <- function(datax,low.num,r.seq,
                            dom.method="greedy",simul=NULL){
  
  # info
  nr <- nrow(datax)
  nc <- ncol(datax)
  ddatax <- as.matrix(dist(datax))
  mac.eps <- .Machine$double.eps
  
  # find all radii, and the dominating set 
  ccd.info <- ccd.Kest(datax,ddatax,low.num,r.seq,simul)

  r <- ccd.info$R
  M <- matrix(as.integer(ddatax < r+mac.eps), length(r))
  score <- rowSums(M)

  # domination method is provided, find the dominating set
  # if not, take all radii
  if(!is.null(dom.method)){
    if(dom.method=="greedy") D <- dominate.mat.greedy(M)
    if(dom.method=="greedy2") D <- dominate.mat.greedy2(M)
    if(dom.method=="ks") D <- dominate.mat.ks(M,ks)
  } else {
    D <- 1:nr
  }
  R <- r[D]

  # find the intersection graph and its dominating set
  # the intersection graph is based on whether two radii catched the same points
  M <- M[D,]
  M.dom <- diag(T,nrow(M))
  for(i in 1:(nrow(M.dom)-1)){
    for(j in (i+1):nrow(M.dom)){
      temp <- (M[i,] & M[j,])
      if(any(temp)) M.dom[i,j] <- M.dom[j,i] <- T
    }
  }
  
  # dominating set of the intersection graph
  D.ind <- dominate.mat.ks(M.dom,score[D])
  
  # keep the 2nd layer of dominating sets 
  Int.D=D[D.ind]
  Int.R=R[D.ind]
  
  return(list(Int.D=Int.D,Int.R=Int.R,D=D,R=R))
}

# ccd clustering that find the dominating set with greedy alg
# translation correction as edge correction
# the boxes are are given with a K function estimation
# old name = ccd3.clustering_correct
# datax, the data set
# low.num, lowest number of a box cardinality 
# r.seq, the number of radii used for each window, proximity region
# simul, provided simulations, if null, compute
rccd.clustering_correct <- function(datax,low.num,r.seq,
                            dom.method="greedy",simul=NULL){
  
  # info
  nr <- nrow(datax)
  nc <- ncol(datax)
  ddatax <- as.matrix(dist(datax))
  mac.eps <- .Machine$double.eps
  
  # find all radii, and the dominating set 
  ccd.info <- ccd.Kest.edge(datax,ddatax,low.num,r.seq,simul)
  
  r <- ccd.info$R
  M <- matrix(as.integer(ddatax < r+mac.eps), length(r))
  score <- rowSums(M)
  
  # domination method is provided, find the dominating set
  # if not, take all radii
  if(!is.null(dom.method)){
    if(dom.method=="greedy") D <- dominate.mat.greedy(M)
    if(dom.method=="greedy2") D <- dominate.mat.greedy2(M)
    if(dom.method=="ks") D <- dominate.mat.ks(M,ks)
  } else {
    D <- 1:nr
  }
  R <- r[D]
  
  # find the intersection graph and its dominating set
  # the intersection graph is based on whether two radii catched the same points
  MD <- M[D,]
  M.dom <- diag(T,nrow(MD))
  for(i in 1:(nrow(M.dom)-1)){
    for(j in (i+1):nrow(M.dom)){
      temp <- (MD[i,] & MD[j,])
      if(any(temp)) M.dom[i,j] <- M.dom[j,i] <- T
    }
  }
  
  # dominating set of the intersection graph
  D.ind <- dominate.mat.ks(M.dom,score[D])
  
  # keep the 2nd layer of dominating sets 
  Int.D=D[D.ind]
  Int.R=R[D.ind]

  # get the catching info of all Dominated graphs 
  MDInt <- M[Int.D,,drop=FALSE]
  # catch <- apply(MDInt,1,function(x){
  #               temp <- ddatax[as.logical(x),as.logical(x),drop = FALSE]
  #               diag(temp) <- 0
  #               return(sum(temp)/(sum(x)^2-sum(x)))
  # })
  catch <- rowSums(MDInt)/(Int.R^2)
  
  
  return(list(Int.D=Int.D,Int.R=Int.R,D=D,R=R, catch=catch))
}

# ccd clustering that find the dominating set with greedy alg
# translation correction as edge correction
# the boxes are are given with a K function estimation
# a quantile is selected from the simulated envelopes 
# old name = ccd3.clustering_correct_quantile
# datax, the data set
# low.num, lowest number of a box cardinality 
# r.seq, the number of radii used for each window, proximity region
# dom.method, method for finding the dominating set
# quan, quantile to be used for confidence interval
# simul, provided simulations, if null, compute
rccd.clustering_correct_quantile <- function(datax,low.num,r.seq,
                                    dom.method="greedy", quan, simul=NULL){
  
  # info
  nr <- nrow(datax)
  nc <- ncol(datax)
  ddatax <- as.matrix(dist(datax))
  mac.eps <- .Machine$double.eps
  
  # find all radii, and the dominating set 
  ccd.info <- ccd.Kest.edge.quantile(datax,ddatax,low.num,r.seq,quan, simul)
  
  r <- ccd.info$R
  M <- matrix(as.integer(ddatax < r+mac.eps), length(r))
  score <- rowSums(M)
  
  # domination method is provided, find the dominating set
  # if not, take all radii
  if(!is.null(dom.method)){
    if(dom.method=="greedy") D <- dominate.mat.greedy(M)
    if(dom.method=="greedy2") D <- dominate.mat.greedy2(M)
    if(dom.method=="ks") D <- dominate.mat.ks(M,ks)
  } else {
    D <- 1:nr
  }
  R <- r[D]
  
  # find the intersection graph and its dominating set
  # the intersection graph is based on whether two radii catched the same points
  MD <- M[D,]
  M.dom <- diag(T,nrow(MD))
  for(i in 1:(nrow(M.dom)-1)){
    for(j in (i+1):nrow(M.dom)){
      temp <- (MD[i,] & MD[j,])
      if(any(temp)) M.dom[i,j] <- M.dom[j,i] <- T
    }
  }
  
  # dominating set of the intersection graph
  D.ind <- dominate.mat.ks(M.dom,score[D])
  
  # keep the 2nd layer of dominating sets 
  Int.D=D[D.ind]
  Int.R=R[D.ind]
  
  # get the catching info of all Dominated graphs 
  MDInt <- M[Int.D,,drop=FALSE]
  catch <- rowSums(MDInt)/(Int.R^2)
  
  # return
  return(list(Int.D=Int.D,Int.R=Int.R,D=D,R=R, catch=catch, M=M))
}

# ccd clustering that find the dominating set with greedy alg
# translation correction as edge correction
# the boxes are given with a K function estimation
# the dominating set is found by maximizing the local intensity
# old name = ccd3.clustering_correct_spatial
# datax, the data set
# low.num, lowest number of a box cardinality 
# r.seq, the number of radii used for each window, proximity region
# r.lim, limit of the maximum radii 
# simul, provided simulations, if null, compute
rccd.clustering_correct_spatial <- function(datax,low.num,r.seq,
                                    dom.method="greedy",simul=NULL, r.lim){
  
  # info
  nr <- nrow(datax)
  nc <- ncol(datax)
  ddatax <- as.matrix(dist(datax))
  mac.eps <- .Machine$double.eps
  
  # find all radii, and the dominating set 
  ccd.info <- ccd.Kest.edge(datax,ddatax,low.num,r.seq,simul,r.lim)
  
  r <- ccd.info$R
  M <- matrix(as.integer(ddatax < r+mac.eps), length(r))
  score <- rowSums(M)
  
  # domination method is provided, find the dominating set
  # if not, take all radii
  if(!is.null(dom.method)){
    if(dom.method=="greedy") D <- dominate.mat.greedy(M)
    if(dom.method=="greedy2") D <- dominate.mat.greedy2(M)
    if(dom.method=="ks") D <- dominate.mat.ks(M,ks)
  } else {
    D <- 1:nr
  }
  R <- r[D]
  
  # find the intersection graph and its dominating set
  # the intersection graph is based on whether two radii catched the same points
  MD <- M[D,]
  M.dom <- diag(T,nrow(MD))
  for(i in 1:(nrow(M.dom)-1)){
    for(j in (i+1):nrow(M.dom)){
      temp <- (MD[i,] & MD[j,])
      if(any(temp)) M.dom[i,j] <- M.dom[j,i] <- T
    }
  }
  
  # dominating set of the intersection graph
  D.ind <- dominate.mat.ks(M.dom,score[D])
  
  # keep the 2nd layer of dominating sets 
  Int.D=D[D.ind]
  Int.R=R[D.ind]
  
  # get the catching info of all Dominated graphs 
  MDInt <- M[Int.D,,drop=FALSE]
  # catch <- apply(MDInt,1,function(x){
  #               temp <- ddatax[as.logical(x),as.logical(x),drop = FALSE]
  #               diag(temp) <- 0
  #               return(sum(temp)/(sum(x)^2-sum(x)))
  # })
  catch <- rowSums(MDInt)/(Int.R^2)
  
  
  return(list(Int.D=Int.D,Int.R=Int.R,D=D,R=R, catch=catch))
}

# ccd clustering that find the dominating set with greedy alg
# translation correction as edge correction
# the boxes are are given with a K function estimation
# old name = ccd3.clustering_correct
#
# fast version with binary search
# 
# datax, the data set
# low.num, lowest number of a box cardinality 
# r.seq, the number of radii used for each window, proximity region
# simul, provided simulations, if null, compute
rccd.clustering_correct_fast <- function(datax,low.num,r.seq,
                                    dom.method="greedy",simul=NULL){
  
  # info
  nr <- nrow(datax)
  nc <- ncol(datax)
  ddatax <- as.matrix(dist(datax))
  mac.eps <- .Machine$double.eps
  
  # find all radii, and the dominating set 
  ccd.info <- ccd.Kest.edge_fast(datax,ddatax,low.num,r.seq,simul)
  
  r <- ccd.info$R
  M <- matrix(as.integer(ddatax < r+mac.eps), length(r))
  score <- rowSums(M)
  
  # domination method is provided, find the dominating set
  # if not, take all radii
  if(!is.null(dom.method)){
    if(dom.method=="greedy") D <- dominate.mat.greedy(M)
    if(dom.method=="greedy2") D <- dominate.mat.greedy2(M)
    if(dom.method=="ks") D <- dominate.mat.ks(M,ks)
  } else {
    D <- 1:nr
  }
  R <- r[D]
  
  # find the intersection graph and its dominating set
  # the intersection graph is based on whether two radii catched the same points
  MD <- M[D,]
  M.dom <- diag(T,nrow(MD))
  for(i in 1:(nrow(M.dom)-1)){
    for(j in (i+1):nrow(M.dom)){
      temp <- (MD[i,] & MD[j,])
      if(any(temp)) M.dom[i,j] <- M.dom[j,i] <- T
    }
  }
  
  # dominating set of the intersection graph
  D.ind <- dominate.mat.ks(M.dom,score[D])
  
  # keep the 2nd layer of dominating sets 
  Int.D=D[D.ind]
  Int.R=R[D.ind]
  
  # get the catching info of all Dominated graphs 
  MDInt <- M[Int.D,,drop=FALSE]
  # catch <- apply(MDInt,1,function(x){
  #               temp <- ddatax[as.logical(x),as.logical(x),drop = FALSE]
  #               diag(temp) <- 0
  #               return(sum(temp)/(sum(x)^2-sum(x)))
  # })
  catch <- rowSums(MDInt)/(Int.R^2)
  
  
  return(list(Int.D=Int.D,Int.R=Int.R,D=D,R=R, catch=catch))
}

# the validation function for the clusters
# old name = ccd3.clustering.valid
# graph is the digraph CCD
# ddatax is the distance matrix
# cls is the actual classes 
rccd.clustering.valid <- function(graph,ddatax,cls){
  
  ncls <- length(unique(cls))
  D <- graph$Int.D
  R <- graph$Int.R
  cl <- cls[D]
  
  ddx <- matrix(ddatax[,D],ncol=length(D))
  
  # which point is in which ball of dominators
    result <- t(apply(ddx,1,function(x){
      ind <- which.min(x/R)
      return(cl[ind])
    }))
    
  if(length(D) < 2) result <- rep(1,nrow(ddx))
  return(result)
}

# the validation function for the clusters, does not need the actual cluster labels
# old name = ccd3.clustering.nonvalid
# graph is the digraph CCD
# ddatax is the distance matrix
# cls is the actual classes 
rccd.clustering.nonvalid <- function(graph,ddatax){
  
  D <- graph$Int.D
  R <- graph$Int.R
  cl <- 1:length(D)
  ddx <- matrix(ddatax[,D],ncol=length(D))
  
  # which point is in which ball of dominators
  ddx <- ddx + .Machine$double.eps
  result <- t(apply(ddx,1,function(x){
    ind <- which.min(x/R)
    return(cl[ind])
  }))
  return(result)
}

# incrementally add cluster points to find the clustering with maximum silhouette
# old name ccd3.silhouette
# graph is the digraph CCD
# ddatax is the distance matrix
# cls is the actual classes, (deprecated)
# ind is the set of indices of clusterings to be checked
# lenDlimit is the maximum index of the dominating point to be checked for silhouette
rccd.silhouette <- function(graph,ddatax,cls,ind=NULL, lenDlimit=NULL){
  
  lenD <- length(graph$Int.D)
  if(!is.null(lenDlimit)){
    if(lenDlimit < lenD) lenD <- lenDlimit
  }
  dgraph <- graph
  
  if(is.null(ind)) ind <- 1:nrow(ddatax)
  
  maxsi <- 0
  si.ind <- 0
  si.list <- rep(0,lenD)
  for(i in 2:lenD){
    dgraph$Int.D <- graph$Int.D[1:i]
    dgraph$Int.R <- graph$Int.R[1:i]
    result <- rccd.clustering.nonvalid(dgraph,ddatax)
    if(length(unique(result[ind])) < 2) datasi <- 0
    else{
      datasi <- silhouette(result[ind],ddatax[ind,ind])
      datasi <- mean(datasi[,3])  
      si.list[i] <- datasi
    }
    if(datasi > maxsi){
      maxsi <- datasi
      si.ind <- i
    } 
  }
  
  return(list(si=maxsi,si.ind=si.ind, si.list = si.list))
}

# rccd clustering for arbitrary shapes of data
# the balls are given with a K function estimation
# old name = ccd4.clustering
# datax, the data set
# low.num, lowest number of a box cardinality 
# r.seq, the number of radii used for each window, proximity region
# method, the method for dominating set algorithm
# simul, the simulation results if provided 
rccd.clustering.connected <- function(datax,low.num,r.seq,dom.method="greedy",simul=NULL){
  
  # info
  nr <- nrow(datax)
  nc <- ncol(datax)
  ddatax <- as.matrix(dist(datax))
  mac.eps <- .Machine$double.eps
  
  # find all radii, and the dominating set 
  ccd.info <- ccd.Kest(datax,ddatax,low.num,r.seq,simul)
  
  r <- ccd.info$R
  M <- matrix(as.integer(ddatax < r+mac.eps), length(r))
  diag(M) <- 1
  score <- rowSums(M)
  
  # domination method is provided, find the dominating set
  # if not, take all radii
  if(!is.null(dom.method)){
    if(dom.method=="greedy") D <- dominate.mat.greedy(M)
    if(dom.method=="greedy2") D <- dominate.mat.greedy2(M)
    if(dom.method=="ks") D <- dominate.mat.ks(M,ks)
  } else {
    D <- 1:nr
  }
  R <- r[D]
  
  # find the intersection graph and its dominating set
  # the intersection graph is based on whether two radii catched the same points
  M <- M[D,]
  M.dom <- diag(T,nrow(M))
  for(i in 1:(nrow(M.dom)-1)){
    for(j in (i+1):nrow(M.dom)){
      temp <- (M[i,] & M[j,])
      if(any(temp)) M.dom[i,j] <- M.dom[j,i] <- T
    }
  }
  
  # find the disconnected components that are the clusters
  D.member <- components.mat(M.dom)
  
  return(list(D.member=D.member,D=D,R=R))
}

# ccd clustering for arbitrary shapes of data
# edge correction included
# the balls are given with a K function estimation
# old name = ccd4.clustering_correct
# datax, the data set
# low.num, lowest number of a box cardinality 
# r.seq, the number of radii used for each window, proximity region
# method, the method for dominating set algorithm
# simul, the simulation results if provided 
rccd.clustering.connected_correct <- function(datax,low.num,r.seq,dom.method="greedy",simul=NULL){
  
  # info
  nr <- nrow(datax)
  nc <- ncol(datax)
  ddatax <- as.matrix(dist(datax))
  mac.eps <- .Machine$double.eps
  
  # find all radii, and the dominating set 
  ccd.info <- ccd.Kest.edge(datax,ddatax,low.num,r.seq,simul)
  
  r <- ccd.info$R
  M <- matrix(as.integer(ddatax < r+mac.eps), length(r))
  diag(M) <- 1
  score <- rowSums(M)
  
  # domination method is provided, find the dominating set
  # if not, take all radii
  if(!is.null(dom.method)){
    if(dom.method=="greedy") D <- dominate.mat.greedy(M)
    if(dom.method=="greedy2") D <- dominate.mat.greedy2(M)
    if(dom.method=="ks") D <- dominate.mat.ks(M,ks)
  } else {
    D <- 1:nr
  }
  R <- r[D]
  
  # find the intersection graph and its dominating set
  # the intersection graph is based on whether two radii catched the same points
  M <- M[D,]
  M.dom <- diag(T,nrow(M))
  for(i in 1:(nrow(M.dom)-1)){
    for(j in (i+1):nrow(M.dom)){
      temp <- (M[i,] & M[j,])
      if(any(temp)) M.dom[i,j] <- M.dom[j,i] <- T
    }
  }
  
  # find the disconnected components that are the clusters
  D.member <- components.mat(M.dom)
  
  return(list(D.member=D.member,D=D,R=R))
}

# the validation function for the clusters, this does not need the cluster labels
# old name = ccd4.clustering.nonvalid
# graph is the digraph CCD
# ddatax is the distance matrix
rccd.clustering.connected.nonvalid <- function(graph,ddatax){
  
  cls <- graph$D.member
  ncls <- length(unique(cls))
  
  D <- graph$D
  R <- graph$R

  ddx <- matrix(ddatax[,D],ncol=length(D))
  
  # which point is in which ball of dominators
  result <- t(apply(ddx,1,function(x){
    ind <- which.min(x/R)
    return(cls[ind])
  }))
  return(result)
}

# the validation function for the clusters
# old name = ccd4.clustering.valid
# graph is the digraph CCD
# ddatax is the distance matrix
# cls is the actual classes 
rccd.clustering.connected.valid <- function(graph,ddatax,cls){
  
  ncls <- length(unique(cls))
  
  D <- graph$D
  R <- graph$R
  cl <- cls[D]
  
  ddx <- matrix(ddatax[,D],ncol=length(D))
  
  # which point is in which ball of dominators
  result <- t(apply(ddx,1,function(x){
    ind <- which.min(x/R)
    return(cl[ind])
  }))
  return(result)
}

# ccd function that the radii of all points with K function, no edge correction, fast version
# low.num is the lowest cardinality of a ball
# r.seq is the number of breaks in the Kest analysis of radii
# simul is the provided simulation for the problem, if null, compute
ccd.Kest <- function(dx,ddx,low.num,r.seq,simul=NULL){
  
  n <- nrow(dx)
  d <- ncol(dx)
  R <- rep(0,n)
  rr <- seq(1/r.seq,1,1/r.seq)
  
  # the new simulation function
  if(!is.null(simul)) Kest.slopes <- simul 
  else Kest.slopes <- Kest.simpois(n,d,r.seq,99)

  for(i in 1:n){
    
    o.d <- order(ddx[i,])
    
    for(j in low.num:n){
      
      # get the K(r) of the data set
      r <- ddx[i,o.d[j]]*rr
      # Kest.obs <- Kest.f(ddx[o.d[1:j],o.d[1:j]],r,d)
      Kest.obs <- Kest.f(ddx[o.d[1:j],o.d[1:j]],r)
      
      # check the values, if rejected, set the R[i] as the radius
      lo.hi <- Kest.slopes$max[j,]
      flag <- (Kest.obs > lo.hi)
      if(any(flag)){
          if(j==low.num) R[i] <- 0
          else R[i] <- ddx[i,o.d[j-1]]
          break
        } else if(j==n){
          R[i] <- 0
          break
      }
    }
  }
  return(list(R=R,KS=NULL))
}

# ccd function that the radii of all points with K function, translation correction as the edge correction, fast version
# old name = ccd.Kest2
# low.num is the lowest cardinality of a ball
# r.seq is the number of breaks in the Kest analysis of radii
# simul is the provided simulation for the problem, if null, compute
ccd.Kest.edge <- function(dx,ddx,low.num,r.seq,simul=NULL){
  
  n <- nrow(dx)
  d <- ncol(dx)
  R <- rep(0,n)
  rr <- seq(1/r.seq,1,1/r.seq)
  
  # the new simulation function
  if(!is.null(simul)) Kest.slopes <- simul 
  else Kest.slopes <- Kest.simpois.edge(n,d,r.seq,99)
  
  # main iteration
  if(requireNamespace("progress"))
    pb <- progress_bar$new(total = n, clear = FALSE)
  for(i in 1:n){
    
    o.d <- order(ddx[i,])
    for(j in low.num:n){

      # get the K(r) of the data set
      r <- ddx[i,o.d[j]]*rr
      Kest.obs <- Kest.f.edge(ddx[o.d[1:j],o.d[1:j]],r,d)
      
      # check the values, if rejected, set the R[i] as the radius
      lo.hi <- Kest.slopes$max[j,]
      flag <- (Kest.obs > lo.hi)
      if(any(flag)){
        if(j==low.num) R[i] <- 0
        else R[i] <- ddx[i,o.d[j-1]]
        break
      } else if(j==n){
        R[i] <- 0
        break
      }
    }
    pb$tick()
  }
  return(list(R=R,KS=NULL))
}

# ccd function that the radii of all points with K function, translation correction as the edge correction, fast version
# incorporates quantiles from the simulated envelopes
# old name = ccd.Kest2.quantile
# low.num is the lowest cardinality of a ball
# r.seq is the number of breaks in the Kest analysis of radii
# quantile to be used as the confidence interval max
# simul is the provided simulation for the problem, if null, compute
ccd.Kest.edge.quantile <- function(dx,ddx,low.num,r.seq, quan, simul=NULL){
  
  n <- nrow(dx)
  d <- ncol(dx)
  R <- rep(0,n)
  rr <- seq(1/r.seq,1,1/r.seq)
  
  # the new simulation function
  if(!is.null(simul)) Kest.slopes <- simul 
  else Kest.slopes <- Kest.simpois.edge.quantile(n,d,r.seq,99)
  
  for(i in 1:n){
    
    o.d <- order(ddx[i,])
    for(j in low.num:n){
      
      # get the K(r) of the data set
      r <- ddx[i,o.d[j]]*rr
      Kest.obs <- Kest.f.edge(ddx[o.d[1:j],o.d[1:j]],r,d)
      
      # check the values, if rejected, set the R[i] as the radius
      lo.hi <- Kest.slopes$quan[[quan]][j,]
      flag <- (Kest.obs > lo.hi)
      if(any(flag)){
        if(j==low.num) R[i] <- 0
        else R[i] <- ddx[i,o.d[j-1]]
        break
      } else if(j==n){
        R[i] <- 0
        break
      }
    }
  }
  return(list(R=R,KS=NULL))
}

# ccd function that the radii of all points with K function, translation correction as the edge correction, fast version
# r is up to 0.5 length
# old name = ccd.Kest3
# low.num is the lowest cardinality of a ball
# r.seq is the number of breaks in the Kest analysis of radii
# simul is the provided simulation for the problem, if null, compute
ccd.Kest.edge.short <- function(dx,ddx,low.num,r.seq,simul=NULL){
  
  n <- nrow(dx)
  d <- ncol(dx)
  R <- rep(0,n)
  rr <- seq(1/r.seq,1,1/r.seq)
  rr <- rr[rr <= 0.5]
  
  # the new simulation function
  if(!is.null(simul)) Kest.slopes <- simul 
  else Kest.slopes <- Kest.simpois.edge(n,d,r.seq,99)

  for(i in 1:n){
    
    o.d <- order(ddx[i,])
    
    for(j in low.num:n){
      
      # get the K(r) of the data set
      r <- ddx[i,o.d[j]]*rr
      Kest.obs <- Kest.f.edge(ddx[o.d[1:j],o.d[1:j]],r,d)
      
      # check the values, if rejected, set the R[i] as the radius
      lo.hi <- Kest.slopes$max[j,1:length(rr)]
      flag <- (Kest.obs > lo.hi)
      if(any(flag)){
        if(j==low.num) R[i] <- 0
        else R[i] <- ddx[i,o.d[j-1]]
        break
      } else if(j==n){
        R[i] <- 0
        break
      }
    }
  }
  return(list(R=R,KS=NULL))
}

# ccd function that the radii of all points with K function, translation correction as the edge correction, fast version
# old name = ccd.Kest2
# low.num is the lowest cardinality of a ball
# r.seq is the number of breaks in the Kest analysis of radii
# simul is the provided simulation for the problem, if null, compute
ccd.Kest.edge_fast <- function(dx,ddx,low.num,r.seq,simul=NULL){
  
  n <- nrow(dx)
  d <- ncol(dx)
  R <- rep(0,n)
  rr <- seq(1/r.seq,1,1/r.seq)
  
  # the new simulation function
  if(!is.null(simul)) Kest.slopes <- simul 
  else Kest.slopes <- Kest.simpois.edge(n,d,r.seq,99)
  
  # main iteration
  if(requireNamespace("progress"))
    pb <- progress_bar$new(total = n, clear = FALSE)
  for(i in 1:n){
    
    # set search params
    o.d <- order(ddx[i,])
    search.flag <- FALSE
    search.list <- list()
    j <- floor((low.num+n)/2)
    search.list[[paste0(n)]] <- "reject"
    search.list[[paste0(low.num)]] <- "reject"
    
    # binary search
    while(!search.flag){

      # get the K(r) of the data set
      r <- ddx[i,o.d[j]]*rr
      Kest.obs <- Kest.f.edge(ddx[o.d[1:j],o.d[1:j]],r,d)
      
      # check the values
      # if rejected, search previous values
      # if accepted, search further values
      lo.hi <- Kest.slopes$max[j,]
      flag <- (Kest.obs > lo.hi)
      limits <- find_closest_search(j, as.numeric(names(search.list)))
      if(any(flag)){
        newj <- floor((j+limits$min)/2)
      } else {
        newj <- floor((j+limits$max)/2)
      }
      
      # get largest accepted radii if search is done
      if(j %in% as.numeric(names(search.list))){
        flags <- unlist(search.list)
        ind <- as.numeric(names(search.list))
        ind <- ind[flags == "accept"]
        if(length(ind) == (length(search.list)-1)) {
          R[i] <- 0
        } else if(length(ind) > 0){
          ind <- max(ind)
          R[i] <- ddx[i,o.d[ind]] 
        } else {
          R[i] <- 0
        }
        search.flag <- TRUE
      } else {
        if(any(flag)){
          search.list[[paste0(j)]] <- "reject"
        } else {
          search.list[[paste0(j)]] <- "accept"
        }
        j <- newj
      }
      
      # progress
    }
    pb$tick()
  }
  return(list(R=R,KS=NULL))
}

find_closest_search <- function(j, search.list){
  tmp <- search.list - j
  search.list.limits <- search.list[order(tmp, decreasing = FALSE)][1:2]
  return(list(min = search.list.limits[1], max = search.list.limits[2]))
}