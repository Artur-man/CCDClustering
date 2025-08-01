# this code is for the simulation of trying CCD codes

# library
library(MASS)
library(cluster)
library(apcluster)

# the main code
set.seed(1)
n <- 200
c1.mu <- c(0,0)
c2.mu <- c(3,0)
inter <- 1
datax <- matrix(c(runif(n,c1.mu[1]-inter,c1.mu[1]+inter),
                  runif(n,c1.mu[2]-inter,c1.mu[2]+inter)),
                byrow=F,nrow=n)
datay <- matrix(c(runif(n,c2.mu[1]-inter,c2.mu[1]+inter),
                  runif(n,c2.mu[2]-inter,c2.mu[2]+inter)),
                byrow=F,nrow=n)
dataf <- rbind(datax,datay)
ddataf <- ddataf <- as.matrix(dist(dataf))

# AP Clustering (q = 0)
apres <- apcluster(negDistMat(r=2), dataf, details=TRUE,q=0)
apres_ind <- unlist(apres@clusters)
apres_clusters <- rep(1:length(apres@clusters),lapply(apres@clusters,length))
names(apres_clusters) <- names(apres_ind)
apres_clusters <- apres_clusters[paste0(1:length(apres_clusters))]

# print the clusters, use 1200 500
par(mar=c(0,0,0,0))
plot(dataf,xlab="",ylim=c(-2.3,1.75),ylab="",pch=16,cex=1.5,axes=FALSE, col = apres_clusters)
box()

# AP Clustering (q = 0.5)
apres <- apcluster(negDistMat(r=2), dataf, details=TRUE,q=0.5)
apres_ind <- unlist(apres@clusters)
apres_clusters <- rep(1:length(apres@clusters),lapply(apres@clusters,length))
names(apres_clusters) <- names(apres_ind)
apres_clusters <- apres_clusters[paste0(1:length(apres_clusters))]

# print the clusters, use 1200 500
par(mar=c(0,0,0,0))
plot(dataf,xlab="",ylim=c(-2.3,1.75),ylab="",pch=16,cex=1.5,axes=FALSE, col = apres_clusters)
box()