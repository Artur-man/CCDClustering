# this code is for the simulation of trying CCD codes

# library
library(MASS)
library(cluster)
library(apcluster)
library(fossil)

# source codes
source("../scripts/apres.R")

# the main code
set.seed(1)
n <- 200
c1.mu <- c(0,0)
c2.mu <- c(5,0)
clsdata <- rep(1:2, c(n,n))
sigma <- diag(2)
datax <- mvrnorm(n,c1.mu,sigma)
datay <- mvrnorm(n,c2.mu,sigma)
dataf <- rbind(datax,datay)
ddataf <- as.matrix(dist(dataf))
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

# calculate rand index
print(adj.rand.index(apres_clusters, clsdata))
print(rand.index(apres_clusters, clsdata))

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

# calculate rand index
print(adj.rand.index(apres_clusters, clsdata))
print(rand.index(apres_clusters, clsdata))