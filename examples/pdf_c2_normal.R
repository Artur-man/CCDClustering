# this code is for the simulation of trying CCD codes

# library
library(MASS)
library(cluster)
library(pdfCluster)

# the main code
set.seed(1)
n <- 200
c1.mu <- c(0,0)
c2.mu <- c(5,0)
sigma <- diag(2)
datax <- mvrnorm(n,c1.mu,sigma)
datay <- mvrnorm(n,c2.mu,sigma)
dataf <- rbind(datax,datay)
ddataf <- as.matrix(dist(dataf))
ddataf <- ddataf <- as.matrix(dist(dataf))

# pdfcluster Clustering
h <- 0.5
clus.data <- pdfCluster(dataf, h = h)
clus <- groups(clus.data)

# print the clusters, use 1200 500
par(mar=c(0,0,0,0))
plot(dataf,xlab="",ylim=c(-2.3,1.75),ylab="",pch=16,cex=1.5,axes=FALSE, col = clus)
box()