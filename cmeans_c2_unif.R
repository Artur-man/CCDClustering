# this code is for the simulation of trying CCD codes

# library
library(MASS)
library(cluster)
library(e1071)

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

# Fuzzy C means clustering
k <- 2
result <- cmeans(dataf,k)
# result <- cmeans.cluster(clus.fuzzy$centers, dataf)

# print the clusters, use 1200 500
par(mar=c(0,0,0,0))
plot(dataf,xlab="",ylim=c(-2.3,1.75),ylab="",pch=16,cex=1.5,axes=FALSE, col = result$cluster)
box()
