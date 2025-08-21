# this code is for the simulation of trying CCD codes

# library
library(MASS)
library(cluster)
library(progress)
library(fossil)

# source codes
source("../scripts/ccd_ripley.R")
source("../scripts/functions.R")
source("../scripts/ccdfunctions.R")
source("../scripts/Kest.R")

# the main code
set.seed(1)
n <- 200
c1.mu <- c(0,0)
c2.mu <- c(3,0)
clsdata <- rep(1:2, c(n,n))
inter <- 1
datax <- matrix(c(runif(n,c1.mu[1]-inter,c1.mu[1]+inter),
                  runif(n,c1.mu[2]-inter,c1.mu[2]+inter)),
                byrow=F,nrow=n)
datay <- matrix(c(runif(n,c2.mu[1]-inter,c2.mu[1]+inter),
                  runif(n,c2.mu[2]-inter,c2.mu[2]+inter)),
                byrow=F,nrow=n)
dataf <- rbind(datax,datay)
ddataf <- ddataf <- as.matrix(dist(dataf))

# the ccd.clustering code
simul <- Kest.simpois.edge(n,2,50,99)
ccd.sim <- rccd.clustering_correct(dataf, low.num=10, r.seq=50, dom.method = "greedy2", simul)

# find the best partitioning
ccd.si <- rccd.silhouette(ccd.sim,ddataf)
lind <- ccd.si$si.ind
ccd.sim$Int.D <- ccd.sim$Int.D[1:lind]
ccd.sim$Int.R <- ccd.sim$Int.R[1:lind]

# print the clusters, use 1200 500
par(mar=c(0,0,0,0))
plot(dataf,xlab="",ylim=c(-2.3,1.75),ylab="",pch=16,cex=1.5,axes=FALSE)
box()
D <- ccd.sim$Int.D
R <- ccd.sim$Int.R
for(i in 1:length(D))
  draw.circle(dataf[D[i],1],dataf[D[i],2],R[i],lwd=2)

# calculate rand index
result <- rccd.clustering.nonvalid(ccd.sim,ddataf)
print(adj.rand.index(result, clsdata))
print(rand.index(result, clsdata))
