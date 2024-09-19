# auxiliary functions for general use

draw.circle <- function (x, y, radius, nv = 100, border = NULL, col = NA, lty = 1, 
                         lwd = 1) 
{
  xylim <- par("usr")
  plotdim <- par("pin")
  angle.inc <- 2 * pi/nv
  angles <- seq(0, 2 * pi - angle.inc, by = angle.inc)
  if (length(col) < length(radius)) 
    col <- rep(col, length.out = length(radius))
  for (circle in 1:length(radius)) {
    xv <- cos(angles) * radius[circle] + x
    yv <- sin(angles) * radius[circle] + y
    polygon(xv, yv, border = border, col = col[circle], lty = lty, 
            lwd = lwd)
  }
  invisible(list(x = xv, y = yv))
}

# multivariate AUC measure for error rate
Mauc <- function(predict,actual){
  if(length(unique(actual)) == 1) return(0.5)
  
  rate <- NULL
  for(i in unique(actual)){
    ind <- which(actual==i)
    rate <- cbind(rate,sum(actual[ind]==predict[ind])/length(ind))
  }
  result <- mean(rate)
  return(result)
}

# multivariate AUC measure for error rate
# adjusted for different possible clusters
# try each permutation to see if which fits
Mauc.clus <- function(predict,actual){
  
  # get all permutations for the partitioning
  uni.clus <- sort(unique(predict))
  perm.label <- rbind(uni.clus,allPerms(uni.clus))
  
  # for all permutations
  max.auc <- 0
  for(j in 1:nrow(perm.label)){
    perm.predict <- perm.label[j,]
    perm.predict <- perm.predict[predict]
    temp.auc <- Mauc(perm.predict,actual)
    if(temp.auc > max.auc) max.auc <- temp.auc
  }
  result <- max.auc
  return(result)
}


# function for providing all possible permutations of a set
allPerms <- function(x){
  
  if(length(x) > 1){
    result <- NULL
    for(i in 1:length(x)){
      result <- rbind(result,cbind(x[i],allPerms(x[-i])))
    }
  } else{
    result <- x
  }
  return(result)
}

# gg_circle
gg_circle <- function(r, xc, yc, color="black", fill=NA, ...) {
  x <- xc + r*cos(seq(0, pi, length.out=100))
  ymax <- yc + r*sin(seq(0, pi, length.out=100))
  ymin <- yc + r*sin(seq(0, -pi, length.out=100))
  annotate("ribbon", x=x, ymin=ymin, ymax=ymax, color=color, fill=fill, ...)
}
