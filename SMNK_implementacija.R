#to run all lines, run next brackets

{
  
########################
# LEAST SQUARES METHOD #
########################

sumX <- function(pX, pY, beg, end) {
  sum <- 0
  
  for (i in beg:end) {
    sum <- sum + pX[i]
  }
  
  return(sum)
}

sumY <- function(pX, pY, beg, end) {
  sum <- 0
  
  for (i in beg:end) {
    sum <- sum + pY[i]
  }
  
  return(sum)
}

sumXY <- function(pX, pY, beg, end) {
  sum <- 0
  
  for (i in beg:end) {
    sum <- sum + pX[i]*pY[i]
  }
  
  return(sum)
}

sumXX <- function(pX, pY, beg, end) {
  sum <- 0
  
  for (i in beg:end) {
    sum <- sum + pX[i]^2
  }
  
  return(sum)
}

slope <- function(pX, pY, beg, end) {
  
  S_x <- sumX(pX, pY, beg, end)
  S_y <- sumY(pX, pY, beg, end)
  S_xy <- sumXY(pX, pY, beg, end)
  S_xx <- sumXX(pX, pY, beg, end)
  
  n <- 1 + end - beg
  numerator <- n * S_xy - S_x * S_y
  denominator <- n * S_xx - S_x^2
  
  return(numerator/denominator)
}

intercept <- function(slope, pX, pY, beg, end) {
  S_x <- sumX(pX, pY, beg, end)
  S_y <- sumY(pX, pY, beg, end)
  
  n <- 1 + end - beg
  numerator <- S_y - slope * S_x
  return(numerator/n)
}

lineModel <- function(pX, pY, beg, end) {
  a <- slope(pX, pY, beg, end)
  b <- intercept(a, pX, pY, beg, end)
  return( c(b,a) )
}

#################
# THE ALGORITHM #
#################

error <- function(intercept, slope, pX, pY, beg, end) {
  sum <- 0
  
  for (i in beg:end) {
    sum <- sum + (pY[i] - pX[i]*slope - intercept)^2
  }
  
  return(sum)
}

minError <- function(pX, pY, beg, end) {
  if(beg==end) return(0)
  
  lm <- lineModel(pX, pY, beg, end)
  b <- lm[1]
  a <- lm[2]
  err <- error(b, a, pX, pY, beg, end)
  
  return(err)
}

findSegments <- function(j, C, OPT, E) {
  min <- .Machine$double.xmax #max floating point number
  index <- NULL
  
  for (i in 1:j) {
    value <- E[i, j] + C + OPT[i]
    if(value < min) {
      min <- value
      optIndex <- i
    }
  }
  
  s2 <- c(optIndex, j)
  if(optIndex==1) return(s2)
  
  s1 <- findSegments(optIndex-1, C, OPT, E)
  return( c(s1, s2) )
}

segmentedLeastSquares <- function(pX, pY, n, C) {
  # OPT[k] = value of optimal solution for p1, p2,..., pk-1
  # our goal is OPT[n+1]
  OPT <- c()
  OPT[1] = 0 
  E <- matrix(0, ncol=n, nrow=n)
  
  for (i in 1:n) {
    for (j in i:n) {
      if (j - i >= 1) E[i,j] <- minError(pX, pY, i, j)
    }
  }
  
  for (j in 1:n) {
    
    candidates <- c()
    for (i in 1:j) candidates[i] <- E[i,j] + C + OPT[i]
    
    OPT[j + 1] <- min(candidates)
  }
  
  partition <- findSegments(n, C, OPT, E)
  return(partition)
}

#################
# VISUALIZATION #
#################

drawLineModel <- function(pX, pY, beg, end) {
  lm <- lineModel(pX, pY, beg, end)
  b <- lm[1]
  a <- lm[2]
  
  x1 <- pX[beg]
  y1 <- a * x1 + b
  x2 <- pX[end]
  y2 <- a * x2 + b
  
  segments(x1,y1,x2,y2,col="blue")
}

drawSegments <- function(pX, pY, partition) {
  k <- length(partition)/2
  for (i in 1:k) {
    beg <- partition[2*i -1]
    end <- partition[2*i]
    drawLineModel(pX, pY, beg, end)
  }
}

########
# DATA #
########

xt <- c(1,2,3,5,6.1,7,8,9,10)
yt <- c(1,2.4,3,6,5.8,5.4,4,4,2)


########
# MAIN #
########

{
plot(xt,yt, pch=19, col="darkorange")
partition <- segmentedLeastSquares(xt,yt,9,0.1)
drawSegments(xt, yt, partition)
}

}
