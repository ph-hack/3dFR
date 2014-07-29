my.icp.2d.v2 <- function(reference, target, maxIter=10, minIter=5, pSample=0.5, threshold=0, isOpt=TRUE, isOpt2=TRUE){
  
  #gets the amount of points, a.k.a. the domain
  m <- length(reference)
  
  if(threshold == 0)
    threshold <- round(m/3)
  
  #checks whether there are at least 2 non-zero points in both curves
  if(!isOpt && (length(which(reference != 0)) < 2 || length(which(target != 0)) < 2))
    return(list(target = target, error = m, energyTotal = m, energyMean = m))
  
  #converts them into 2D matrices
  reference <- matrix(c(1:m, reference), nrow=m)
  target <- matrix(c(1:m, target), nrow=m)
  
  #takes out the null parts of both curves
  reference <- takeNoneFaceOut(reference)
  target <- takeNoneFaceOut(target)
  
  lr <- length(reference)
  lt <- length(target)
  if(isOpt && (lr < 5 || lt < 5)){
    
    if(lr < 3 || lt < 3)
      return(list(target = target, error = m, energyTotal = m, energyMean = m))
    
    if((reference[1,2] == 0 && reference[2,2] == 0) || (target[1,2] == 0 && target[2,2] == 0))
      return(list(target = target, error = m, energyTotal = m, energyMean = m))
  }
  
  if(commonDomain(reference, target, isOpt) == 0)
    return(list(target = target, error = m, energyTotal = m, energyMean = m))
  
  #remembers the prime target
  primeTarget <- target
  
  #computes the descriptor lines of both curves
  referenceLine <- linearInterpolation(reference, isOpt)
  targetLine <- linearInterpolation(target, isOpt)
  
  #computes the angle between them
  angle <- atan(referenceLine$a - targetLine$a)
  #performs the rotation in the target to make it closer to the reference
  target <- rotateCurve(target, 0, angle, isOpt)
  
  #interpolates the points in order to obtain interger coordinates in X
  target <- interpolateXinteger(target, isOpt)
  
  #checks whether the curves got too far
  if(commonDomain(reference, target, isOpt) >= threshold){
    
    #if they didn't, measures the distances for each point
    distances <- dist.p2p(reference, target, isOpt)
    #cat("they are common enougth\n")
    #cat("common domain = ", commonDomain(reference, target), "; threshold = ", threshold, "\n")
    #remembers the prime target
    primeTarget <- target
  }
  else{
    #otherwise...
    #retrieves the prime target
    target <- primeTarget
    #measures the distances for each point
    distances <- dist.p2p(reference, target, isOpt)
    #computes the mean error
    error <- mean(abs(distances))
    
    return(list(target = primeTarget, error = error, energyTotal = error, energyMean = error))
  }
  
  #computes the mean error
  error <- mean(abs(distances))
  
  #initializes the prime error that will always be equal or less than error
  primeError <- error
  primeDistances <- distances
  #initializes the iteration index with 1
  i <- 1
  #initializes the energy with 0
  energy <- 0
  
  #as long as the error keeps decreasing and the the maximum number of
  #iterations hasn't been reached ...
  while((error < primeError || i <= minIter) && i <= maxIter){
    
    #only if the error is less than the primeError...
    if(error < primeError){
      #remembers the prime error
      primeError <- error
      #remembers the prime target, that will always have the least error
      primeTarget <- target
      primeDistances <- distances
    }
    #sums the error into the energy
    energy <- energy + primeError
    
    m <- length(target[,1])
    nSamples <- round(m * pSample)
    dX <- round(1/pSample)
    
    samples <- 0
    if(m %% 2 == 0)
      samples <- 0:(nSamples - 1) * dX + sample(1:dX, 1)
    else
      samples <- 0:(nSamples - 1) * dX + sample(1:(dX-1), 1)
    
    translationFactorX <- 0
    translationFactorY <- 0
    
    if(FALSE){
      
      data <- list(target = target, reference = reference, translationFactorX = 0, translationFactorY = 0)
      
      data <- Reduce(function(data, sample){
        
        dists <- as.matrix(dist(rbind(target[sample,], reference)))[1,-1]
        k <- which.min(dists)
        
        data$translationFactorX <- data$translationFactorX + reference[k,1] - target[sample,1]
        data$translationFactorY <- data$translationFactorY + reference[k,2] - target[sample,2]
        
        return (data)
      }, samples, data)
      
      translationFactorX <- data$translationFactorX
      translationFactorY <- data$translationFactorY
    }
    else
      for(j in samples){
        
        dists <- as.matrix(dist(rbind(target[j,], reference)))[1,-1]
        k <- which.min(dists)
        
        translationFactorX <- translationFactorX + reference[k,1] - target[j,1]
        translationFactorY <- translationFactorY + reference[k,2] - target[j,2]
      }
    
    translationFactorX <- translationFactorX/nSamples
    translationFactorY <- translationFactorY/nSamples
    
    #performs the translation in X
    if(isOpt2)
      target[,1] <- target[,1] + round(translationFactorX)
    else
      target[,1] <- target[,1] + translationFactorX
    #performs the translation in Y
    target[,2] <- target[,2] + translationFactorY
    
    #measures the distances for each point
    if(!isOpt2)
      target <- interpolateXinteger(target, isOpt)
    distances <- dist.p2p(reference, target, isOpt)
    
    #checks whether the curves got too far
    if(commonDomain(reference, target, isOpt) >= threshold)
      #if they didn't, measures the error
      error <- mean(abs(distances))
    else
      #otherwise, sets the erro to the prime error plus 1
      error <- primeError + 1
    
    #cat("Iteration ", i, "; error = ", error, "\n")
    #increasing the iteration index
    i <- i + 1
  }
  #returns the informations
  (list(target = primeTarget, dist = distances, error = primeError, energyTotal = energy, energyMean = (energy/(i - 1))))
}

# Compares 2 lines (curves) by better matching them throught linear transformations
# and returns the mean and last error. The error is measured point-to-point.
# input:
#   reference = a vector of numbers
#   target = a vector of numbers
#   maxIter = the maximum number of iterations allowed, by default 10
#   threshold = the minimum amount of points to be considered. By default
#               it is equal to third part of the number of points of the reference
#               If there is no such amount of points available for measurement, the
#               error measurement will be returned
# output:
#   a list containing:
#     'target' = the target line after the applied transformations
#     'error' = the last error computed
#     'energyTotal' = the sum of the errors of all iterations
#     'energyMean' = the mean of the errors of all iterations
my.icp.2d <- function(reference, target, by="mean", maxIter=10, threshold=0){
  
  #gets the amount of points, a.k.a. the domain
  m <- length(reference)
  
  if(threshold == 0)
    threshold <- m/3
  
  #checks whether there are at least 2 non-zero points in both curves
  if(length(which(reference != 0)) < 2 && length(which(target != 0)) < 2)
    return(list(target = target, error = m, energyTotal = m, energyMean = m))
  
  #converts them into 2D matrices
  reference <- matrix(c(1:m, reference), nrow=m)
  target <- matrix(c(1:m, target), nrow=m)
  
  #takes out the null parts of both curves
  reference <- takeNoneFaceOut(reference)
  target <- takeNoneFaceOut(target)
  
  if(commonDomain(reference, target) == 0)
    return(list(target = target, error = m, energyTotal = m, energyMean = m))
  
  #remembers the prime target
  primeTarget <- target
  
  #computes the descriptor lines of both curves
  referenceLine <- linearInterpolation(reference)
  targetLine <- linearInterpolation(target)
  
  #computes the angle between them
  angle <- atan(referenceLine$a - targetLine$a)
  #performs the rotation in the target to make it closer to the reference
  target <- rotateCurve(target, 0, angle)
  
  #computes the distance
  distances <- dist.p2p(reference, target)
  
  #computes the initial translation parameters
  translationFactorX <- reference[which.max(reference[,2]),1] - target[which.max(target[,2]),1]
  translationFactorY <- mean(distances)
  
  #performs the translation
  target[,1] <- target[,1] + translationFactorX
  target[,2] <- target[,2] + translationFactorY
  
  #interpolates the points in order to obtain interger coordinates in X
  #target <- interpolateXinteger(target)
  
  #checks whether the curves got too far
  if(commonDomain(reference, target) >= threshold){
    
    #if they didn't, measures the distances for each point
    distances <- dist.p2p(reference, target)
    #cat("they are common enougth\n")
    #cat("common domain = ", commonDomain(reference, target), "; threshold = ", threshold, "\n")
  }
  else{
    #otherwise...
    #retrieves the prime target
    target <- primeTarget
    #measures the distances for each point
    distances <- dist.p2p(reference, target)
    #computes the mean error
    #error <- m
    #if(by == "mean")
    error <- mean(abs(distances))
    #else if(by == "cosine")
    #error <- cosineDist(reference[,2], target[,2])
    
    return(list(target = primeTarget, error = error, energyTotal = error, energyMean = error))
  }
  
  #computes the mean error
  error <- mean(abs(distances))
  
  #initializes the prime error bigger than the 1st computed error
  primeError <- error + 1
  #initializes the iteration index with 1
  i <- 1
  #initializes the energy with 0
  energy <- 0
  
  #as long as the error keeps decreasing and the the maximum number of
  #iterations hasn't been reached ...
  while(error < primeError && i <= maxIter){
    
    #remembers the prime error
    primeError <- error
    #remembers the prime target
    primeTarget <- target
    #sums the error into the energy
    energy <- energy + primeError
    
    #computes the scale factor Y
    factors <- factor.p2p(reference, target, distances)
    scaleFactorY <- (max(factors) + min(factors[which(factors != 0)])) /2
    
    #computes the scale factor X
    refXvar <- getXvariation(reference)
    tarXvar <- getXvariation(target)
    #attenpting to consider only a peace of the reference
    #     if(refXvar$min < tarXvar$max && refXvar$min > tarXvar$min && refXvar$max > tarXvar$max)
    #       tarXvar$min <- refXvar$min
    #     else if(refXvar$max < tarXvar$max && refXvar$max > tarXvar$min && refXvar$min < tarXvar$min)
    #       tarXvar$max <- refXvar$max
    scaleFactorX <- (refXvar$max - refXvar$min)/(tarXvar$max - tarXvar$min)
    
    #performs the scalling
    target[,2] <- target[,2] * (1 + scaleFactorY)
    target[,1] <- target[,1] * scaleFactorX
    
    #if the X coordinates changed, interpolates the points in order to obtain interger coordinates in X
    if(scaleFactorX != 1)
      target <- interpolateXinteger(target)
    
    #performs the translation in X
    translationFactorX <- reference[which.max(reference[,2]),1] - target[which.max(target[,2]),1]
    target[,1] <- target[,1] + translationFactorX
    
    #computes the distance
    distances <- dist.p2p(reference, target)
    #performs the translation in Y
    translationFactorY <- mean(distances)
    target[,2] <- target[,2] + translationFactorY
    
    #measures the distances for each point
    target <- interpolateXinteger(target)
    distances <- dist.p2p(reference, target)
    
    #checks whether the curves got too far
    if(commonDomain(reference, target) >= threshold)
      #if they didn't, measures the error
      error <- mean(abs(distances))
    else
      #otherwise, sets the erro to the prime error plus 1
      error <- primeError + 1
    
    #cat("Iteration ", i, "; error = ", error, "\n")
    #increasing the iteration index
    i <- i + 1
  }
  #returns the informations
  #if(by == "mean")
  #primeError <- mean(abs(distances))
  #else if(by == "cosine")
  #primeError <- cosineDist(reference[,2], target[,2])
  (list(target = primeTarget, error = primeError, energyTotal = energy, energyMean = (energy/(i - 1))))
}

# Computes the common domain between two lines/curves.
# This is given by the number of points whose domains belongs to
# both lines/curves.
# input:
#   reference = a 2D matrix, the reference line/curve
#   target = a 2D matrix, the target line/curve
# output:
#   a integer corresponding to the number of target points whose domain
#   (1st column) are common to both lines/curves
commonDomain <- function(reference, target, isOpt=FALSE){
  
  #gets the number of points of the target
  n <- length(target[,1])
  #initializes the result with 0
  k <- 0
  
  #if it is to use the optimal algorithm
  if(isOpt){
    
    cMin <- max(reference[1,1], target[1,1])
    cMax <- min(reference[(length(reference[,1])),1], target[n,1])
    k <- cMax - cMin + 1
  }
  else
    #for each point...
    for(i in 1:n){
      #if the there is at least one reference point with the same domain ...
      if(length(which(reference[,1] == target[i,1])) > 0)
        #adds 1 into 'k'
        k <- k + 1
    }
  
  #returns the result
  (k)
}

# Computes the distance in 'y' for each point by the target's coordinate
dist.p2p <- function(reference, target, isOpt=FALSE){
  
  n <- length(target[,1])
  
  distances <- rep(0, n)
  notPresent <- 0
  
  if(isOpt){
    
    m <- length(reference[,1])
    
    xmin <- max(target[1,1], reference[1,1])
    xmax <- min(reference[m,1], target[n,1])
    
    ref <- (Position(function(x){return(x == xmin)}, reference[,1]):Position(function(x){return(x == xmax)}, reference[,1]))
    tar <- (Position(function(x){return(x == xmin)}, target[,1]):Position(function(x){return(x == xmax)}, target[,1]))
    
    distances <- mapply(difference, reference[ref,2], target[tar,2])
  }
  else
    for(i in 1:n){
      
      x <- target[i,1]
      
      if(length(which(reference[,1] == x)) > 0)
        distances[i] <- reference[which(reference[,1] == x), 2] - target[i,2]
      else
        notPresent <- c(notPresent, i)
    }
  
  notPresent <- notPresent[-1]
  if(length(notPresent) > 0)
    distances <- distances[-notPresent]
  
  (distances)
}

# Removes the points with black/zero value.
# input:
#   data = a 2D matrix containing a line/curve
# output:
#   a 2D matrix containing the biggest part of 'data' which
#   doesn't contain black/zero value.
takeNoneFaceOut <- function(data){
  
  #finds the domain whose image is zero
  noface <- data[which(data[,2] == 0),1]
  #adds the first and the last points
  noface <- c(data[1,1], noface, data[length(data[,1]),1])
  #finds where is the first biggest interval
  biggestInterval <- which.max(differenceVector(noface))
  #finds the interval itself
  interval <- c(noface[biggestInterval]+1, noface[biggestInterval+1]-1)
  
  #removes all points but the found interval
  data <- data[interval[1]:interval[2],]
  
  #returns the cropped data
  (data)
}

# Rotates a given curve/line based on a given reference coordinate
# and a given angle. The transformation is perfomed by multiplying matrices.
# input:
#   curve = a 2D matrix of numbers with 2 columns, one for each dimension
#   referenceX = the coordinate of the 1st column dimension where to fix the rotation
#   angle = the angle of the rotation
# output:
#   a 2D matrix of numbers, corresponding the curve/line rotated
rotateCurve <- function(curve, referenceX, angle, isOpt=FALSE){
  
  #gets the number of points of the curve/line
  m <- length(curve[,1])
  
  #translates the curve so the referenceX point will be at the origin
  curve[,1] <- curve[,1] - referenceX
  
  #builds the rotation matrix
  rotationMatrix <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), ncol=2)
  
  if(isOpt){
    curve <- t(rotationMatrix %*% t(curve))
  }
  else
    #for each point...
    for(i in 1:m)
      #computes the rotation by multiplying the rotation matrix with the point
      curve[i,] <- rotationMatrix %*% matrix(curve[i,], nrow=2)
  
  #translates back the curve as the original one
  curve[,1] <- curve[,1] + referenceX
  
  #interpolates the curve to make sure all 1st column coordinates will be integers
  curve <- interpolateXinteger(curve)
  
  #returns the curve rotated
  (curve)
}

# Interpolates the signal m, a 2D matrix, in a way that all 'x's (1st column)
# will be integers.
# input:
#   m = a curve, represented as a 2D matrix with 2 columns, one for each dimension
# output:
#   a 2D matrix containing the curve with all 1st column values as integers
interpolateXinteger <- function(m, isOpt=FALSE){
  
  #gets the number of points (rows) of 'm'
  n <- length(m[,2])
  
  m2 <- 0
  if(isOpt){
    
    x <- floor(m[1,1])
    N <- floor(m[n,1]) - x
    
    if(AND(m[,1] != floor(m[,1])))
      
      #m2 <- t(mapply(discretizePoint, m[1:(n-1),1], m[1:(n-1),2], m[2:n,1], m[2:n,2]))
      m2 <- t(mapply(discretizePoint, x:(x + N - 1), MoreArgs=list(m = m)))
    #return(m2)
    else
      return (m)
  }
  else
    #for each point, but the last...
    for(i in 1:(n-1)){
      
      #computes the line which passes through the 'i'th point and its successor
      line <- findLine(m[i,], m[i+1,])
      
      #if the difference between the domain values of the 'i'th point and its successor
      #is equal or greater than a half, interpolates the 'i'th point
      if(abs(m[i,1] - m[i+1,1]) >= 0.5){
        
        #finds the closest integer
        x <- round(m[i,1])
        #finds the correspoding value (2nd column) to this new domain value
        y <- appLinear(line, x)
        #applies the results
        m[i,1] <- x
        m[i,2] <- y
      }
    }
  
  #does the same thing to the last point
  line <- findLine(m[n-1,], m[n,])
  
  if(abs(m[n-1,1] - m[n,1]) >= 0.5){
    
    x <- floor(m[n,1])
    y <- appLinear(line, x)
    m[n,1] <- x
    m[n,2] <- y
  }
  
  #returns the interpolated curve
  if(isOpt)
    m <- rbind(m2, m[n,])
  
  #removes the remaining floating domains
  #m <- m[which(m[,1] == ceiling(m[,1])),]
  
  (m)
}