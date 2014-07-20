
my.icp.2d.v2 <- function(reference, target, maxIter=10, minIter=5, pSample=0.5, threshold=0, isOpt=TRUE){
  
  #gets the amount of points, a.k.a. the domain
  m <- length(reference)
  
  if(threshold == 0)
    threshold <- round(m/3)
  
  #checks whether there are at least 2 non-zero points in both curves
  if(length(which(reference != 0)) < 2 && length(which(target != 0)) < 2)
    return(list(target = target, error = m, energyTotal = m, energyMean = m))
  
  #converts them into 2D matrices
  reference <- matrix(c(1:m, reference), nrow=m)
  target <- matrix(c(1:m, target), nrow=m)
  
  #takes out the null parts of both curves
  reference <- takeNoneFaceOut(reference)
  target <- takeNoneFaceOut(target)
  
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
    
    if(isOpt){
      
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
    target[,1] <- target[,1] + translationFactorX
    #performs the translation in Y
    target[,2] <- target[,2] + translationFactorY
    
    #measures the distances for each point
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
  (list(target = primeTarget, error = primeError, energyTotal = energy, energyMean = (energy/(i - 1))))
}

# Performs outlier correction on a curve
# input:
#   curve = a vector of numbers
curveCorrection3 <- function(curve, meanCurve, smooth=0, progress=FALSE){
  
  n <- length(curve)
  
  #checks whether there are at least 2 non-zero points in both curves
  if(length(which(curve != 0)) < 2)
    return(curve)
  
  imgVarX <- getXvar(matrix(c(1:n, curve), ncol=2))
  
  diff <- meanCurve[imgVarX$min:imgVarX$max] - curve[imgVarX$min:imgVarX$max]
  diff <- abs(diff - mean(diff))
  
  threshold <- 3*mean(diff)
  
  xout <- which(diff >= threshold)
  
  xs <- (imgVarX$min:imgVarX$max)[-xout]
  ys <- curve[xs]
  
  xout <- xout + imgVarX$min - 1
  
  if(length(which(xout == imgVarX$min)) > 0){
    #xout <- xout[-which(xout + imgVarX$min - 1 == imgVarX$min)]
    xs <- c(xout[which(xout == imgVarX$min)] - 1, xs)
    ys <- c(meanCurve[xout[which(xout == imgVarX$min)]], ys)
  }
  
  if(length(which(xout == imgVarX$max)) > 0){
    #xout <- xout[-which(xout + imgVarX$min - 1 == imgVarX$max)]
    xs <- c(xs, xout[which(xout == imgVarX$max)] + 1)
    ys <- c(ys, meanCurve[xout[which(xout == imgVarX$max)]])
  }
  
  if(length(xout) > 1){
  
    newCurve <- stinterp(xs, ys, xout)
    curve[newCurve$x] <- newCurve$y
  }
  
  #cat("Executed\n")
  if(smooth > 0)
    curve <- gaussianSmooth(curve, c(smooth))
  
  (curve)
}

# Performs outlier correction on a curve
# input:
#   curve = a vector of numbers
curveCorrection <- function(curve, progress=FALSE){
  
  n <- length(curve)
  
  imgVarX <- getXvar(matrix(c(1:n, curve), ncol=2))
  
  curvature <- curvatureVector(curve)
  on <- onset(curvature[which(curvature != 0)] + max(curvature))
  #on <- onset(curve)
  on <- abs(on - mean(on))
  c2 <- curvatureVector(on)
  above <- which(c2 > mean(c2))  
#  above <- which(c2 > max(c2)/4)

  if(length(above) == 0)
    return(curve)
  
  tb <- getTopBottom2(on[above], above)
  
  i <- 1
  while(i <= length(tb)){
        
    if(i == length(tb)){
      tb <- tb[-i]
    }
    else
    if(above[tb[[i+1]][1]] - above[tb[[i]][2]] <= 5){
      
      tb[[i]][2] <- tb[[i+1]][2]
      tb <- tb[-(i+1)]
      i <- i + 1
    }
    else{
      
      tb <- tb[-i]
    }
  }
  
  if(length(tb) > 0){
    
    xout <- list()
    index <- 1
    
    for(interval in tb){
      
      for(i in above[interval[1]]:above[interval[2]]){
        
        xout[[index]] <- i + imgVarX$min - 1
        index <- index + 1
      }
    }
    xout <- list2vector(xout)
    
    if(length(which(xout - imgVarX$min + 1 == 1)) > 0)
      xout <- xout[-which(xout - imgVarX$min + 1 == 1)]
    
    xs <- (imgVarX$min:imgVarX$max)[-(xout - imgVarX$min + 1)]
    ys <- curve[xs]
    
    newCurve <- stinterp(xs, ys, xout)
    curve[newCurve$x] <- newCurve$y
    
    curve <- curveCorrection(curve)
  }
  
  #cat("Executed\n")
  
  (curve)
}

# Performs outlier correction on a curve
# input:
#   curve = a vector of numbers
curveCorrection2 <- function(curve, smooth=0, progress=FALSE){
  
  n <- length(curve)
  
  imgVarX <- getXvar(matrix(c(1:n, curve), ncol=2))
  
  on <- onset(curve[imgVarX$min:imgVarX$max])
#  on <- onset(curvatureVector(curve))
  on <- abs(on - mean(on))
  c2 <- curvatureVector(on)
  above <- which(c2 > max(c2)/3)
  above <- above + imgVarX$min -1
  
  if(length(above) > 0){
      
    if(length(which(above == imgVarX$min)) > 0)
      above <- above[-which(above == imgVarX$min)]

    if(length(which(above == imgVarX$max)) > 0)
      above <- above[-which(above == imgVarX$max)]
    
    xs <- (imgVarX$min:imgVarX$max)[-(above - imgVarX$min + 1)]
    ys <- curve[xs]
    
    if(length(xs) > 1){
      
      newCurve <- stinterp(xs, ys, above)
      curve[newCurve$x] <- newCurve$y
    }
    
    if(smooth > 0)
      curve <- gaussianSmooth(curve, c(smooth))
  }
  
  #cat("Executed\n")
  
  (curve)
}

onset <- function(data){
  
  n <- length(data)
  
  on <- rep(0, n-1)
  
  for(i in 2:n){
    
    if(data[i-1] == 0)
      data[i-1] <- min(data[which(data != 0)])
    
    on[i-1] <- data[i]/data[i-1]
  }
  
  (on)
}

getTopBottom <- function(data){
  
  n <- length(data)
  
  tb <- list()
  
  s <- "-"
  if(data[2] - data[1] >= 0)
    s <- "+"
  
  top <- 0
  index <- 1
  bottom <- 0
  if(s == "-")
    top <- 1
  
  for(i in 3:n){
    
    if(data[i] - data[i-1] >= 0){
      if(s == "-")
        bottom <- i-1
      
      s <- "+"
    }
    else
    if(data[i] - data[i-1] < 0){
      if(s == "+")
        top <- i-1
      
      s <- "-"
    }
    
    if(top != 0 && bottom != 0){
      
      tb[[index]] <- c(top, bottom)
      index <- index + 1
      top <- 0
      bottom <- 0
    }
  }
  
  (tb)
}

getTopBottom2 <- function(data, above){
  
  n <- length(data)
  
  tb <- list()
  
  s <- "-"
  if(data[2] - data[1] >= 0)
    s <- "+"
  
  top <- 0
  index <- 1
  bottom <- 0
  if(s == "-")
    top <- 1
  
  for(i in 3:n){
    
    if(data[i] - data[i-1] >= 0){
      s <- "+"
    }
    else
      if(data[i] - data[i-1] < 0){
        if(s == "+")
          if(top == 0)
            top <- i-1
          else
            if(above[i-1] - above[top] <= 10)
              bottom <- i-1
            else
              top <- 0
        
        s <- "-"
      }
    
    if(top != 0 && bottom != 0){
      
      tb[[index]] <- c(top, bottom)
      index <- index + 1
      top <- 0
      bottom <- 0
    }
  }
  
  (tb)
}

getXvar <- function(m){
  
  variation <- getXvariation(m)
  
  while(m[variation$min,2] == 0)
    variation$min <- variation$min + 1
  
  while(m[variation$max, 2] == 0)
    variation$max <- variation$max - 1
  
  (list(max=variation$max, min=variation$min))
}
