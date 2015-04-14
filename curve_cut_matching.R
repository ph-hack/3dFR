my.dtwBasedDistance <- function(reference, target, smooth=0, W=0, range=0, threshold=0, tol=5){
  
  M <- length(reference)
  N <- length(target)
  
  if(threshold == 0)
    threshold <- round(M/3)

  if(smooth > 0){
    reference <- gaussianSmooth(c(rep(reference[1], 2*smooth), reference, rep(reference[M], 2*smooth)), c(smooth))[(1 + 2*smooth):(M + 2*smooth)]
    target <- gaussianSmooth(c(rep(target[1], 2*smooth), target, rep(target[N], 2*smooth)), c(smooth))[(1 + 2*smooth):(N + 2*smooth)]
  }
  
  #converts them into 2D matrices
  reference <- my.as.matrix(reference)
  target <- my.as.matrix(target)
  
  #takes out the null parts of both curves
  reference <- takeNoneFaceOut(reference, TRUE)
  target <- takeNoneFaceOut(target, TRUE)
  
  #computes the angle between them
  #angle <- angleBetween(reference, target)[1]
  #performs the rotation in the target to make it closer to the reference
  #target <- rotateCurve(target, 0, -angle, TRUE)
  #interpolates the points in order to obtain interger coordinates in X
  #target <- interpolateXinteger(target, TRUE)
  
#   M <- length(reference[,2])
#   N <- length(target[,2])
#   
#   dtw <- matrix(rep(100, (N+1)*(M+1)), nrow=N+1)
#   dtw[1,1] <- 0
  
  #gRef <- gradient(reference[,2], FALSE)
  gRef <- computeCurvature(reference, FALSE)
  #gTar <- gradient(target[,2], FALSE)
  gTar <- computeCurvature(target, FALSE)
  
  M <- length(gRef)
  N <- length(gTar)
  
  dtw <- matrix(rep(100, (N+1)*(M+1)), nrow=N+1)
  dtw[1,1] <- 0
  
  if(W == 0)
    W <- M
  else
    W <- max(W, abs(N-M))
  
  for(i in 2:(N+1)){
    
    for(j in (max(2, i-W)):(min(M+1, i+W))){
      
      d <- my.gradDistance(gRef, gTar, j-1, i-1, range)
      dmin <- min(dtw[i-1,j], dtw[i,j-1], dtw[i-1,j-1])
      
      #if(d <= dtw[i,j-1])
      #  dtw[i,j] <- d
      #else
        dtw[i,j] <- d #+ dmin
    }
  }

  dtw <- dtw[-1,-1]
  
  image(dtw, col = topo.colors(M*N), xlab = "Target", ylab = "Reference")
  
  startPtsTar <- mapply(function(x){
    
    if(x == 1 && dtw[1,x] < dtw[1,x+1])
      return(c(1,x))
    
    else if(x == M && dtw[1,x] < dtw[1,x-1])
      return(c(1,x))
    
    else if(dtw[1,x] < dtw[1,x-1] && dtw[1,x] < dtw[1,x+1])
      return(c(1,x))
    
  }, 1:M)

  startPtsRef <- mapply(function(x){
    
    if(x == 1 && dtw[x,1] < dtw[x+1,1])
      return(c(x,1))
    
    else if(x == N && dtw[x,1] < dtw[x-1,1])
      return(c(x,1))
    
    else if(dtw[x,1] < dtw[x-1,1] && dtw[x,1] < dtw[x+1,1])
      return(c(x,1))
    
  }, 1:N)

  startPtsRef <- removeNullFromList(startPtsRef)
  startPtsTar <- removeNullFromList(startPtsTar)

  startPoints <- concatenateList(list(startPtsTar, startPtsRef))

  changeStart <- TRUE
  i <- 1
  j <- 1
  pts <- list()
  Np <- 0
  tolerance <- 0
  sPoint <- 0

  while(changeStart){
    
    sPoint <- sPoint + 1
    cat("sPoint: ", sPoint, " : ", startPoints[[sPoint]], "\n")
    startPoint <- startPoints[[sPoint]]
    
    i <- startPoint[1]
    j <- startPoint[2]
    
    Np <- Np + 1
    pts[[Np]] <- c(i,j)
    i <- i + 1
    
    findPoints <- TRUE
    
    while(findPoints){
    
      if(i <= N){
        
        if(j < M){
          
          if(dtw[i,j+1] <= dtw[i,j] && dtw[i,j+1] <= dtw[i-1,j+1]){
            
            Np <- Np + 1
            j <- j + 1
            pts[[Np]] <- c(i,j)
            i <- i + 1
          }
          else if(dtw[i,j] < dtw[i,j+1] && dtw[i,j] <= dtw[i-1,j+1]){
            
            if(tolerance + 1 <= tol){
              
              Np <- Np + 1
              pts[[Np]] <- c(i,j)
              i <- i + 1
              tolerance <- tolerance + 1
            }
            else{
              
              findPoints <- FALSE
              if(sPoint == length(startPoints))
                changeStart <- FALSE
              Np <- 0
              tolerance <- 0
              pts <- list()
            }
          }
          else if(dtw[i-1,j+1] < dtw[i,j+1] && dtw[i-1,j+1] < dtw[i,j]){
            
            if(tolerance - 1 >= -tol){
              
              Np <- Np + 1
              j <- j + 1
              pts[[Np]] <- c(i-1,j)
              tolerance <- tolerance - 1
            }
            else{
              
              findPoints <- FALSE
              if(sPoint == length(startPoints))
                changeStart <- FALSE
              Np <- 0
              tolerance <- 0
              pts <- list()
            }
          }
        }
        else{
          
          if(Np < threshold){
            
            Np <- 0
            tolerance <- 0
            pts <- list()
            
            if(sPoint == length(startPoints))
              changeStart <- FALSE
          }
          else{
            
            changeStart <- FALSE
          }
          findPoints <- FALSE
        }
      }
      else{
        
        if(Np < threshold){
          
          Np <- 0
          tolerance <- 0
          pts <- list()
          
          if(sPoint == length(startPoints))
            changeStart <- FALSE
        }
        else{
          
          changeStart <- FALSE
        }
        findPoints <- FALSE
      }
    }
  }

  cat("Tolerance: ", tolerance, "\n")
  
  if(length(pts) > 0){
    pts <- list2matrix(pts)
    p <- pts
    p[,1] <- p[,1]/N
    p[,2] <- p[,2]/M
    lines(p, col="red")
    
    minPoints <- apply(dtw, 1, which.min)
    mp <- my.as.matrix(minPoints)
    mp[,1] <- mp[,1]/N
    mp[,2] <- mp[,2]/M
    #points(mp, col="red")
  
    lastPoint <- pts[length(pts[,1]),]
    
    return(list(dtw=dtw, distance=dtw[lastPoint[1], lastPoint[2]], points=pts, minPoints=minPoints))
  }
  else{
    
    cat("No match found!\n")
    return(list(dtw=dtw, distance=100, points=c(), minPoints=c()))
  }
}

my.gradDistance <- function(C1, C2, x1, x2, range=0){
  
  #g1 <- gradient(C1, FALSE)
  #g2 <- gradient(C2, FALSE)
  g1 <- C1
  g2 <- C2
  
  w1 <- getWindow(g1, x1, range)
  w2 <- getWindow(g2, x2, range)
  
  D <- abs(g1[x1] - g2[x2])
  d <- 0
  
  if(length(w1) > 0 && length(w2) > 0){
    
    weights <- mapply(function(x){
      
      return(x*1/(2*range))
    }, 1:range)
    weights <- c(weights[range:1], weights)
    
    d <- mapply(function(k,w,x,y,X,Y){
      
      if(is.na(X[x+w]) || is.na(Y[y+w]) || length(X[x+w]) != 1 || length(Y[y+w]) != 1)
        return(100)
      
      else{
        
        return(abs(k * (X[x+w] - Y[y+w])))
      }
    }, weights, c(-range:-1,1:range), MoreArgs = list(x=x1, y=x2, X=g1, Y=g2))
    
    d <- list2vector(d)
    
    noNeighbor <- which(d == 100)
    if(length(noNeighbor) > 0)
      d <- d[-noNeighbor]
    
  }
  return(D + sum(d))
}

computeCurvature <- function(data, absolute=TRUE){
  
  if(length(dim(data)) < 2)
    data <- my.as.matrix(data)
  
  N <- length(data[,2])
  
  C <- mapply(function(i, X, Y){
    
    a <- 0
    
    if(i == 1){
      
      a <- pi/2
      
      if(X[i+1] - X[i] != 0)
        a <- atan(angularCoeff(X[i], Y[i], X[i+1], Y[i+1]))
      
      if(Y[i] > Y[i+1])
        a <- pi/2 + a
    }
    else if(i == N){
      
      a <- pi/2
      
      if(X[i-1] - X[i] != 0)
        a <- atan(angularCoeff(X[i], Y[i], X[i-1], Y[i-1]))
      
      if(Y[i-1] > Y[i])
        a <- pi/2 + a
    }
    else{
      
      a1 <- pi/2
      
      if(Y[i] > Y[i+1])
        a1 <- -a1
      
      if(X[i+1] - X[i] != 0)
        a1 <- atan(angularCoeff(X[i], Y[i], X[i+1], Y[i+1]))
      
      a2 <- pi/2
      
      if(Y[i-1] > Y[i])
        a2 <- -a2
      
      if(X[i-1] - X[i] != 0)
        a2 <- atan(angularCoeff(X[i], Y[i], X[i-1], Y[i-1]))
      
      if(a1 >= 0 && a2 >= 0){
        
        a1 <- pi/2 - a1
        a <- a1 + a2 + pi/2
      }
      else if(a1 < 0 && a2 < 0){
        
        a2 <- pi/2 + a2
        a <- -a1 + a2 + pi/2
      }
      else if(a1 >= 0 && a2 < 0){
        
        a <- pi + a2 - a1
      }
      else{
        
        a <- pi + a2 - a1
      }
    }
    
    return(a)
    
  }, 2:(N-1), MoreArgs = list(X=data[,1], Y=data[,2]))
  
  if(absolute)
    return(abs(C))
  else
    return(C)
}

my.icp.2d.v2 <- function(reference, target, maxIter=10, minIter=5, pSample=0.5, threshold=0, isOpt=TRUE, isOpt2=TRUE){
  
  #gets the amount of points, a.k.a. the domain
  m <- length(reference)
  
  if(threshold == 0)
    threshold <- round(m/3)
  
  #checks whether there are at least 2 non-zero points in both curves
  if(!isOpt && (length(which(reference != 0)) < 2 || length(which(target != 0)) < 2))
    return(list(target = target, error = m, energyTotal = m, energyMean = m, dist = matrix(c(1:m, rep(0, m)), ncol=2)))
  
  #converts them into 2D matrices
  reference <- matrix(c(1:m, reference), nrow=m)
  target <- matrix(c(1:(length(target)), target), nrow=length(target))
  
  #takes out the null parts of both curves
  reference <- takeNoneFaceOut(reference, TRUE)
  target <- takeNoneFaceOut(target, TRUE)
  
  lr <- length(reference)
  lt <- length(target)
  if(isOpt && (lr < 5 || lt < 5)){
    
    if(lr < 3 || lt < 3)
      return(list(target = target, error = m, energyTotal = m, energyMean = m, dist = matrix(c(1:m, rep(0, m)), ncol=2)))
    
    if((reference[1,2] == 0 && reference[2,2] == 0) || (target[1,2] == 0 && target[2,2] == 0))
      return(list(target = target, error = m, energyTotal = m, energyMean = m, dist = matrix(c(1:m, rep(0, m)), ncol=2)))
  }
  
  if(commonDomain(reference, target, isOpt) == 0)
    return(list(target = target, error = m, energyTotal = m, energyMean = m, dist = matrix(c(1:m, rep(0, m)), ncol=2)))
  
  #remembers the prime target
  primeTarget <- target
  
  #computes the descriptor lines of both curves
  #referenceLine <- linearInterpolation(reference, isOpt)
  #targetLine <- linearInterpolation(target, isOpt)
  
  #computes the angle between them
  #angle <- atan(referenceLine$a - targetLine$a)
  angle <- angleBetween(reference, target)[1]
  #performs the rotation in the target to make it closer to the reference
  target <- rotateCurve(target, 0, -angle, isOpt)
  
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
    error <- mean(abs(distances[,2]))
    
    return(list(target = primeTarget, error = error, energyTotal = error, energyMean = error, dist = distances))
  }
  
  #computes the mean error
  #error <- mean(abs(distances[,2]))
  error <- 100
  
  #initializes the prime error that will always be equal or less than error
  primeError <- error
  primeDistances <- distances
  #initializes the iteration index with 1
  i <- 1
  #initializes the energy with 0
  energy <- 0
  
  samples <- 0
  refSamples <- 0
  
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
    nSamples <- floor(m * pSample)
    dX <- round(1/pSample)
    
    samples <- 0
    refSamples <- 1:length(reference[,2])
    
    if(pSample > 0){
      if(m %% 2 == 0)
        samples <- 0:(nSamples - 1) * dX + sample(1:dX, 1)
      else{
        if(dX > 1)
          samples <- 0:(nSamples - 1) * dX + sample(1:(dX-1), 1)
        else
          samples <- 0:(nSamples - 1) * dX + 1
      }
    }
    else{
      
      samples <- detectPeaks(target[,2], 3)
      refSamples <- detectPeaks(reference[,2], 3)
      
      newSamples <- getCorrespondents(target[samples,], reference[refSamples,])
      
      samples <- samples[newSamples[[1]]]
      refSamples <- refSamples[newSamples[[2]]]
      
      nSamples <- length(samples)
      #cat(nSamples, "and", length(refSamples), "samples\n")
    }
    
    #cat("m: ", m, "dX:", dX, "samples:", samples, "\n")
    
    translationFactorX <- 0
    translationFactorY <- 0
    
    if(isOpt){
      
      data <- list(target = target, reference = reference[refSamples,], translationFactorX = 0, translationFactorY = 0)
      
      data <- Reduce(function(data, sample){
        
        dists <- as.matrix(dist(rbind(data$target[sample,], data$reference)))[1,-1]
        k <- which.min(dists)
        
        data$translationFactorX <- data$translationFactorX + data$reference[k,1] - data$target[sample,1]
        data$translationFactorY <- data$translationFactorY + data$reference[k,2] - data$target[sample,2]
        
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
        
        cat(translationFactorX, "; ", translationFactorY, "\n")
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
    #distances <- dist.p2p(reference, target[samples,], isOpt, 2)
    distances <- dist.p2p(reference, target, isOpt)
    #distances <- dist.p2p(reference[refSamples,], target[samples,], isOpt, 2)
    
    #checks whether the curves got too far
    if(commonDomain(reference, target, isOpt) >= threshold)
      #if they didn't, measures the error
      error <- mean(abs(distances[,2]))
    else
      #otherwise, sets the erro to the prime error plus 1
      error <- primeError + 1
    
    #cat("Iteration ", i, "; error = ", error, "\n")
    #increasing the iteration index
    i <- i + 1
  }
  #returns the informations
  (list(target = primeTarget, dist = primeDistances, error = primeError, energyTotal = energy,
        energyMean = (energy/(i - 1)), samples=samples, refSamples=refSamples))
}

getCorrespondents <- function(P1, P2){
  
  reference <- 1
  target <- 1
  
  n1 <- length(P1[,1])
  n2 <- length(P2[,1])
  
  if(n1 > n2){
    
    target <- P1
    reference <- P2
  }
  else if(n1 < n2){
    target <- P2
    reference <- P1
  }
  else{
    return(list(1:n1, 1:n2))
  }
  
  x <- rep(0, length(reference[,1]))
  
  minRef = min(reference[,2])
  maxRef = max(reference[,2])
  minTar = min(target[,2])
  maxTar = max(target[,2])
  
  target[,2] = ((target[,2] - minTar)*(maxRef-minRef))/(maxTar - minTar) + minRef

#   xMeanRef = mean(reference[,1])
#   yMeanRef = mean(reference[,2])
#   xMeanTar = mean(target[,1])
#   yMeanTar = mean(target[,2])
# 
#   target[,2] <- target[,2] + yMeanRef - yMeanTar
#   target[,1] <- target[,1] + xMeanRef - xMeanTar
  
  Ds <- mapply(function(x,y,target){
    
    ds <- as.matrix(dist(rbind(matrix(c(x,y), nrow=1), target)))[-1,1]
    return(min(ds))
    
  }, reference[,1], reference[,2], MoreArgs = list(target=target))
  
  Ds <- mean(Ds)

  tarInit <- 1
  tarN <- length(target[,1])
  
  for(i in 1:(length(reference[,1]))){
    
    ds <- as.matrix(dist(rbind(reference[i,], target[tarInit:tarN,])))[1,-1]
      
#     ok <- TRUE
#     
#     while(ok){
#       
#       tarI <- which.min(ds) + tarInit - 1
#       wTar <- getWindow(target[,1], tarI, 1)
#       wRef <- getWindow(reference[,1], i, 1)
#       
#       if(min(length(wTar), length(wRef)) < 2){
#         
#         aTar <- 0
#         aRef <- 0
#         
#         if(length(wTar) < 2){
#           if(tarI == tarInit)
#             aTar <- findLine(target[tarI,], target[wTar,])$a
#           else
#             aTar <- findLine(target[wTar,], target[tarI,])$a
#           }
#         else{
#           
#           if(i == 1)
#             aTar <- findLine(target[tarI,], target[wTar[2],])$a
#           else
#             aTar <- findLine(target[wTar[1],], target[tarI,])$a
#         }
#         
#         if(length(wRef) < 2){
#           if(i == 1)
#             aRef <- findLine(reference[i,], reference[wRef,])$a
#           else
#             aRef <- findLine(reference[wRef,], reference[i,])$a
#         }
#         else{
#           
#           if(tarI == tarInit)
#             aRef <- findLine(reference[i,], reference[wRef[2],])$a
#           else
#             aRef <- findLine(reference[wRef[1],], reference[i,])$a
#         }
#         
#         if(aTar > aRef*0.66 && aTar < aRef*1.5){
#           
#           if(min(ds) <= Ds){
#             x[i] <- tarI
#             tarInit <- tarI
#             target[x[i],2] <- 1000000
#           }
#           else
#             reference[i,1] <- 0.5
#           
#           ok <- FALSE
#         }
#         else{
#           
#           ds[which.min(ds)] <- 10000000
#         }
#       }
#       else{
#         
#         a1Tar <- findLine(target[wTar[1],], target[tarI,])$a
#         a2Tar <- findLine(target[tarI,], target[wTar[2],])$a
#         a1Ref <- findLine(reference[wRef[1],], reference[i,])$a
#         a2Ref <- findLine(reference[i,], reference[wRef[2],])$a
#         
#         if(a1Tar > a1Ref*0.66 && a1Tar < a1Ref*1.5 && a2Tar > a2Ref*0.66 && a2Tar < a2Ref*1.5){
#           
#           if(min(ds) <= Ds){
#             x[i] <- tarI
#             tarInit <- tarI
#             target[x[i],2] <- 1000000
#           }
#           else
#             reference[i,1] <- 0.5
#           
#           ok <- FALSE
#         }
#         else{
#           
#           ds[which.min(ds)] <- 10000000
#         }
#       }
#       if(min(ds) == 10000000){
#         ok <- FALSE
#         reference[i,1] <- 0.5
#       }
#     }
    
    if(min(ds) <= Ds){
      
      x[i] <- which.min(ds)
      target[x[i],2] <- 100000
    }
    else{
      
      reference[i,1] <- 0.5
    }
  }
  
  toRemove <- which(reference[,1] == 0.5)
  y <- 1:(length(reference[,1]))
  
  if(length(toRemove) > 0){
    y <- y[-toRemove]
    x <- x[-toRemove]
  }
  
  if(length(y) < 2){
    
    return(list(1:(length(reference)), 1:(length(reference))))
  }
  
  if(n1 > n2){
    
    return(list(x, y))
  }
  else if(n1 < n2){
    
    return(list(y, x))
  }
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
  return(limit(k, 0, "lower"))
}

difference <- function(x, y){
  
  return (x - y);
}

difference2 <- function(x, y, ix, iy, mx=0, my=0){
  
  dif <- x - y
  #ix <- which(mx[,1] == ix)
  #iy <- which(my[,1] == iy)
  if(is.matrix(mx) && is.matrix(my)){
    
    if(abs(dif) <= 1)
      return(dif)
    else{
      
      n <- floor(abs(dif))
      candidates <- limit(ix-n, 1, "lower"):limit(ix+n, length(mx[,2]), "upper")
      cand <- candidates
      #cat("candidates: ", candidates, "mx dims:", dim(mx), "\n")
      
      if(dif < 0){
        candidates <- candidates[which(mx[candidates,2] >= x)]
	
	if(length(candidates) == 1)
	   return(dif)

        dists <- as.matrix(dist(rbind(my[iy,], mx[candidates,])))[1,-1]
        return(-min(dists))
      }
      else{
        candidates <- candidates[which(mx[candidates,2] <= x)]

	if(length(candidates) == 1)
	   return(dif)

        dists <- as.matrix(dist(rbind(my[iy,], mx[candidates,])))[1,-1]
        return(min(dists))
      }
    }
  }
  else
    return(dif)
}

# Computes the distance in 'y' for each point by the target's coordinate
dist.p2p <- function(reference, target, isOpt=FALSE, type=1){
  
  n <- length(target[,1])
  
  distances <- rep(0, n)
  notPresent <- 0
  xs <- 0
  
  if(isOpt){
    
    m <- length(reference[,1])
    
    xmin = 1
    xmax = length(target[,1])
    
    ref <- xmin:xmax
    tar <- xmin:xmax
    
    if(type == 1){
    
      xmin <- max(target[1,1], reference[1,1])
      xmax <- min(reference[m,1], target[n,1])
      
      ref <- (Position(function(x){return(x == xmin)}, reference[,1]):Position(function(x){return(x == xmax)}, reference[,1]))
      tar <- (Position(function(x){return(x == xmin)}, target[,1]):Position(function(x){return(x == xmax)}, target[,1]))
    }
    else{
      
      ref <- mapply(function(x, i, y){
        
        if(length(which(x == y)) > 0)
          
          return(i)
        else
          return(0)
        
      }, reference[,1], c(1:length(reference[,1])), MoreArgs = list(y=target[,1]))
      
      ref <- ref[-which(ref == 0)]
      xmin <- min(reference[ref,1])
      xmax <- max(reference[ref,1])
      
      toRemove <- which(target[tar,1] > xmax | target[tar,1] < xmin)
      
      if(length(toRemove) > 0)
        tar <- tar[-toRemove]
    }
    
    #distances <- mapply(difference, reference[ref,2], target[tar,2])
    #distances <- mapply(difference2, reference[ref,2], target[tar,2], reference[ref,1], target[tar,1], MoreArgs=list(reference, target))
    distances <- mapply(difference2, reference[ref,2], target[tar,2], ref, tar, MoreArgs=list(reference, target))
#     distances <- mapply(function(reference, target, ref, tar){
#       
#       m = matrix(c(ref, tar, reference, target), ncol=2)
#       d = as.matrix(dist(m))[-1,1]
#       
#       if(target > reference)
#         return(-d)
#       else
#         return(d)
#       
#     }, reference[ref,2], target[tar,2], reference[ref,1], target[tar,1])
    distances <- matrix(c(target[tar,1], distances), ncol=2)
  }
  else{
    for(i in 1:n){
      
      x <- target[i,1]
      
      if(length(which(reference[,1] == x)) > 0){
        distances[i] <- reference[which(reference[,1] == x), 2] - target[i,2]
        xs <- c(xs, x)
      }
      else
        notPresent <- c(notPresent, i)
    }
    
    notPresent <- notPresent[-1]
    if(length(notPresent) > 0)
      distances <- distances[-notPresent]
    
    distances <- matrix(c(xs, distances), ncol=2)
  }
  
  (distances)
}

# Removes the points with black/zero value.
# input:
#   data = a 2D matrix containing a line/curve
# output:
#   a 2D matrix containing the biggest part of 'data' which
#   doesn't contain black/zero value.
takeNoneFaceOut <- function(data, takeOffSteep=FALSE){
  
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
  
  if(takeOffSteep){
    
    m <- length(data[,1])
    angularCoefficients <- mapply(angularCoeff, data[1:(m-1),1], data[1:(m-1),2], data[2:m,1], data[2:m,2])
    angularCoefficients <- abs(angularCoefficients)
    threshold <- 3 * mean(angularCoefficients)
    
    xout <- which(angularCoefficients > threshold)
    init <- 1
    while(length(xout) > 0 && xout[1] == init){
      
      data <- data[-1,]
      xout <- xout[-1]
      init <- init + 1
    }
     
    init <- m-1
    while(length(xout) > 0 && xout[length(xout)] == init){
      
      data <- data[-(length(data[,1])),]
      xout <- xout[-(length(xout))]
      init <- init - 1
    }
  }
  
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
  #curve <- interpolateXinteger(curve)
  
  #returns the curve rotated
  (curve)
}

#discretizePoint <- function(x1, y1, x2, y2, m){
discretizePoint <- function(xs, m){
  
  i <- Position(function(x){return(x > xs)}, m[,1])
  
  if(is.na(i))
    #i <- length(m[,1])
    #i <- which.max(m[,1])
    return(c(max(m[,1]), NA))
  else if(i == 1)
    i <- 2
  
  #cat("i:", i, " mi-1:", m[i-1,], " mi:", m[i,], "\n")
  
  if(xs == m[i,1])
    return(m[i,])
  
  line <- findLine(m[i-1,], m[i,])
  
  y <- appLinear(line, xs)
  
  return(c(xs, y))
  
  #   if(x1 + 0.5 < x2){
  #     #computes the line which passes through the 'i'th point and its successor
  #     line <- findLine(c(x1, y1), c(x2, y2))
  #     
  #     #finds the correspoding value (2nd column) to this new domain value
  #     x <- round(x1)
  #     y <- appLinear(line, x)
  #     #applies the results
  #     return (c(x,y))
  #   }
  #   else
  #     return (c(x1,y1))
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
    #N <- floor(m[n,1]) - x
    
    if(!AND(m[,1] == x:(x + n - 1))){
      
      #m2 <- t(mapply(discretizePoint, m[1:(n-1),1], m[1:(n-1),2], m[2:n,1], m[2:n,2]))
      m2 <- t(mapply(discretizePoint, x:(x + n - 1), MoreArgs=list(m = m)))
      m2 <- m2[which(!is.na(m2[,2])),]
      return(m2)
    }
    else
      return (m)
  }
  else{
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
    
#     x <- 0
#     if(isOpt)
#       x <- m2[n-1,1] + 1
#     else
      x <- ceiling(m[n,1])
    
    y <- appLinear(line, x)
    m[n,1] <- x
    m[n,2] <- y
  }}
  
  #returns the interpolated curve
#   if(isOpt)
#     m <- rbind(m2, m[n,])
  
  #removes the remaining floating domains
  #m <- m[which(m[,1] == ceiling(m[,1])),]

  (m)
}
