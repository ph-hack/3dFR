angleBetween <- function(c1, c2){
  
  n1 <- length(c1[,1])
  n2 <- length(c2[,1])
  
  d1 <- mapply(difference, c1[1:(n1-1),2], c1[2:n1,2])
  d2 <- mapply(difference, c2[1:(n2-1),2], c2[2:n2,2])
  
#   angles <- mapply(function(x,y){
#     
#     return(atan(x - y))
#     
#   }, d1, d2)
#   
#   result <- mean(angles)

  a1 <- mean(atan(d1))
  a2 <- mean(atan(d2))

  result <- a1 - a2
  
#   #computes the descriptor lines of both curves
#   line1 <- linearInterpolation(c1, TRUE)
#   line2 <- linearInterpolation(c2, TRUE)
#   
#   #computes the angle between them
#   angle <- atan(line1$a - line2$a)
#   
#   return(c(result, angle))

    return(c(result, 1))
}

detectPeaks <- function(curve, smooth=0, isToPlot=FALSE){
  
  n <- length(curve)
  M <- 1
  g <- rep(0, n)
  
  while(M >= 0.09){
    if(smooth > 0 && smooth < n/2)
    scurve <- gaussianSmooth(c(rep(curve[1], 2*smooth), curve, rep(curve[n], 2*smooth)), c(smooth))[(1 + 2*smooth):(n + 2*smooth)]
    
    g <- gradient(scurve)
    g <- gradient(g)
    M <- mean(g)
    smooth <- smooth + 1
  }
  
  P <- mapply(function(x){
    
    w <- getWindow(g, x, 1)
    
    if(length(w) == 1){
        return (g[x] > g[w])
    }
    else{
      return ((g[w[1]] < g[x] && g[x] > g[w[2]]))# ||
                #(g[w[1]] > g[x] && g[x] < g[w[2]]))
    }
      
  }, 1:n)
  
  P <- which(P)
  
  if(isToPlot){
    
    curve <- setRange(curve)
    g <- setRange(g)
    
    plot(curve, type="l")
    lines(g, col="blue")
    
    for(p in P)
      lines(rep(p, 2), c(0,1), col="red")
  }
  
  return(P)
}

gradient <- function(curve, absolute=TRUE){
  
  N <- length(curve)
  G <- mapply(function(x, y){
    
    w <- getWindow(curve, y, 1)
    
    g <- sum(abs(x - curve[w]))/length(w)
    
    if(!absolute){
      
      if(w[1] > y){
        
        if(x - curve[w[1]] > 0){
          
          g = -g
        }
      }
      else{
        
        if(x - curve[w[1]] < 0){
          
          g = -g
        }
      }
    }
    
    return(g)
    
  }, curve, 1:N)
  
  return(G)
}

getWindow <- function(signal, idx, size){
  
  N <- length(signal)
  
  W <- c()
  
  if(idx > N || idx < 1 || size <= 0)
    return(W)
  
  if(idx - size > 0)
    W <- c(W, (idx - size):(idx - 1))
  
  if(idx + size <= N)
    W <- c(W, (idx + 1):(idx + size))
  
  return(W)
}

setRange <- function(x, range=c(0,1)){
  
  minX <- min(x)
  maxX <- max(x)
  
  x <- ((x - minX)*(range[2]-range[1]))/(maxX - minX) + range[1]
  
  return(x)
}

my.cosineDistance <- function(reference, target, method = "dct", coeffRange=c()){
  
  N <- length(reference)
  
  ref = dtt(reference, method)
  tar = dtt(target, method)
    
  if(length(coeffRange) > 0){
    
    ref = ref[coeffRange]
    tar = tar[coeffRange]
  }
  
  error = 1 - sum(ref * tar)/sqrt(sum(ref^2) * sum(tar^2))
  
  dist = apply(matrix(c(ref, tar), nrow=2, byrow=TRUE), 2, function(x){ return(x[1] - x[2])})
  
  return(list(error=error, dist=my.as.matrix(dist), target=matrix(c(1:(length(target)), target), ncol=2)))
}

applyDCT <- function(listData, meanCurves=list()){
  
  N <- length(listData)
  
  colors <- c('black', 'red', 'blue', 'green', 'cyan', 'orange', 'brown', 'gray', 'yellow')
  
  if(N > 9)
    N <- 9
  
  plot(matrix(c(0,20,-1000,500), ncol=2), col="white")
  
  DCT <- list()
  
  for(i in 1:N){
    
    M <- dim(listData[[i]])[1]
    
    dcts <- apply(listData[[i]], 1, function(x){
      
      #return(dtt(curveCorrection3(x, meanCurves[[i]], 1), "dct"))
      return(dtt(x, "dct"))
    })
    
    for(j in 1:M){
      
      lines(dcts[2:20, j], col=colors[i])
    }
    
    DCT[[i]] <- dcts
  }
  
  return(DCT)
}

hierarchicalFeatureBasedPrediction3 <- function(model, testDir="", testing=character(0), subset=integer(0),
                                                useErrorRange=TRUE, logFile="", errorFunctions=list(my.icp.2d.v2),
                                                errorParams=list(list(minIter=4, pSample=0.2)), weightLimit=0.5,
                                                evaluate=FALSE, weights=c(), isToPlot=FALSE){
  
  #gets the files' names
  if(length(testing) == 0)
    testing <- dir(testDir)
  
  if(length(subset) > 0)
    testing <- testing[subset]
  
  #retrieves the classes information
  classes <- getClassFromFiles(files=testing)
  
  #gets the descriptors' values for all test samples
  descriptors <- lapply(concatenate(list(testDir, testing)), readMainLines, "list")
  
  #the number of testing sample
  M <- length(testing)
  #the number of different error functions
  Nef = length(errorFunctions)
  #the number of models, one for each descriptor
  N <- length(model)/Nef
  
  tests <- list()
  
  #for each descriptor, ...
  for(i in 1:N){
    
    #separates only the vectors for the ith descriptor
    samples <- getAllFieldFromList(descriptors, i, 2)
    #puts them into a matrix
    samples <- list2matrix(samples)
    #puts this matrix into tests list
    if(i < 4){
      tests[[i]] <- samples[,1:100]
      dim(tests[[i]]) <- c(length(tests[[i]])/100,100)
    }
    else
      tests[[i]] <- samples
  }
  
  corrects <- 0
  
  votesByDescriptor <- list()
  
  cat("Testing", M, "samples!\n", file=logFile, append=FALSE)
  
  #for each test sample, ...
  for(m in 1:M){
    
    cat("\npredicting test", m, ":", file=logFile, append=TRUE)
    
    start <- getTime()
    
    comparisons <- 0
    
    #initializes the votes as an empty list
    votes <- list()
    votesByDescriptor[[testing[m]]] <- list()
    
    #for each error function, ...
    for(f in 1:Nef){
      
      #initiates a vote list for each error function
      votes[[f]] <- list()
      
      #for each descriptor, ...
      for(i in 1:N){
        
        cat("\nWith predictor", i, file=logFile, append=TRUE)
        
        #initializes this descriptor's votes as a matrix with zeros (these zeros will be ignored later)
        #votes[[(f-1)*N + i]] <- matrix(c(0,0), nrow=1)
        votes[[f]][[i]] <- matrix(c(0,0), nrow=1)
        voteWeights <- c()
        maxWeightErrors <- c()
        
        #gets the ith descriptor of the mth test sample
        test <- tests[[i]][m,]
        
        #gets the upper level's representants
        level <- list2matrix(getAllFieldFromList(model[[i]], "representant", 2))
        
        levelIndex <- list(c(0, 1))
        
        levelQueue <- list(model[[(f-1)*N + i]])
        
        #for each node from second level which matched the test, ...
        while(length(levelQueue) > 0){
          
          branch <- levelQueue[[1]]
          
          if(levelIndex[[1]][1] != 0)
            cat("\n", concatenate(rep("   ", levelIndex[[1]][1])), levelIndex[[1]][1], "level, node", levelIndex[[1]][2], file=logFile, append=TRUE)
          
          #gets the first level's representants
          representants <- list2matrix(getAllFieldFromList(branch, "representant", 2))
          dists <- 1
          targets <- 1
          refPoints <- 1
          tarPoints <- 1
          errors <- c()
          
          passed <- 1
          
          #if(i == 5 && f == 1){
          #  cat("\ngot it\n")
          #  cat("\n")
          #}
          
          #if there is more than one representant at this level
          if(!is.null(dim(representants))){
            #computes the errors for the first level
            icpResults <- apply(representants, 1, function(reference, target){
              
              return (do.call(errorFunctions[[f]], merge.list(list(reference, curveCorrection3(target, reference, 1)), errorParams[[f]])))
              #return (dtw(reference, curveCorrection3(target, reference, 1)))
              
            }, test)
            
            errors <- list2vector(getAllFieldFromList(icpResults, "error", 2))
            targets <- getAllFieldFromList(icpResults, "target", 2)
            #if(f == 1){
            #  refPoints <- getAllFieldFromList(icpResults, "refSamples", 2)
            #  tarPoints <- getAllFieldFromList(icpResults, "samples", 2)
            #}
            #errors <- list2vector(getAllFieldFromList(icpResults, "normalizedDistance", 2))
            #maxErrors <- list2vector(getAllFieldFromList(branch, "maxError", 2))
            ws <- list2vector(getAllFieldFromList(branch, "weight", 2))
            
            rangeCheck <- rep(TRUE, length(errors))
            
            if(useErrorRange){
              
              dists <- getAllFieldFromList(icpResults, "dist", 2)
              errorRange <- getAllFieldFromList(branch, "errorRange", 2)
              rangeCheck <- checkErrorRange(errorRange, dists)
            }
            #retrieves which nodes of the first level matched the test
            #passed <- which(errors <= maxErrors & rangeCheck)
            wFailed = which(ws > rep(weightLimit, length(ws)))
            if (length(wFailed) > 0)
              cat(" [Removed by weight: ", wFailed, "]", file = logFile, append=TRUE)
            rFailed = which(!rangeCheck)
            if(length(rFailed) > 0)
              cat(" [Removed by error range: ", rFailed, "]", file = logFile, append=TRUE)
            
            passed <- which(ws <= rep(weightLimit, length(ws)) & rangeCheck)
          }
          else{
            #computes the error with the single first level's representant
            icpResults <- do.call(errorFunctions[[f]], merge.list(list(representants, curveCorrection3(test, representants, 1)), errorParams[[f]]))
            #icpResults <- dtw(representants, curveCorrection3(test, representants, 1))
            targets <- list(icpResults$target)
            errors <- icpResults$error
            #if(f == 2){
            #  refPoints <- list(icpResults$refSamples)
            #  tarPoints <- list(icpResults$samples)
            #}
            #maxError <- list2vector(getAllFieldFromList(branch, "maxError", 2))
            ws <- list2vector(getAllFieldFromList(branch, "weight", 2))
            
            rangeCheck <- TRUE
            
            if(useErrorRange){
              
              errorRange <- getAllFieldFromList(branch, "errorRange", 2)
              rangeCheck <- checkErrorRange(errorRange, list(icpResults$dist))
              dists <- list(icpResults$dist)
            }
            
            #checks whether the representant matched
            #passed <- which(icpResults$error <= maxError & rangeCheck)
            if (ws > weightLimit)
              cat(" [Removed by weight: ", 1, "]", file = logFile, append=TRUE)
            
            if(!rangeCheck)
              cat(" [Removed by error range: ", 1, "]", file = logFile, append=TRUE)
            
            passed <- which(ws <= weightLimit & rangeCheck)
          }
          
          comparisons <- comparisons + dim(representants)[1]
          
          #if this is a leaf, ...
          if(is.null(branch[[1]]$children)){
            
#             for(v in 1:(length(branch))){
#               if(isToPlot && length(passed) == 0 && names(branch)[v] == classes$fileClasses[m]){
#                 
#                 ymax = max(c(branch[[v]]$errorRange[,3], dists[[v]][,2]))
#                 ymin = min(c(branch[[v]]$errorRange[,2], dists[[v]][,2]))
#                 xmax = max(branch[[v]]$errorRange[,1], dists[[v]][,1])
#                 xmin = min(branch[[v]]$errorRange[,1], dists[[v]][,1])
#                 
#                 #plots the error range graph
#                 plot(c(1,xmax), c(ymin, ymax), col="white", main=concatenate(c("Error Range predictor ", f, ", ", i)))
#                 lines(c(1,xmax), c(0,0), col="gray")
#                 sorted <- sort.int(branch[[v]]$errorRange[,1], index.return = TRUE)$ix
#                 meanRange = rowMeans(branch[[v]]$errorRange[(sorted),-1])
#                 lines(x = branch[[v]]$errorRange[sorted,1], y = meanRange, col="green")
#                 lines(x = branch[[v]]$errorRange[sorted,1], y = branch[[v]]$errorRange[sorted,2], col="red")
#                 lines(x = branch[[v]]$errorRange[sorted,1], y = branch[[v]]$errorRange[sorted,3], col="red")
#                 lines(dists[[v]], col="black")
#                 
#                 plot(branch[[v]]$representant, type="l", col="red", main=concatenate(c("Curves predictor ", f, ", ", i)))
#                 lines(test, col="blue")
#                 lines(x = targets[[v]][,1], y = targets[[v]][,2], col="black")
#                 if(f == 1){
#                   points(refPoints[[v]], branch[[v]]$representant[refPoints[[v]]], col="red")
#                   points(targets[[v]][tarPoints[[v]],], col="black")
#                 }
#               }
#             }
            
            for(v in passed){
              
              cat("\n", concatenate(rep("   ", levelIndex[[1]][1])), " leaf(", names(branch)[v], ")", file=logFile, append=TRUE)
              
              #gets the leaf's samples
              #samples <- branch[[v]]$samples
              #computes the errors for each leaf sample
              #icpResults <- apply(samples, 1, function(reference, target){
                
              #  return (do.call(errorFunctions[[f]], merge.list(list(curveCorrection3(reference, representants[v,], 1),
              #                                                       curveCorrection3(target, representants[v,], 1)),
              #                                                  errorParams[[f]])))
                #return (dtw(curveCorrection3(reference, representants[v,], 1), curveCorrection3(target, representants[v,], 1)))
                
              #}, test)
              #gets the minimum computed error
              #minErrorIndex <- which.min(list2vector(getAllFieldFromList(icpResults, "error", 2)))
              #minErrorIndex <- which.min(list2vector(getAllFieldFromList(icpResults, "normalizedDistance", 2)))
              #minError <- list2vector(getAllFieldFromList(icpResults, "error", 2))[minErrorIndex]
              minError <- errors[v]
              #minError <- mean(list2vector(getAllFieldFromList(icpResults, "error", 2))) #/branch[[v]]$maxError
              #minError <- list2vector(getAllFieldFromList(icpResults, "normalizedDistance", 2))[minErrorIndex]
              
              #cat(" -------", minErrorIndex, "------", file=logFile, append=TRUE)
              cat(" E =", minError, file=logFile, append=TRUE)
              
              maxWeightErrors <- c(maxWeightErrors, branch[[v]]$deviation + branch[[v]]$meanError)
              
              if(length(weights) > 0)
                if(!is.null(weights[[names(branch)[v]]]))
                  minError <- applyWeight(minError, weights[[names(branch)[v]]][i], branch[[v]]$deviation + branch[[v]]$meanError)
              else
                minError <- minError * 20
              #else
              #  minError <- applyWeight(minError, branch[[v]]$weight, branch[[v]]$deviation + branch[[v]]$meanError)
              
              #cat(" Ew =", minError, file=logFile, append=TRUE)
              
              #adds a vote for this leaf's class with the weight as the minimum error value
              #cat("leaf:", v, " descriptor:", i, "test:", m, "first level:", k, "second level:", j, "\n")
              if(minError <= branch[[v]]$maxError){
                
                #votes[[(f-1)*N + i]] <- rbind(votes[[(f-1)*N + i]], matrix(c(as.numeric(names(branch)[v]), minError), nrow=1))
                votes[[f]][[i]] <- rbind(votes[[f]][[i]], matrix(c(as.numeric(names(branch)[v]), minError), nrow=1))
                voteWeights <- c(voteWeights, branch[[v]]$weight)
              }
              else{
                
                cat(" Failed due to max error", file=logFile, append=TRUE)
              }
            }
          }
          else{
            
            for(v in passed){
              
              #levelQueue[[length(levelQueue) + 1]] <- branch[[v]]$children
              
              if(length(levelQueue) >= 2){
                levelQueue <- merge.list(list(levelQueue[[1]], branch[[v]]$children), levelQueue[2:length(levelQueue)])
                levelIndex <- merge.list(list(levelIndex[[1]], c(levelIndex[[1]][1] + 1, v)), levelIndex[2:length(levelIndex)])
              }
              else{
                levelQueue[[length(levelQueue) + 1]] <- branch[[v]]$children
                levelIndex[[length(levelIndex) + 1]] <- c(levelIndex[[1]][1] + 1, v)
              }
              
              #levelIndex[[length(levelIndex) + 1]] <- c(levelIndex[[1]][1] + 1, v)
            }
          }
          
          levelQueue <- levelQueue[-1]
          levelIndex <- levelIndex[-1]
        }
        
        #removes the initialization value
        #votes[[(f-1)*N + i]] <- votes[[(f-1)*N + i]][-1,]
        votes[[f]][[i]] <- votes[[f]][[i]][-1,]
        
        #Normalizatioin by the max of this descriptors votes
#         if(length(votes[[(f-1)*N + i]]) > 2)
#           votes[[(f-1)*N + i]][,2] <- votes[[(f-1)*N + i]][,2]/max(votes[[(f-1)*N + i]][,2])
#         else
#           if(length(votes[[(f-1)*N + i]]) > 1)
#             votes[[(f-1)*N + i]][2] <- votes[[(f-1)*N + i]][2]/max(votes[[(f-1)*N + i]][2])
        
        #dim(votes[[i]]) <- c(length(votes[[i]])/2, 2)
        
        #applies the weights to the votes
        #if(length(votes[[i]][,1]) > 1){
        
        #  votes[[i]][,2] <- mapply(applyWeight, votes[[i]][,2], voteWeights, MoreArgs = list(maxError = mean(votes[[i]][,2])))
        #}
        #else{
        
        #  votes[[i]][,2] <- mapply(applyWeight, votes[[i]][,2], voteWeights, maxWeightErrors)
        #}
        
        
        #votesByDescriptor[[testing[m]]][[(f-1)*N + i]] <- votes[[(f-1)*N + i]]
        votesByDescriptor[[testing[m]]][[(f-1)*N + i]] <- votes[[f]][[i]]
      }

#       votes[[f]] <- Reduce(rbind, votes[[f]], matrix(c(0,0), nrow=1))[-1,]
#       dim(votes[[f]]) <- c(length(votes[[f]])/2, 2)
#   
#       if(is.null(dim(votes[[f]])))
#         dim(votes[[f]]) <- c(1,2)
    }
    
    #counts the votes
    #votes <- Reduce(rbind, votes, matrix(c(0,0), nrow=1))[-1,]
    
    #dim(votes) <- c(length(votes)/2, 2)
    
    cat("\nvotes:\n", file=logFile, append=TRUE)
#     votesSizes <- rep(0, Nef)
#     for(f in 1:Nef){
#       cat.matrix(votes[[f]], file=logFile, append=TRUE)
#       votesSizes[f] <- length(votes[[f]][,1])
#       cat('------------- number: ', votesSizes[f], '\n', file=logFile, append=TRUE)
#     }
#     cat("number of votes:", sum(votesSizes), "\n", file=logFile, append=TRUE)

    finalVotes <- votes[[1]]

#     for(i in 1:N){
#       
#       if(length(finalVotes[[i]]) == 2)
#         dim(finalVotes[[i]]) <- c(1,2)
#       
#       if(length(votes[[2]][[i]]) == 2)
#         dim(votes[[2]][[i]]) <- c(1,2)
#       
#       Nv <- length(finalVotes[[i]][,2])
#       
#       if(Nv > 0){
#         for(j in 1:Nv){
#           
#           p <- Position(function(x){ return(x == finalVotes[[i]][j,1])}, votes[[2]][[i]][,1], nomatch = 0)
#           
#           if(p > 0){
#             
#             cat("desc", i, "class", finalVotes[[i]][j,1], ":", finalVotes[[i]][j,2], "+", votes[[2]][[i]][p,2], "\n", file = logFile, append = TRUE)
#             finalVotes[[i]][j,2] <- finalVotes[[i]][j,2] + votes[[2]][[i]][p,2]
#           }
#           else{
#             
#             cat("desc", i, "class", finalVotes[[i]][j,1], ":", finalVotes[[i]][j,2], "+ ", "1\n", file = logFile, append = TRUE)
#             finalVotes[[i]][j,2] <- finalVotes[[i]][j,2] + 1
#           }
#         }
#       }
#     }

    finalVotes <- Reduce(rbind, finalVotes, matrix(c(0,0), nrow=1))[-1,]
    dim(finalVotes) <- c(length(finalVotes)/2, 2)
    
    if(is.null(dim(finalVotes)))
      dim(finalVotes) <- c(1,2)
    
    cat("votes number:", length(finalVotes[,1]), "\n", file=logFile, append=TRUE)
    
    #if(is.null(dim(votes)))
    #  dim(votes) <- c(1,2)

#     finalVotes <- matrix(rep(0,Nef*2), nrow=Nef)
# 
#     for(f in 1:Nef){
#       
#       if(length(votes[[f]]) > 0){
#         result <- ponderateVote(votes[[f]], "min", "value")
#         finalVotes[f,1] <- result[1]
#         finalVotes[f,2] <- result[3]
#       }
#     }

    cat("final votes:\n", file=logFile, append=TRUE)
    cat.matrix(finalVotes, file=logFile, append=TRUE)
    cat("comparisons:", comparisons, "\n", file=logFile, append=TRUE)
    cat("test", m, ". ", file=logFile, append=TRUE)
    
    if(length(finalVotes) > 1){
      
      result <- ponderateVote(finalVotes, by="min", type="value")
      #checks the result
      if(paste("0", as.character(result[1]), sep="") == classes$fileClasses[m]){
        
        cat("Found with", result[2], "\n", file=logFile, append=TRUE)
        corrects <- corrects + 1
      }
      else
        cat("Missed with", result[2], "\n", file=logFile, append=TRUE)
      
      cat("Result:", paste("0", as.character(result[1]), sep=""), file=logFile, append=TRUE)
    }
    else
      cat("Unknown!\n", file=logFile, append=TRUE)
    
    cat(" Expected:", classes$fileClasses[m], "\n", file=logFile, append=TRUE)
    cat("time: ", crono.end(start), "\n", file=logFile, append=TRUE)
    
    cat("\n", file=logFile, append=TRUE)
  }
  
  cat("--------Accuracy:", corrects/M*100, "----------\n", file=logFile, append=TRUE)
  
  if(evaluate)
    return(votesByDescriptor)
}

testHierarchicalFamiliarityWeight <- function(smaller, bigger, type=1){
  
  N <- length(smaller[,1])
  
  for(i in 1:N){
    
    s <- hierarchicalFamiliarityWeight(smaller[i,1], smaller[i, 2], smaller[i, 3], type)
    b <- hierarchicalFamiliarityWeight(bigger[i,1], bigger[i, 2], bigger[i, 3], type)
    
    if(s >= b){
      
      cat("\nError!\n")
      cat("smaller: ", smaller[i,], "->", s, "\n")
      cat("bigger: ", bigger[i,], "->", b, "\n")
    }
  }
  cat("\nDone!")
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
