library(adabag)

hierarchicalFeatureBasedClassifier <- function(trainingDir, training=c(), groupNumbers=c(7,3),
                                               features=c(1:11), errorParams=list(range=20, smooth=5, tol=10)){
  
  if(length(training) == 0)
    training <- dir(trainingDir)
  
  #retrieves the classes information
  classes <- getClassFromFiles(files=training)
  
  descriptors <- lapply(concatenate(list(trainingDir, training)), readMainLines, "list")
  
  N <- length(features) #descriptors[[1]])
  C <- length(classes$classes)
  
  model <- list()
  samples <- list()
  
  for(i in features){
    #separates only the vectors for the ith descriptors
    samples[[i]] <- getAllFieldFromList(descriptors, i, 2)
    #puts them into a matrix
    samples[[i]] <- list2matrix(samples[[i]])
    
    if(i < 4){
      
      samples[[i]] <- samples[[i]][,1:100]
    }
  }
  
  #creates the first C nodes, where C = number of classes
  currentLevel <- computeNodes(samples, classes$fileClasses, features, errorParams = errorParams, progress=TRUE)
  
  fittingResults <- fitWholeLevel(currentLevel, errorParams, TRUE)
  
  nGroups <- length(groupNumbers)
  
  for(g in groupNumbers){
    
    #divides the leafs into groups
    gResult <- computeGrouping(fittingResults, g, progress=TRUE)
    currentLevel <- gResult$level
    #mounts the first level
    
    levelSamples <- lapply(1:N, function(x, y, z){
      
      #return(list2matrix(y[(1+((x-1)*z)):(z*x)]))
      return(list2matrix(y[,x]))
      
    }, list2matrix(getAllFieldFromList(currentLevel, "representant", 2)), C)
    
    
    
    currentLevel <- computeNodes(levelSamples, gResult$groups, features, currentLevel, errorParams,TRUE)
    
    fittingResults <- fitWholeLevel(currentLevel, errorParams, TRUE)
    
  }
  
  model <- fittingResults$nodes
  
  return(model)
}

computeNodes <- function(samples, groups, features=1:11, children=0, errorParams=list(range=20, smooth=5, tol=10), maxErrorFactor=1.15, progress=FALSE){
  
  g <- unique(groups)
  G <- length(g)
  
  levels <- list()
  
  for(j in 1:G){
    
    node <- list()
    
    thisClassSamplesIndex <- which(groups == g[j])
    thisClassSamples <- list()
    
    for(i in features)
      thisClassSamples[[i]] <- samples[[i]][thisClassSamplesIndex,]
    
    node[["samples"]] <- thisClassSamples
    
    if(length(thisClassSamplesIndex) == 1){
      
      node[["representant"]] <- thisClassSamples
      node[["meanError"]] <- children[[thisClassSamplesIndex]]$meanError
      node[["maxError"]] <- children[[thisClassSamplesIndex]]$maxError
      node[["variation"]] <- children[[thisClassSamplesIndex]]$variation
      node[["errors"]] <- children[[thisClassSamplesIndex]]$errors
      node[["children"]] <- children[thisClassSamplesIndex]
    }
    else{
      
      node[["representant"]] <- list()
      node[["meanError"]] <- list()
      node[["variation"]] <- list()
      node[["errors"]] <- list()
      node[["maxError"]] <- list()
      
      for(i in features){
        
        meanClassSample <- colMeans(thisClassSamples[[i]])
        
        node[["representant"]][[i]] <- meanClassSample
        
        #compute the mean error and the mean error kind
        icpResults <- apply(thisClassSamples[[i]], 1, function(target, reference){
          
          return (do.call(my.dtwBasedDistance2, merge.list(list(reference, curveCorrection3(target, reference, 1)), errorParams)))
          
        }, meanClassSample)
        
        errors <- list2vector(getAllFieldFromList(icpResults, "error", 2))
        dists <- getAllFieldFromList(icpResults, "dist", 2)
        
        node[["meanError"]][[i]] <- mean(errors) #meanOfInterval(errors)
        node[["variation"]][[i]] <- sum(abs(errors - node[["meanError"]][[i]]))/length(errors) #sd(errors)
        node[["errors"]][[i]] <- errors
        node[["maxError"]][[i]] <- max(errors)
        node[["maxError"]][[i]] <- node[["maxError"]][[i]] * maxErrorFactor + node[["variation"]][[i]]
        
        if(is.list(children)){
          node[["children"]] <- children[thisClassSamplesIndex]
          
          #computes the mean max error of the children
          maxErrors <- getAllFieldFromList(node[["children"]], "maxError", 2)
          maxErrors <- getAllFieldFromList(maxErrors, i, 2)
          maxErrors <- list2vector(maxErrors)
          maxErrors <- mean(maxErrors)
          if(node[["maxError"]][[i]] < maxErrors)
            node[["maxError"]][[i]] <- maxErrors
        }
      }
    }
    
    levels[[g[j]]] <- node
    
    if(progress)
      cat("Computing node:", j*100/G, "%\n")
  }
  
  return(levels)
}

computeGrouping <- function(nodesAndErrorMatrix, nGroups=0, threshold=0, progress=FALSE){
  
  nodes <- nodesAndErrorMatrix$nodes
  errorMatrix <- nodesAndErrorMatrix$similarityMatrix
  
  C <- length(nodes)
  
  #determine the maximum number of groups for the level 1
  if(nGroups <= 0)
    nGroups <- floor(C/6)
  
  #determines the thresholding on the similarity index for grouping
  if(threshold <= 0)
    threshold <- C/nGroups
  
  #initiates the similarity matrix
  similarityMatrix <- matrix(rep(0, C*C), nrow=C)
  n <- computeNumberOfCombinations(C, 2)
  if(progress)
    cat(n, " combinations\n")
  m <- 1
  features <- length(errorMatrix)
  
  #computes the similarity matrix with the C representants
  for(j in 1:(C-1)){
    for(k in (j+1):C){
      
      data <- mapply(function(x){
        
        return(errorMatrix[[x]][j,k])
        
      }, 1:features)
      
      data <- data.frame(matrix(data, nrow=1))
      
      similarityMatrix[j,k] <- attr(predict(nodes[[j]][["classifier"]], data, decision.values=TRUE), "decision.values")
       
      similarityMatrix[k,j] <- similarityMatrix[j,k]
      
      if(progress){
        cat("computing similarity:", m*100/n, "%\n")
        m <- m + 1
      }
    }
  }
  
  #Transforming the values in the similarityMatrix so it will lie in the interval [0..inf]
  zeros <- which(similarityMatrix == 0)
  similarityMatrix[zeros] <- max(similarityMatrix) + 1
  similarityMatrix <- (similarityMatrix * -1) + min(similarityMatrix)
  
  #print(similarityMatrix)
  #computes the rank of similarity for the whole matrix
  similarityMatrix <- t(apply(similarityMatrix, 1, rank, ties.method = "random"))
  similarityMatrix <- similarityMatrix - 1
  #computes the similarity index instead of the ranking
  for(j in 1:(C-1)){
    for(k in (j+1):C){
      
      similarityMatrix[j,k] <- mean(c(similarityMatrix[j,k], similarityMatrix[k,j]))
      similarityMatrix[k,j] <- similarityMatrix[j,k]
    }
  }
  print(similarityMatrix)
  
  groups <- rep(0, C)
  
  #puts the first representant into the first group
  groups[1] <- 1
  groupIndex <- 2
  
  #exchanges all zeros by the greatest possible value
  similarityMatrix[which(similarityMatrix == 0)] <- C + 1
  
  for(j in 2:C){
    
    #while the jth node has no group, ...
    while(groups[j] == 0){
      #gets the closest representant
      closest <- which.min(similarityMatrix[j,])
      #gets the similarity index for the closest
      value <- similarityMatrix[j,closest]
      
      #if this value is smaller or equal than the threshold, ...
      if(value <= threshold){
        #and if closest already has a group, ...
        if(groups[closest] != 0){
          
          #assigns the group of the closest to it
          groups[j] <- groups[closest]
        }
        else{
          
          #if it is allowed to create another group, ...
          if(groupIndex <= nGroups){
            
            #creates a new group and assigns it to this representant
            groups[j] <- groupIndex
            #and to its closest
            groups[closest] <- groupIndex
            #updates group index
            groupIndex <- groupIndex + 1
          }
          else{
            
            similarityMatrix[j, closest] <- C + 1
          }
        }
      }
      else{
        
        #if it is allowed to create another group
        if(groupIndex <= nGroups){
          #creates a new group and assigns it to this representant
          groups[j] <- groupIndex
          #updates the group index
          groupIndex <- groupIndex + 1
        }
        else{
          
          #similarityMatrix[j, closest] <- C + 1
          similarityMatrix[j,] <- similarityMatrix[j,] - 1
        }
      }
    }
    
    if(progress)
      cat("choosing groups:", j*100/C, "%\n")
  }
  
  return(list(groups=groups, level=nodes))
}

fitWholeLevel <- function(nodes, errorParams=list(range=20, smooth=5, tol=10), progress=FALSE){
  
  C <- length(nodes)
  
  #initiates the similarity matrix
  similarityMatrix <- list()
  n <- computeNumberOfCombinations(C, 2)
  if(progress)
    cat(n, " combinations\n")
  m <- 1
  
  features <- length(nodes[[1]][["errors"]])
  
  #computes the similarity matrix with the C representants
  for(f in 1:features){
    
    similarityMatrix[[f]] <- matrix(rep(0, C*C), nrow=C)
    
    for(j in 1:(C-1)){
      for(k in (j+1):C){
        
        similarityMatrix[[f]][j,k] <- do.call(my.dtwBasedDistance2, merge.list(list(nodes[[j]][["representant"]][[f]],
                                                                        nodes[[k]][["representant"]][[f]]), errorParams))$error
        similarityMatrix[[f]][k,j] <- similarityMatrix[[f]][j,k]
        
        if(progress){
          cat("computing similarity:", m*100/n, "%\n")
          m <- m + 1
        }
      }
    }
  }
  
  for(i in 1:C){
    
    otherClasses <- c(1:C)[-i]
    
    otherErrors <- lapply(similarityMatrix, function(x, k, o){
      
      return(x[k,o])
      
    }, i, otherClasses)
    
    
    nodes[[i]] <- fitNodeOfTrees(nodes[[i]], otherErrors)
    
    if(progress)
      cat("Fitting models:", i*100/C, "%\n")
  }
  
  return(list(nodes=nodes, similarityMatrix=similarityMatrix))
}

fitNodeOfTrees <- function(node, otherErrors){
  
  #Gets the number of descriptors
  N <- length(node[["errors"]])
  
  #gets the number of samples of this class
  M <- length(node[["errors"]][[1]])
  
  #gets the number of other classes
  Mo <- length(otherErrors[[1]])
  
  data <- data.frame(matrix(ncol=N + 1, nrow=M + Mo, dimnames = list(c(), c(paste(1:N), "Y"))))
  
  for(i in 1:M){
    
    for(j in 1:N){
      
      data[i,j] <- node[["errors"]][[j]][i]
    }
    data[i,"Y"] <- "Yes"
  }
  
  for(i in 1:Mo){
    
    for(j in 1:N){
      
      data[i + M,j] <- otherErrors[[j]][i]
    }
    data[i + M,"Y"] <- "No"
  }
  
  #node[["classifier"]] <- boosting(Y~., data, FALSE, N, control = rpart.control(minsplit = 2, maxdepth = N))
  node[["classifier"]] <- rpart(Y~., data, control = rpart.control(minsplit = 2, maxdepth = N))
  #node[["classifier"]] <- svm(formula = Y~., data = data, kernel="radial", type = "C-classification")
  
  return(node)
}

hierarchicalFeatureBasedPrediction3 <- function(model, testDir="", testing=character(0), subset=integer(0), logFile="",
                                                errorParams=list(range=20, smooth=5, tol=10), isToPlot=FALSE){
  
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
  #the number of models, one for each descriptor
  N <- length(model[[1]][["errors"]])
  
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
  
  cat("Testing", M, "samples!\n", file=logFile, append=FALSE)
  
  #for each test sample, ...
  for(m in 1:M){
    
    cat("\npredicting test", m, ":", file=logFile, append=TRUE)
    
    start <- getTime()
    
    comparisons <- 0
    
    #initializes the votes as a matrix with zeros (these zeros will be ignored later)
    votes <- matrix(c(0,0), nrow=1)
    
    levelIndex <- list(c(0, 1))
    
    levelQueue <- list(model)
        
    #for each node from second level which matched the test, ...
    while(length(levelQueue) > 0){
      
      branch <- levelQueue[[1]]
      
      if(levelIndex[[1]][1] != 0)
        cat("\n", concatenate(rep("   ", levelIndex[[1]][1])), levelIndex[[1]][1], "level, node", levelIndex[[1]][2], file=logFile, append=TRUE)
      
      Nnodes <- length(branch)
      
      #gets the first level's representants
      #representants <- list2matrix(getAllFieldFromList(branch, "representant", 2))
      representants <- lapply(1:N, function(x, y, z){
        
        #return(list2matrix(y[(1+((x-1)*z)):(z*x)]))
        return(list2matrix(y[,x]))
        
      }, list2matrix(getAllFieldFromList(branch, "representant", 2)), Nnodes)
      
      targets <- 1
      errors <- data.frame(matrix(nrow=Nnodes, ncol=N))
      
      okError <- rep(TRUE, Nnodes)
      passed <- c(1:(Nnodes))
      
      maxErrors <- lapply(1:N, function(x, y, z){
        
        #return(list2matrix(y[(1+((x-1)*z)):(z*x)]))
        return(list2matrix(y[,x]))
        
      }, list2matrix(getAllFieldFromList(branch, "maxError", 2)), Nnodes)
      
      for(i in 1:N){
        
        test <- tests[[i]][m,]
        
        #if there is more than one representant at this level
        if(!is.null(dim(representants[[i]]))){
          #computes the errors for the first level
          icpResults <- apply(representants[[i]], 1, function(reference, target){
            
            return (do.call(my.dtwBasedDistance2, merge.list(list(reference, curveCorrection3(target, reference, 1)), errorParams)))
            #return (dtw(reference, curveCorrection3(target, reference, 1)))
            
          }, test)
          
          errors[,i] <- list2vector(getAllFieldFromList(icpResults, "error", 2))
          targets <- getAllFieldFromList(icpResults, "target", 2)
          
          
          #okError <- okError & (errors[,i] < maxErrors[[i]])
        }
        else{
          #computes the error with the single first level's representant
          icpResults <- do.call(my.dtwBasedDistance2, merge.list(list(representants[[i]], curveCorrection3(test, representants[[i]], 1)), errorParams))
          #icpResults <- dtw(representants, curveCorrection3(test, representants, 1))
          targets <- list(icpResults$target)
          errors[,i] <- icpResults$error
          #maxError <- branch[[1]]$maxError
          
          #okError <- okError & (errors[,i] < maxError[[i]])
        }
      }
      
      #classificationResult <- mapply(function(x){
        
        #return(attr(predict(branch[[x]][["classifier"]], errors[x,], decision.values=TRUE), "decision.values"))
        #return(median(as.numeric(errors[x,])))
        
      #}, 1:Nnodes)
      
      passed <- makeDecision(errors, 10)
      
      #passed <- which((classificationResult > 0) & okError)
      
      #rmByError <- which(!okError)
      #rmByClassf <- which(classificationResult < 0)
      
      #if(length(rmByError) > 0)
      #  cat(" [Removed by maxError: ", rmByError, "]", file = logFile, append=TRUE)
      
      #if(length(rmByClassf) > 0)
      #  cat(" [Removed by classifier: ", rmByClassf, "]", file=logFile, append=TRUE)
      
      comparisons <- comparisons + Nnodes
      
      #if this is a leaf, ...
      if(is.null(branch[[1]]$children)){
        
#         for(v in 1:(length(branch))){
#           if(isToPlot && length(passed) == 0 && names(branch)[v] == classes$fileClasses[m]){
#             
#             ymax = max(c(branch[[v]]$errorRange[,3], dists[[v]][,2]))
#             ymin = min(c(branch[[v]]$errorRange[,2], dists[[v]][,2]))
#             xmax = max(branch[[v]]$errorRange[,1], dists[[v]][,1])
#             xmin = min(branch[[v]]$errorRange[,1], dists[[v]][,1])
#             
#             #plots the error range graph
#             plot(c(1,xmax), c(ymin, ymax), col="white", main=concatenate(c("Error Range predictor ", f, ", ", i)))
#             lines(c(1,xmax), c(0,0), col="gray")
#             sorted <- sort.int(branch[[v]]$errorRange[,1], index.return = TRUE)$ix
#             meanRange = rowMeans(branch[[v]]$errorRange[(sorted),-1])
#             lines(x = branch[[v]]$errorRange[sorted,1], y = meanRange, col="green")
#             lines(x = branch[[v]]$errorRange[sorted,1], y = branch[[v]]$errorRange[sorted,2], col="red")
#             lines(x = branch[[v]]$errorRange[sorted,1], y = branch[[v]]$errorRange[sorted,3], col="red")
#             lines(dists[[v]], col="black")
#             
#             plot(branch[[v]]$representant, type="l", col="red", main=concatenate(c("Curves predictor ", f, ", ", i)))
#             lines(test, col="blue")
#             lines(x = targets[[v]][,1], y = targets[[v]][,2], col="black")
#             if(f == 1){
#               points(refPoints[[v]], branch[[v]]$representant[refPoints[[v]]], col="red")
#               points(targets[[v]][tarPoints[[v]],], col="black")
#             }
#           }
#         }
        
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
          #minError <- classificationResult[v]
          minError <- median(as.numeric(errors[passed,]))
          
          #cat(" -------", minErrorIndex, "------", file=logFile, append=TRUE)
          cat(" E =", minError, "\n", file=logFile, append=TRUE)
          
          #adds a vote for this leaf's class with the weight as the minimum error value
          #cat("leaf:", v, " descriptor:", i, "test:", m, "first level:", k, "second level:", j, "\n")
          votes <- rbind(votes, matrix(c(as.numeric(names(branch)[v]), minError), nrow=1))
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
    
    #removes the zeros from the votes matrix
    votes <- votes[-1,]
    
    if(is.null(dim(votes)))
      dim(votes) <- c(1,2)
    
    cat("\nvotes:", length(votes[,1]), "\n", file=logFile, append=TRUE)
    cat.matrix(votes, file=logFile, append=TRUE)

    cat("comparisons:", comparisons, "\n", file=logFile, append=TRUE)

    cat("test", m, ". ", file=logFile, append=TRUE)
    
    if(length(votes) > 1){
      
      result <- votes[which.min(votes[,2]),]
      #result <- votes[which.max(votes[,2]),]
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
}

print.hierarchicalModel <- function(model, file=""){
  
  N <- length(model)
  
  cat("", file=file, append=FALSE)
  print.hierarchicalModelAux(model, 1, file)
}

print.hierarchicalModelAux <- function(model, level, file){
  
  N <- length(model)
  
  for(i in 1:N){
    
    if(!is.null(names(model)[i]))
      cat("\n", concatenate(rep("\t", level)), ">", names(model)[i], file=file, append=TRUE)
    else
      cat("\n", concatenate(rep("\t", level)), ">", i, file=file, append=TRUE)
    
    if(!is.null(model[[i]]$children))
      print.hierarchicalModelAux(model[[i]]$children, level + 1, file)
    
    #cat("\n", concatenate(rep("\t", level)), ")", file=file, append=TRUE)
  }
}

makeDecision <- function(results, K=2){
  
  medians <- apply(results, 1, median)
  
  sortedIndex <- sort(medians, index.return = TRUE, decreasing = FALSE)[["ix"]][1:K]
  
  N <- length(results[1,])
  
  Y <- mapply(function(x){
    
    return(which.min(results[sortedIndex,x]))
    
  }, 1:N)
  
  counts <- mapply(function(x){
    
    return(length(which(Y == x)))
    
  }, 1:K)
  
  decision = sortedIndex[which.max(counts)]
  
  return(decision)
}
