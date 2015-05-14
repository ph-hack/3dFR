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
  
  currentLevel <- fitNodeTrees(currentLevel)
  
  nGroups <- length(groupNumbers)
  
  for(j in 1:nGroups){
    
    #divides the leafs into groups
    gResult <- computeGrouping(currentLevel, groupNumbers[j], errorParams = errorParams, progress=TRUE)
    currentLevel <- gResult$level
    #mounts the first level
    currentLevel <- computeNodes(list2matrix(getAllFieldFromList(currentLevel, "representant", 2)),
                                 gResult$groups, currentLevel, errorParams,TRUE)
    
    currentLevel <- fitNodeTrees(currentLevel)
  }
  
  model <- currentLevel
  
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
          maxErrors <- getAllFieldFromList(children, "maxError", 3)
          maxErrors <- list2vector(maxErrors[[i]])
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

computeGrouping <- function(nodes, nGroups=0, threshold=0, errorParams=list(range=20, smooth=5, tol=10), progress=FALSE){
  
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
  #computes the similarity matrix with the C representants
  for(j in 1:(C-1)){
    for(k in (j+1):C){
      
      similarityMatrix[j,k] <- do.call(errorFunction, merge.list(list(nodes[[j]][["representant"]],
                                                                      nodes[[k]][["representant"]]), errorParams))$error
      similarityMatrix[k,j] <- similarityMatrix[j,k]
      
      if(progress){
        cat("computing similarity:", m*100/n, "%\n")
        m <- m + 1
      }
    }
  }
  
  nodes <- computeFamiliarityWeight(nodes, similarityMatrix)
  
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

fitNodeTrees <- function(node, otherErrors){
  
  #Gets the number of descriptors
  N <- length(node[["errors"]])
  
  #gets the number of samples of this class
  M <- length(node[["errors"]][[1]])
  
  #gets the number of other classes
  Mo <- length(otherErrors)
  
  data <- data.frame(matrix(ncol=N + 1, nrow=M + Mo, dimnames = list(c(), c(paste(1:N), "Y"))))
  
  for(i in 1:M){
    
    for(j in 1:N){
      
      data[i,j] <- node[["errors"]][[j]][i]
    }
    data[i,"Y"] <- 1
  }
  
  for(i in 1:Mo){
    
    for(j in 1:N){
      
      data[i + M,j] <- otherErrors[[j]][i]
    }
    data[i + M,"Y"] <- 0
  }
  
  node[["classifier"]] <- boosting(Y ~., data, FALSE, N, control = rpart.control(minsplit = 2, maxdepth = N))
  
  return(node)
}