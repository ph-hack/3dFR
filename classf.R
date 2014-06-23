library(class)
library(mmand)
library(e1071)

separateDataSet <- function(dataDir, sulfix=".unholed", proportions=c(0.7,0.3)){
  
  faceNames <- unique(getFaceID(dir(dataDir)))
  n <- length(faceNames)
  faceClasses <- getPersonID(faceNames)
  classes <- unique(faceClasses)
  classesNumbers <- tabulate(match(faceClasses, classes))
  
  singles <- which(classesNumbers == 1)
  singleClasses <- classes[singles]
  classes <- classes[-singles]
  classesNumbers <- classesNumbers[-singles]
  nClasses <- length(classes)
  
  training <- list()
  trainingIndex <- 1
  test <- list()
  testIndex <- 1
  
  for(i in 1:nClasses){
    
    trainingNumber <- round(proportions[1] * classesNumbers[i])
    trainingSample <- sample(1:classesNumbers[i], trainingNumber)
    testSample <- c(1:classesNumbers[i])[-trainingSample]
    
    for(j in 1:trainingNumber){
      
      training[[trainingIndex]] <- faceNames[which(faceClasses == classes[i])[trainingSample[j]]]
      trainingIndex <- trainingIndex + 1
    }
    
    for(j in 1:length(testSample)){
      
      test[[testIndex]] <- faceNames[which(faceClasses == classes[i])[testSample[j]]]
      testIndex <- testIndex + 1
    }
  }
  
  training <- list2vector(training)
  training <- concatenate(list(dataDir, training, sulfix, ".jpg.dat"))
  test <- list2vector(test)
  test <- concatenate(list(dataDir, test, sulfix, ".jpg.dat"))
  
  nSingle <- length(singleClasses)
  singles <- rep(0, nSingle)
  for(i in 1:nSingle){
    
    singles[i] <- concatenate(c(dataDir, faceNames[which(faceClasses == singleClasses[i])], sulfix, ".jpg.dat"))
  }
  
  (list(training=training, test=test, singles=singles))
}

# Trains a set of 'C' SVM models where C = number of classes ########################################################
# input:
#   imgsDir = the folder of the training vectors, which must be jpg files
#   weighted = wheather the weight to each pair of classes will be proportional to the number of samples of each class
#   kernel = a string specifying the kernel, e.g.= {radial, polynomial, sigmoid, linear}
#   tol = the tolerance
#   g = the gamma parameter setup, e.g. = {1, 2}
#       if the 'g' = 1, gamma = 1/m*dev^2; otherwise, gamma = 1/m
#   type = a string describing the type of data ('image', 'mainLines', 'icp')
#   maxIter = only for 'icp' type
# output:
#   a list with the set of SVM and a vector with the classes' names
svmTraining <- function(imgsDir, weighted=FALSE, kernel="radial", tol=10^(-3), g=0, type="image"){
  
  imgs <- dir(imgsDir)
  if(type == "image")
    imgs <- imgs[-which(my.strsplit(imgs, "[.]", 0) == "jpg")]
    
  m <- length(imgs)
  
  trainingImg <- 0
  if(type == "image")
    trainingImg <- readImageData(paste(imgsDir, imgs[1], sep=""))
  
  else if(type == "mainLines")
    trainingImg <- readMainLines(paste(imgsDir, imgs[1], sep=""))
  
  trainingMatrix <- matrix(trainingImg, nrow=1)
  
  imgsClasses <- rep(0, m)
  imgsClasses[1] <- getPersonID(imgs[1])
  
  for(j in 2:m){
    
    if(type == "image")
      trainingImg <- readImageData(paste(imgsDir, imgs[j], sep=""))
    
    else if(type == "mainLines")
      trainingImg <- readMainLines(paste(imgsDir, imgs[j], sep=""))
    
    #appends it to the matrix
    trainingMatrix <- rbind(trainingMatrix, matrix(trainingImg, nrow=1))
    #gets the classes' IDs
    imgsClasses[j] <- getPersonID(imgs[j])
    
    cat("building training matrix: ", (j*100/m), "%\n")
  }
  
  classes <- list()
  i <- 0
  while(length(imgsClasses > 0)){
    i <- i + 1
    ci <- imgsClasses[1]
    imgsClasses <- imgsClasses[-which(imgsClasses == ci)]
    classes[[i]] <- ci
  }
  
  imgsClasses <- getPersonID(imgs)
  
  mc <- length(classes) #the number of classes
  
  cat("Number of classes: ", mc, "\n")
  
  models <- list()
  n <- length(trainingMatrix[1,])
  
  #sets the gamma parameter
  #g <- 0
  if(g == 1){
    dev <- sapply(as.data.frame(trainingMatrix), sd)
    g <- 1/(m*dev^2)  
  }
  else{
    g <- 1/m 
  }
  
  for(i in 1:mc){
    cMatrix <- rep("-1", m)
    #sets the positive samples
    cMatrix[which(imgsClasses == classes[i])] <- "1"
    cMatrix <- list2vector(cMatrix)
    #computes the number of samples to the ith class
    mci <- length(which(imgsClasses == classes[i]))
    
    if(weighted){
      cMatrix <- factor(cMatrix)
      #sets the proportinal weights
      weights <- c("-1" = (m/(m - mci)), "1" = (m/mci))
      models[[i]] <- svm(trainingMatrix, cMatrix, cachesize=((n*m)/1024), scale=FALSE, gamma=g, class.weights=weights, kernel=kernel, tolerance=tol)
    }
    else{
      models[[i]] <- svm(trainingMatrix, factor(cMatrix), cachesize=((n*m)/1024), kernel=kernel, gamma=g, tolerance=tol)
    }
    
    cat("Training the SVMs: ", (i*100/mc), "%\n")
  }
  (list(models = models,
        classes = classes))
}

# Tests the classification of the SVM classifier #################################################################
# input:
#   testDir = the folder with the test vectors
#   models = a list with a list of the SVM models and a vector with the classes' names
#   t = (optional) a threshold
#   type = (optional) the mode of the operation
#          if 't' = 0, then 'type' = 1
# output:
#   a list with the accuracy and the SVM models
svmClassf <- function(testDir, models, t="a", type=1, typeOfInput="image"){
  
  tests <- dir(testDir)
  if(typeOfInput == "image")
    tests <- tests[-which(my.strsplit(tests, "[.]", 0) == "jpg")]
  
  M <- length(tests)
  
  if(typeOfInput == "image")
    testImg <- readImageData(paste(testDir, tests[1], sep=""))
  
  else if(typeOfInput == "mainLines")
    testImg <- readMainLines(paste(testDir, tests[1], sep=""))
  
  testMatrix <- matrix(testImg, nrow=1)
  
  for(j in 2:M){
    
    if(typeOfInput == "image")
      testImg <- readImageData(paste(testDir, tests[j], sep=""))
    
    else if(typeOfInput == "mainLines")
      testImg <- readMainLines(paste(testDir, tests[j], sep=""))
    
    #appends it to the matrix
    testMatrix <- rbind(testMatrix, matrix(testImg, nrow=1))
    
    cat("building test matrix: ", (j*100/M), "%\n")
  }
  
  corrects <- 0
  mc <- length(models[[1]])
  
  for(i in 1:M){
    values <- rep(0, mc)
    
    for(j in 1:mc){
      value <- attr(predict(models[[1]][[j]], matrix(testMatrix[i,], nrow=1), decision.values=TRUE), "decision.values")
      posClass <- strsplit(dimnames(value)[[2]][1], "[/]")[[1]][1]
      negClass <- strsplit(dimnames(value)[[2]][1], "[/]")[[1]][2]
      
      if(negClass == "-1"){ #if the positive class is the persons,
        #we want the value
        values[j] <- value[1]
      }
      else{ #if the negative class is the persons
        #we want the inverse value
        values[j] <- -value[1]
      }
    }
    
    k <- which.max(values)
    
    if(is.character(t)){ #if the threshold isn't setup
      #just check the greatest value
      if(type == 1){
        if((values[k] >= 0) && (models[[2]][k][[1]] == getPersonID(tests[i]))){
          corrects <- corrects + 1
          cat("- found \t")
        }
        else
          cat("- missed\t")
      }
      else{
        if(models[[2]][k][[1]] == getPersonID(tests[i])){
          corrects <- corrects + 1
          cat("- found \t")
        }
        else
          cat("- missed\t")
      }
    }
    else{
      if(values[k] >= t){ #checks the threshold
        if(models[[2]][k][[1]] == getPersonID(tests[i])){
          corrects <- corrects + 1
          cat("- found \t")
        }
        else
          cat("- missed\t")
      }
      else{
        #if the test vector class ins't in the training set
        if(length(models[[2]][which(models[[2]] == getPersonID(tests[i]))]) == 0){
          corrects <- corrects + 1
          cat("- unknown\t")
        }
        else
          cat("- missed\t")
      }
    }
    
    cat(models[[2]][k][[1]], "==", getPersonID(tests[i]), "\tdist = ", values[k], "\t")
      
    cat(i*100/M, "%\n")
  }
  #returns the correctness of the classification and the SVM models
  (list(acc = corrects*100/M,
        models = models[[1]]))
}

# tests the classification of vectors with k-neighbours and the specified distance ###########################################################
# it needs the directory of the training images ('trainingDir')
#, the directory of the test images ('testDir')
#, the distance method ('method') e.g. = {cosine, dice, avg, manhattan, euclidean}
#, the amount of neighbours to be considered ('top'), default = 1
# and the inferior threshold
# when top = 1, is the same as just looking for the closest training vector
# for top > 1, only the cosine and dice distances can be used
# the verification is made with the 1st part (numeric) of the filename
# the training and test sets must be composed of vectors, not images
# it returns the percentage of correctly classified images
kNeigbourClassf <- function(trainingDir, testDir, method="euclidean", top=1, t=0.41, type="image"){
  
  training <- dir(trainingDir)
  if(type == "image")
    training <- training[-which(my.strsplit(training, "[.]", 0) == "jpg")]
  
  test <- dir(testDir)
  if(type == "image")
    test <- test[-which(my.strsplit(test, "[.]", 0) == "jpg")]
  
  m <- length(training)
  M <- length(test)
  corrects <- 0
  
  for(i in 1:M){
    
    #reads the ith test image
    if(type == "image")
      testImg <- readImageData(paste(testDir, test[i], sep=""))
    
    else if(type == "mainLines")
      testImg <- readMainLines(paste(testDir, test[i], sep=""))
    
    dists <- rep(0, m) #initializes the vector with the distances
    
    for(j in 1:m){ #for each training image
      
      #reads the jth training image
      if(type == "image")
        trainingImg <- readImageData(paste(trainingDir, training[j], sep=""))
      
      else if(type == "mainLines")
        trainingImg <- readMainLines(paste(trainingDir, training[j], sep=""))
      
      #computes the distance between the ith test sample and the jth training sample
      if(method == "cosine"){
        dists[j] = sum(testImg * trainingImg)/sqrt(sum(testImg^2) * sum(trainingImg^2))
      }
      else{
        if(method == "harmonic"){
          dists[j] = abs(sum(testImg) - (2 * sum(testImg * trainingImg/(testImg + trainingImg))))
        }
        else{
          if(method == "dice"){
            dists[j] = (2 * sum(testImg * trainingImg))/(sum(testImg^2) + sum(trainingImg^2))
          }
          else{
            if(method == "avg"){
              dists[j] = sum(abs(testImg - trainingImg) + max(abs(testImg - trainingImg)))/2
            }
            else{
              distMatrix <- matrix(c(testImg, trainingImg), nrow=2, byrow=TRUE)
              dists[j] = dist(distMatrix, method=method)[1]
            }
          }
        }
      }
    }
    
    if(top == 1){ #if it to just look for the closests vector
      #takes the index of the closest image vector
      k <- which.min(dists)
      if((method == "cosine") || (method == "dice")){
        k <- which.max(dists)
      }
      
      #if the classifier is correct
      if(dists[k] >= t){ #for the cosine distance only
        if(getPersonID(test[i]) == getPersonID(training[k])){
          corrects <- corrects + 1
          cat("- found \t")
        }
        else
          cat("- missed\t")
      }
      else{
        #if the test vector's class isn't in the training set
        if(length(getPersonID(training)[which(getPersonID(training) == getPersonID(test[i]))]) == 0){
          corrects <- corrects + 1
          cat("- unknown\t")
        }
        else
          cat("- missed\t")
      }
      cat(getPersonID(test[i]), "==", getPersonID(training[k]), "\tdist = ", dists[k], "\t")
    }
    else{
      k <- rep(0, top)
      classes <- list()
      cDists <- list()
      c <- 0
      #computes the acumulative values to each of the 'top' closest classes
      for(j in 1:top){
        k[j] <- which.max(dists)
        if(length(classes[which(classes == getPersonID(training[k[j]]))]) == 0){
          c <- c + 1
          classes[[c]] <- getPersonID(training[k[j]])
          cDists[[c]] <- dists[k[j]]
        }
        else{
          cDists[[c]] <- cDists[[c]] + dists[k[j]]
        }
        dists[k[j]] <- -1
      }
      
      c <- which.max(cDists)
      
      #if the classifier is correct
      if(cDists[[c]] >= t){
        if(getPersonID(test[i]) == classes[[c]]){
          corrects <- corrects + 1
          cat("- found \t")
        }
      }
      else{
        #if the test vector isn't in the training set
        if(length(getPersonID(training)[which(getPersonID(training) == getPersonID(test[i]))]) == 0){
          corrects <- corrects + 1
          cat("- unknown\t")
        }
      }
      cat(getPersonID(test[i]), "==", classes[[c]], "\tdist = ", cDists[[c]], "\t")
    }
    
    cat(i*100/M, "%\n") #prints the progress
  }
  (corrects*100/M) #returns the percentage of correct classifications
}

# tests the classification of vectors with k-neighbours and the specified distance ###########################################################
# it needs the directory of the training images ('trainingDir')
#, the directory of the test images ('testDir')
#, the distance method ('method') e.g. = {cosine, dice, avg, manhattan, euclidean}
#, the amount of neighbours to be considered ('top'), default = 1
# and the inferior threshold
# when top = 1, is the same as just looking for the closest training vector
# for top > 1, only the cosine and dice distances can be used
# the verification is made with the 1st part (numeric) of the filename
# the training and test sets must be composed of vectors, not images
# it returns the percentage of correctly classified images
ICPClassf <- function(trainingDir, testDir, closest="", method="election", by="error", weights=0, t=0.5, maxIter=10, minIter=5, pSample=0.33, 
                      range=0, nClosest=0, outlierCorrection=FALSE, smooth=0, meanDir="", ldmkDir="", logFile="", append=FALSE, exceptClasses=0){
  
  training <- trainingDir
  if(length(trainingDir) == 1)
    training <- dir(trainingDir)
  else
    trainingDir <- ""
  
  test <- testDir
  if(!is.file(testDir))
    test <- dir(testDir)
  else
    testDir <- ""
  
  if(length(range) == 2)
    test <- test[range[1]:range[2]]
  
  #takes out the unwanted classes given by 'exceptClasses'
  if(exceptClasses[1] != 0){
    
    for(i in 1:(length(exceptClasses))){
    
      training <- training[-which(getPersonID(training) == exceptClasses[i])]
      test <- test[-which(getPersonID(test) == exceptClasses[i])]
    }
  }
  
  classTraining <- unique(getPersonID(training))
  
  m <- length(training)
  M <- length(test)
  corrects <- 0
  
  trC <- length(classTraining)
  confusion <- matrix(rep(0, M*trC), ncol=trC)
  
  cat("Classifying ", M, " faces...\n", file=logFile, append=append)
  
  for(i in 1:M){
    
    begin <- getTime()
    
    #reads the ith test image
    testImg <- readMainLines(paste(testDir, test[i], sep=""), "list")
    n <- length(testImg)
    
    #cat("Testing ", test[i], "\n", file=logFile, append=TRUE)
    
    if(closest != ""){
      
      training <- readLines(concatenate(c(closest, "cl__", getFaceID(test[i]), ".txt")))
      training <- concatenate(list(training, ".lines"))
      if(nClosest != 0)
        training <- training[1:nClosest]
      m <- length(training)
      #cat("m = ", m, "\n", trainingDir, "\n")
    }
    
    #initializes the matrix with the distances, one column for each line
    dists <- matrix(rep(0, length(testImg) * m), nrow=m)
    
    start <- getTime()
    
    for(j in 1:m){ #for each training image
      
      #cat("with the training ", training[j], "\n", file=logFile, append=TRUE)
      
      #reads the jth training image
      trainingImg <- readMainLines(paste(trainingDir, training[j], sep=""), "list")
      
      if(method != "cosine"){
        #computes the distance between the ith test sample and the jth training sample
        if(by == "error"){
          
          #initializes a vector to manage the order of the curves (descriptor vectors)
          sortedIndex <- 1:n
          #if the weights are given, then the order will be changed for each training sample
          if(length(weights) > 1){
            
            #ws <- discretize(weights[[getPersonID(training[j])]], 11, c(0,1.1))
            ws <- weights[[getPersonID(training[j])]]
            maxW <- max(ws) + 1
            
            #cat(getPersonID(training[j]), "\n")
            
            for(k in 1:n){
              #print(ws)
              if(length(which.min(ws)) == 0){
                
                cat("FUCK!\n", ws, "\n", weights[[getPersonID(training[j])]], "\n", training[j], "\n")
                return()
              }
              sortedIndex[k] <- which.min(ws)
              ws[sortedIndex[k]] <- maxW
            }
            #cat("done\n")
          }
          #cat("\n")
          
          for(k in 1:n){
            K <- sortedIndex[k]
            
            if(outlierCorrection){
              
              #cat("getting lines of mean image\n")
              meanImg <- getMainLines(readImageData(concatenate(c(meanDir, getPersonID(training[j]), "d000.mean.jpg.dat"))),
                                      readLandmark(concatenate(c(ldmkDir, getFaceID(training[j]), ".ldmk"))))
              #cat("applying my ICP\n")
              dists[j, k] <- my.icp.2d.v2(curveCorrection3(testImg[[K]], meanImg[[K]], smooth), curveCorrection3(trainingImg[[K]], meanImg[[K]], smooth), maxIter, minIter, pSample)$error 
              #cat("Done\n")
            }
            else
              dists[j, k] <- my.icp.2d.v2(testImg[[K]], trainingImg[[K]], maxIter, minIter, pSample)$error
          }
          
        }
        else if(by == "energy"){
          
          for(k in 1:n)
            dists[j, k] <- my.icp.2d.v2(testImg[[k]], trainingImg[[k]], maxIter, minIter, pSample)$energyMean
          
        }
      }
      else{
        
        tr <- my.icp.2d.v2(testImg[[1]], trainingImg[[1]], maxIter, minIter, pSample)$target[,2]
        for(k in 2:n)
          tr <- c(tr, my.icp.2d.v2(testImg[[k]], trainingImg[[k]], maxIter, minIter, pSample)$target[,2])
        
        te <- testImg[[1]]
        for(k in 2:n)
          te <- c(te, testImg[[k]])
                
        dists[j,1] = sum(te * tr)/sqrt(sum(te^2) * sum(tr^2))
      }
      
#       if(length(weights) > 1)
#         dists[j,] <- dists[j,] + weights[[getPersonID(training[j])]]
      
#       for(k in 1:n)
#         cat("line ", k, "\t", dists[j, k], "\n")
    }
    
    cat("Error measurement: ", crono.end(start), "\n", file=logFile, append=TRUE)
    
    start <- getTime()
    
    #takes the index of the closest image vector
    k <- 0
    value <- 0
    
    if(method == "election"){
      cat("i:", i, "\n", file=logFile, append=TRUE)
      votes <- matrix(rep(0, 2*n), ncol=2)
      
      values <- rep(0, n)
      if(length(weights) > 1)
        values <- 0:10/10
        
      for(k in 1:n){
        #d <- 0
        #if(icpBy == "mean")
          d <- which.min(dists[,k])
        #else
          #d <- which.max(dists[,k])
        #cat("training: ", training[d], "\n")
        
        tr <- which(classTraining == getPersonID(training[d]))
        confusion[i,tr] <- confusion[i,tr] + 1
        
        #votes[k,1] <- which(getPersonID(training) == getPersonID(training[d]))[1]
        #votes[k,1] <- d
	      votes[k,1] <- which(classTraining == getPersonID(training[d]))
        votes[k,2] <- dists[d,k] + values[k]
#         if(length(weights) > 1){
#           ws <- discretize(weights[[getPersonID(training[d])]], 11, c(0,1.1))
#           votes[k,2] <- votes[k,2] + ws[k]
#           
#           cat("--- vote in line ", k, ":\t", classTraining[votes[k,1]], ",", tr, "\twith ", votes[k,2], "\t weigth:", ws[k], "\t", file=logFile, append=TRUE)
#         }
#         else{
          
          cat("--- vote in line ", k, ":\t", classTraining[votes[k,1]], ",", tr, "\twith ", votes[k,2], "\t", file=logFile, append=TRUE)
          #cat("<icpBy = ", icpBy, "\n")
#         }
        dists[d,k] <- 100
        d <- which.min(dists[,k])
        cat(getPersonID(training[d]), ", ", file=logFile, append=TRUE)
        dists[d,k] <- 100
        d <- which.min(dists[,k])
        cat(getPersonID(training[d]), ", ", file=logFile, append=TRUE)
        dists[d,k] <- 100
        d <- which.min(dists[,k])
        cat(getPersonID(training[d]), ", ", file=logFile, append=TRUE)
        dists[d,k] <- 100
        d <- which.min(dists[,k])
        cat(getPersonID(training[d]), ", ", file=logFile, append=TRUE)
        dists[d,k] <- 100
        d <- which.min(dists[,k])
        cat(getPersonID(training[d]), "\n", file=logFile, append=TRUE)
      }
      electionBy <- "min"
      #if(icpBy == "cosine")
        #electionBy <- "max"
      election <- ponderateVote(votes, electionBy)
      k <- election[1]
      value <- 1/election[2]
      
#       minDists <- rep(0, n)
#       for(k in 1:n)
#         minDists[k] <- which.min(dists[,k])
#       
#       k <- getMode(minDists)
#       
#       if(is.integer(k)){
#         value <- dists[k, which.min(dists[k,])]
#         cat("by vote    \t")
#       }
#     
#       else{
#         
#         less <- min(dists[minDists,])
#         #less <- min(dists[which.min(dists[,1]),1], dists[which.min(dists[,2]),2], dists[which.min(dists[,3]),3], dists[which.min(dists[,4]),4])
#         
#         for(l in 1:n)
#           if(length(which(dists[,l] == less)) > 0){
#             k <- which(dists[,l] == less)
#             value <- dists[k,l]
#           }
#         
#         cat("by the less\t")
#       }
    }
    else if(method == "mean"){
      
      sumOfDists <- 0
      for(l in 1:n){
        sumOfDists <- sumOfDists + dists[,l]
      }
      
#       for(j in 1:m)
#         cat("sum of dists for training ", training[j], ":\t", sumOfDists[j], "\n")
      
      k <- which.min(sumOfDists)
      
#       sumOfDists <- 0
#       for(l in 1:n)
#         sumOfDists <- dists[k,l] 
      
      value <- sumOfDists[k]/n
      #cat("by the mean ")
    }
    else if(method == "cosine"){
      
      k <- which.max(dists[,1])
      value <- dists[k,1]
      t=1
      #cat("by the cosine ")
    }
    
    #if the classifier is correct
    if(value <= t){
      if(getPersonID(test[i]) == classTraining[k]){
        corrects <- corrects + 1
        cat("- found  \t", file=logFile, append=TRUE)
      }
      else
        cat("- missed \t", file=logFile, append=TRUE)
    }
    else{
      #if the test vector's class isn't in the training set
      if(length(classTraining[which(classTraining == getPersonID(test[i]))]) == 0){
        corrects <- corrects + 1
        cat("- unknown\t", file=logFile, append=TRUE)
      }
      else
        cat("- missed \t", file=logFile, append=TRUE)
    }
    cat(getPersonID(test[i]), "==", classTraining[k], "\tdist = ", value, file=logFile, append=TRUE)
    
    cat("\t", i*100/M, "%\n", file=logFile, append=TRUE) #prints the progress
    
    cat("Classification finished: ", crono.end(start), "\n", file=logFile, append=TRUE)
    cat("Total ", crono.end(begin), "\n", file=logFile, append=TRUE)
  }
  (list(accuracy=corrects*100/M, corrects=corrects, confusion=confusion)) #returns the percentage of correct classifications
}

kNeigbourSelector <- function(trainingDir, testDir, training=0, test=0, toFile="", method="euclidean", amount=200, range=0, logFile="", useFisher=TRUE){
  
  if(length(training) == 1){
    training <- dir(trainingDir)
    if(!useFisher)
      training <- training[-which(my.strsplit(training, "[.]", 0) == "jpg")]
  }
  else{
    
    training <- getFaceID(training)
    training <- concatenate(list(training, ".jpg.dat"))
  }
  
  if(length(test) == 1){
    test <- dir(testDir)
    if(!useFisher)
      test <- test[-which(my.strsplit(test, "[.]", 0) == "jpg")]
  }
  else{
    
    test <- getFaceID(test)
    test <- concatenate(list(test, ".jpg.dat"))
  }
  
  if(length(range) == 2)
    test <- test[range[1]:range[2]]
  
  m <- length(training)
  M <- length(test)
  closest <- list()
  
  cat("Analysing ", M, " faces\n", file=logFile, append=FALSE)
  
  for(i in 1:M){
    
    start <- getTime()
    
    #reads the ith test image
    testImg <- 0
    if(useFisher)
      testImg <- scan(paste(testDir, test[j], sep=""), quiet=TRUE)
    else
      testImg <- readImageData(paste(testDir, test[i], sep=""))
    
    dists <- rep(0, m) #initializes the vector with the distances
    
    for(j in 1:m){ #for each training image
      
      #reads the jth training image
      trainingImg <- 0
      if(useFisher)
        trainingImg <- scan(paste(trainingDir, training[j], sep=""), quiet=TRUE)
      else
        trainingImg <- readImageData(paste(trainingDir, training[j], sep=""))
      
      #computes the distance between the ith test sample and the jth training sample
      if(method == "cosine"){
        dists[j] = sum(testImg * trainingImg)/sqrt(sum(testImg^2) * sum(trainingImg^2))
      }
      else{
        if(method == "harmonic"){
          dists[j] = abs(sum(testImg) - (2 * sum(testImg * trainingImg/(testImg + trainingImg))))
        }
        else{
          if(method == "dice"){
            dists[j] = (2 * sum(testImg * trainingImg))/(sum(testImg^2) + sum(trainingImg^2))
          }
          else{
            if(method == "avg"){
              dists[j] = sum(abs(testImg - trainingImg) + max(abs(testImg - trainingImg)))/2
            }
            else{
              distMatrix <- matrix(c(testImg, trainingImg), nrow=2, byrow=TRUE)
              dists[j] = dist(distMatrix, method=method)[1]
            }
          }
        }
      }
    }
    
    cl <- rep("", amount)
    
    for(j in 1:amount){ #if it to just look for the closests vector
      #takes the index of the closest image vector
      k <- which.min(dists)
      if((method == "cosine") || (method == "dice")){
        k <- which.max(dists)
      }
      
      cl[j] <- getFaceID(training[k])
      
      if((method == "cosine") || (method == "dice"))
        dists[k] <- -100
      else
        dists[k] <- 100
    }
    
    closest[[i]] <- cl
    
    if(toFile != "")
      write(cl, concatenate(c(toFile, "_", getFaceID(test[i]), ".txt")), 1)
    
    cat(i*100/M, "%\n") #prints the progress
    cat("Testing ", test[i], ": ", crono.end(start), ", ", i*100/M, "%\n", file=logFile, append=TRUE)
  }
  (closest) #returns the 'amount' closest images
}

is.file <- function(string){
  
  chars <- strsplit(string, "")[[1]]
  
  if(chars[length(chars)] == "/")
    (FALSE)
  else
    (TRUE)
}
