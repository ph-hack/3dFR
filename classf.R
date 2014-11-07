library(class)
library(mmand)
library(e1071)

# Retrieves the face's person ID from the fileName, which can be the ABS, JPG or LDMK file(s).
# input:
#   fileName = either a string with the name of a 3D face file(ABS, JPG or LDMK)
#              or a vector with all file names
# output:
#   either a string containing only the ID of the face or a vector with all IDs.
#   e.g. "02463d452.abs" -> "02463"
getPersonID <- function(fileName){
  
  name <- ""
  n <- length(fileName)
  if(n > 1){
    
    name <- rep(0, n)
    for(i in 1:n){
      
      aux <- strsplit(fileName[i], "[/.]")[[1]]
      name[i] <- aux[which(regexpr(text=aux, pattern="[0-9]+d[0-9]+") == 1)]
      name[i] <- strsplit(name[i], "d")[[1]][1]
    }
  }
  else{
    name <- strsplit(fileName, "[/.]")[[1]]
    name <- name[which(regexpr(text=name, pattern="[0-9]+d[0-9]+") == 1)]
    name <- strsplit(name, "d")[[1]][1]
  }
  
  (name)
}

getFaceID <- function(fileName){
  
  n <- length(fileName)
  name <- ""
  
  if(n > 1){
    
    name <- rep("", n)
    
    for(i in 1:n){
      
      aux <- strsplit(fileName[i], "[/.]")[[1]]
      name[i] <- aux[which(regexpr(text=aux, pattern="[0-9]+d[0-9]+") == 1)]
    }
  }
  else{
    name <- strsplit(fileName, "[/.]")[[1]]
    name <- name[which(regexpr(text=name, pattern="[0-9]+d[0-9]+") == 1)]
  }
  
  (name)
}

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
    trainingDir <- getDirectory(training[1])
  
  test <- testDir
  if(length(testDir) == 1)
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
      
      cat("with the ", j, " training ", training[j], "\n", file=logFile, append=TRUE)
      
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
  
  trainingNames <- training
  testNames <- test
  if(useFisher){
    trainingNames <- my.strsplit(trainingNames, "mapped", 2)
    testNames <- my.strsplit(testNames, "mapped", 2)
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
      testImg <- scan(paste(testDir, test[i], sep=""), quiet=TRUE)
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
      
      #cat("k=", k, " training[k]=", trainingNames[k], "\n")
      
      cl[j] <- getFaceID(trainingNames[k])
      
      if((method == "cosine") || (method == "dice"))
        dists[k] <- -100
      else
        dists[k] <- 100
    }
    
    closest[[i]] <- cl
    
    if(toFile != "")
      write(cl, concatenate(c(toFile, "_", getFaceID(testNames[i]), ".txt")), 1)
    
    cat(i*100/M, "%\n") #prints the progress
    cat("Testing ", testNames[i], ": ", crono.end(start), ", ", i*100/M, "%\n", file=logFile, append=TRUE)
  }
  (closest) #returns the 'amount' closest images
}

ponderateVote <- function(votes, by="min", type="number"){
  
  #gets the number of votes
  n <- length(votes[,1])
  
  #gets the unique candidates
  uVotes <- unique(votes[,1])
  #gets the number of unique candidates
  m <- length(uVotes)
  
  #gets the number of votes for each class
  nVotes <- tabulate(match(votes[,1], uVotes))
  #gets the bigger number of votes for all classes
  maxK <- max(nVotes)
  
  maxV <- max(votes[,2]) + 1
  
  results <- rep(0,m)
  
  for(i in 1:m){
    
    v <- which(votes[,1] == uVotes[i])
    k <- length(v)
    
    for(j in 1:k){
      
      if(votes[v[j], 2] == 0)
        votes[v[j],2] <- 0.0000001
      
      if(by == "min" && type == "number")
        results[i] <- results[i] + maxV - votes[v[j],2]
      else if(by == "max" || type == "value")
        results[i] <- results[i] + votes[v[j],2]
    }
    #results[i] <- results[i]/k
    
    if(type == "value")
      results[i] <- results[i]/k - k / (2 * maxK) * (maxV-1)
      #results[i] <- results[i]/k - (k^1.8)/results[i]
  }
  
  if(type == "number" || by == "max")
    return(c(uVotes[which.max(results)], max(results)))
  else if(type == "value" && by == "min")
    return(c(uVotes[which.min(results)], min(results)))
}

hierarchicalFeatureBasedPrediction <- function(model, testDir="", subset=integer(0), useErrorRange=TRUE, logFile=""){
  
  #gets the files' names
  testing <- dir(testDir)
  
  if(length(subset) > 0)
    testing <- testing[subset]
  
  #retrieves the classes information
  classes <- getClassFromFiles(files=testing)
  
  #gets the descriptors' values for all test samples
  descriptors <- lapply(concatenate(list(testDir, testing)), readMainLines, "list")
  
  #the number of models, one for each descriptor
  N <- length(model)
  #the number of testing sample
  M <- length(testing)
  
  tests <- list()
  
  #for each descriptor, ...
  for(i in 1:N){
    
    #separates only the vectors for the ith descriptor
    samples <- getAllFieldFromList(descriptors, i, 2)
    #puts them into a matrix
    samples <- list2matrix(samples)
    #puts this matrix into tests list
    tests[[i]] <- samples
  }
  
  corrects <- 0
  
  cat("Testing", M, "samples!\n", file=logFile, append=FALSE)
  
  #for each test sample, ...
  for(m in 1:M){
    
    cat("\npredicting test", m, ":\n", file=logFile, append=TRUE)
    
    start <- getTime()
    
    comparisons <- 0
    
    #initializes the votes as an empty list
    votes <- list()
    
    #for each descriptor, ...
    for(i in 1:N){
      
      cat("With predictor", i, "\n", file=logFile, append=TRUE)
      
      #initializes this descriptor's votes as a matrix with zeros (these zeros will be ignored later)
      votes[[i]] <- matrix(c(0,0), nrow=1)
      
      #gets the ith descriptor of the mth test sample
      test <- tests[[i]][m,]
      
      #gets the second level's representants
      secondLevel <- list2matrix(getAllFieldFromList(model[[i]], "representant", 2))
      #computes the errors for the second level
      icpResults <- apply(secondLevel, 1, function(reference, target){
        
        return (my.icp.2d.v2(reference, curveCorrection3(target, reference, 1), pSample=0.2, minIter=4))
        
      }, test)
      
      errors <- list2vector(getAllFieldFromList(icpResults, "error", 2))
      maxErrors <- list2vector(getAllFieldFromList(model[[i]], "maxError", 2))
      
      rangeCheck <- rep(TRUE, length(errors))
      
      if(useErrorRange){
        
        dists <- getAllFieldFromList(icpResults, "dist", 2)
        errorRange <- getAllFieldFromList(model[[i]], "errorRange", 2)
        rangeCheck <- checkErrorRange(errorRange, dists)
      }
      
      #retrieves which nodes of the second level matched the test
      passed2 <- which(errors <= maxErrors & rangeCheck)
      
      comparisons <- comparisons + dim(secondLevel)[1]
      
      #for each node from second level which matched the test, ...
      for(j in passed2){
        
        cat("2ndL(", j, "){", file=logFile, append=TRUE)

        #gets the first level's representants
        firstLevel <- list2matrix(getAllFieldFromList(model[[i]][[j]]$children, "representant", 2))
        
        passed1 <- 1
        
        #if there is more than one representant at this level
        if(!is.null(dim(firstLevel))){
          #computes the errors for the first level
          icpResults <- apply(firstLevel, 1, function(reference, target){
            
            return (my.icp.2d.v2(reference, curveCorrection3(target, reference, 1), pSample=0.2, minIter=4))
            
          }, test)
          
          errors <- list2vector(getAllFieldFromList(icpResults, "error", 2))
          maxErrors <- list2vector(getAllFieldFromList(model[[i]][[j]]$children, "maxError", 2))
          
          rangeCheck <- rep(TRUE, length(errors))
          
          if(useErrorRange){
            
            dists <- getAllFieldFromList(icpResults, "dist", 2)
            errorRange <- getAllFieldFromList(model[[i]][[j]]$children, "errorRange", 2)
            rangeCheck <- checkErrorRange(errorRange, dists)
          }
          #retrieves which nodes of the first level matched the test
          passed1 <- which(errors <= maxErrors & rangeCheck)
        }
        else{
          #computes the error with the single first level's representant
          icpResults <- my.icp.2d.v2(firstLevel, curveCorrection3(test, firstLevel, 1), pSample=0.2, minIter=4)
          maxError <- list2vector(getAllFieldFromList(model[[i]][[j]]$children, "maxError", 2))
          
          rangeCheck <- TRUE
          
          if(useErrorRange){
          
            errorRange <- getAllFieldFromList(model[[i]][[j]]$children, "errorRange", 2)
            rangeCheck <- checkErrorRange(errorRange, list(icpResults$dist))
          }
          
          #checks whether the representant matched
          passed1 <- which(icpResults$error <= maxError & rangeCheck)
        }
        
        comparisons <- comparisons + dim(firstLevel)[1]
        
        #for each node from first level which matched the test, ...
        for(k in passed1){
          
          cat(" 1stL(", k, "){", file=logFile, append=TRUE)
          
          #gets the nodes of the leaf level
          nodes <- model[[i]][[j]]$children[[k]]$children
          
          #gets the leafs' representants
          leafs <- list2matrix(getAllFieldFromList(nodes, "representant", 2))
          
          passed <- 1
          
          #if there is more than one representant at this level, ...
          if(!is.null(dim(leafs))){
            #computes the errors for the leafs
            icpResults <- apply(leafs, 1, function(reference, target){
              
              return (my.icp.2d.v2(reference, curveCorrection3(target, reference, 1), pSample=0.2, minIter=4))
              
            }, test)
            
            errors <- list2vector(getAllFieldFromList(icpResults, "error", 2))
            maxErrors <- list2vector(getAllFieldFromList(nodes, "maxError", 2))
            
            rangeCheck = rep(TRUE, length(errors))
            
            if(useErrorRange){
              
              dists <- getAllFieldFromList(icpResults, "dist", 2)
              errorRange <- getAllFieldFromList(nodes, "errorRange", 2)
              rangeCheck <- checkErrorRange(errorRange, dists)
            }
            
            #retrieves which nodes of the leafs level matched the test
            passed <- which(errors <= maxErrors & rangeCheck)
          }
          else{
            #computes the error with the single leafs' representant
            icpResults <- my.icp.2d.v2(leafs, curveCorrection3(test, firstLevel, 1), pSample=0.2, minIter=4)
            maxError <- list2vector(getAllFieldFromList(nodes, "maxError", 2))
            
            rangeCheck <- TRUE
            
            if(useErrorRange){
              
              errorRange <- getAllFieldFromList(nodes, "errorRange", 2)
              rangeCheck <- checkErrorRange(errorRange, list(icpResults$dist))
            }
            
            #checks whether the representant matched
            passed <- which(icpResults$error <= maxError & rangeCheck)
          }
          
          comparisons <- comparisons + dim(leafs)[1]
          
          #cat(names(nodes), "\n")
          
          #for each leaf that matched the test, ...
          for(v in passed){
            
            cat(" leaf(", names(nodes)[v], ")", file=logFile, append=TRUE)
            
            #gets the leaf's samples
            samples <- nodes[[v]]$samples
            #computes the errors for each leaf sample
            icpResults <- apply(samples, 1, function(reference, target){
              
              return (my.icp.2d.v2(curveCorrection3(reference, leafs[v,], 1), curveCorrection3(target, leafs[v,], 1), pSample=0.2, minIter=4))
              
            }, test)
            #gets the minimum computed error
            minError <- min(list2vector(getAllFieldFromList(icpResults, "error", 2)))
            
            #adds a vote for this leaf's class with the weight as the minimum error value
            #cat("leaf:", v, " descriptor:", i, "test:", m, "first level:", k, "second level:", j, "\n")
            votes[[i]] <- rbind(votes[[i]], matrix(c(as.numeric(names(nodes)[v]), minError), nrow=1))
          }
          cat("}", file=logFile, append=TRUE)
        }
        cat("}\n", file=logFile, append=TRUE)
      }
      #removes the initialization value
      votes[[i]] <- votes[[i]][-1,]
    }
    
    #counts the votes
    votes <- Reduce(rbind, votes, matrix(c(0,0), nrow=1))[-1,]
    
    cat("votes:\n", file=logFile, append=TRUE)
    cat.matrix(votes, file=logFile, append=TRUE)
    cat("comparisons:", comparisons, "\n", file=logFile, append=TRUE)
    cat("test", m, ". ", file=logFile, append=TRUE)
    
    if(is.null(dim(votes)))
      dim(votes) <- c(1,2)
    
    if(length(votes) > 1){
      
      result <- ponderateVote(votes, by="min", type="value")
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

hierarchicalFeatureBasedClassifier <- function(trainingDir, training=c()){
  
  if(length(training) == 0)
    training <- dir(trainingDir)
  
  #retrieves the classes information
  classes <- getClassFromFiles(files=training)
  
  descriptors <- lapply(concatenate(list(trainingDir, training)), readMainLines, "list")
  
  N <- length(descriptors[[1]])
  C <- length(classes$classes)
  
  model <- list()
  
  for(i in 1:N){
    
    cat("Computing tree for the", i, "th descriptor-------\n")
    
    #separates only the vectors for the ith descriptors
    samples <- getAllFieldFromList(descriptors, i, 2)
    #puts them into a matrix
    samples <- list2matrix(samples)
  
    #creates the first C nodes, where C = number of classes
    leafs <- computeNodes(samples, classes$fileClasses, progress=TRUE)
    
    #divides the leafs into groups
    groups <- computeGrouping(leafs, "brute", 7, progress=TRUE)
    #mounts the first level
    firstLevel <- computeNodes(list2matrix(getAllFieldFromList(leafs, "representant", 2)), groups, leafs, TRUE)
    
    #divides the first level into groups
    groups <- computeGrouping(firstLevel, "brute", 3, progress=TRUE)
    #mounts the second level
    secondLevel <- computeNodes(list2matrix(getAllFieldFromList(firstLevel, "representant", 2)), groups, firstLevel, TRUE)
    
    model[[i]] <- secondLevel
  }
  
  return(model)
}

computeNodes <- function(samples, groups, children=0, progress=FALSE){
  
  g <- unique(groups)
  G <- length(g)
  
  levels <- list()
  
  for(j in 1:G){
    
    node <- list()
    
    thisClassSamplesIndex <- which(groups == g[j])
    thisClassSamples <- samples[thisClassSamplesIndex,]
    
    node[["samples"]] <- thisClassSamples
    
    if(is.null(dim(thisClassSamples))){
      
      node[["representant"]] <- thisClassSamples
      node[["meanError"]] <- children[[thisClassSamplesIndex]]$meanError
      node[["maxError"]] <- children[[thisClassSamplesIndex]]$maxError
      node[["deviation"]] <- children[[thisClassSamplesIndex]]$deviation
      node[["errorRange"]] <- children[[thisClassSamplesIndex]]$errorRange
      node[["children"]] <- children[thisClassSamplesIndex]
    }
    else{
    
      meanClassSample <- colMeans(thisClassSamples)
      
      node[["representant"]] <- meanClassSample
      
      #compute the mean error and the mean error kind
      icpResults <- apply(thisClassSamples, 1, function(target, reference){
        
        return (my.icp.2d.v2(reference, curveCorrection3(target, reference, 1), pSample=0.2, minIter=4))
        
      }, meanClassSample)
      
      errors <- list2vector(getAllFieldFromList(icpResults, "error", 2))
      dists <- getAllFieldFromList(icpResults, "dist", 2)
      
      node[["meanError"]] <- mean(errors)
      node[["deviation"]] <- sd(errors)
      node[["maxError"]] <- max(errors)
      node[["maxError"]] <- node[["maxError"]] * 1.2
      node[["errorRange"]] <- computeErrorRanges(dists, 1.5)
      
      if(is.list(children)){
        node[["children"]] <- children[thisClassSamplesIndex]
        
        #computes the mean max error of the children
        maxErrors <- list2vector(getAllFieldFromList(children, "maxError", 2))
        maxErrors <- mean(maxErrors)
        if(node[["maxError"]] < maxErrors)
          node[["maxError"]] <- maxErrors
        
        errorRange <- getAllFieldFromList(children, "errorRange", 2)
        dists <- list()
        dists[[1]] <- matrix(c(node[["errorRange"]][,1], node[["errorRange"]][,2]), ncol=2)
        dists[[2]] <- matrix(c(node[["errorRange"]][,1], node[["errorRange"]][,3]), ncol=2)
        
        for(i in 1:(length(errorRange))){
          dists[[length(dists) + 1]] <- matrix(c(errorRange[[i]][,1], errorRange[[i]][,2]), ncol=2)
          dists[[length(dists) + 1]] <- matrix(c(errorRange[[i]][,1], errorRange[[i]][,3]), ncol=2)
        }
        
        node[["errorRange"]] <- computeErrorRanges(dists)
      }
    }
    
    levels[[g[j]]] <- node
    
    if(progress)
      cat(j*100/G, "%\n")
  }
  
  return(levels)
}

computeGrouping <- function(nodes, mode, nGroups=0, threshold=0, progress=FALSE){
  
  C <- length(nodes)
  
  #determine the maximum number of groups for the level 1
  if(nGroups <= 0)
    nGroups <- floor(C/6)
  
  #determines the thresholding on the similarity index for grouping
  if(threshold <= 0)
    threshold <- C/nGroups
  
  if(mode == "brute"){
    
    return( computeGroupingByBrute(nodes, nGroups, threshold, progress))
  }
  else if(mode == "Kmeans"){
    
    return( computeGroupingByKmeans(nodes, nGroups, threshold))
  }
}

computeGroupingByBrute <- function(nodes, nGroups=0, threshold=0, progress=FALSE){
  
  C <- length(nodes)
  
  #initiates the similarity matrix
  similarityMatrix <- matrix(rep(0, C*C), nrow=C)
  n <- computeNumberOfCombinations(C, 2)
  if(progress)
    cat(n, " combinations\n")
  m <- 1
  #computes the similarity matrix with the C representants
  for(j in 1:(C-1)){
    for(k in (j+1):C){
      
      similarityMatrix[j,k] <- my.icp.2d.v2(nodes[[j]][["representant"]],
                                            nodes[[k]][["representant"]], minIter=4, pSample=0.2)$error
      similarityMatrix[k,j] <- similarityMatrix[j,k]
      
      if(progress){
        cat("computing similarity:", m*100/n, "%\n")
        m <- m + 1
      }
    }
  }
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
  
  return(groups)
}

#' Not finished!
computeGroupingByKmeans <- function(nodes, nGroups, threshold){
  
  C <- length(nodes)
  
  #chooses the first K candidates as group ****
  Ks <- sample(1:C, nGroups)
  
  #initiates the similarity matrix
  similarityMatrix <- matrix(rep(0, nGroups*C), nrow=C)
  
  #computes the similarity matrix with the C representants
  for(j in 1:nGroups){
    
    representants <- c(1:C)[-Ks[j]]
    #sets the similiarity measure as 0 for the identity case
    similarityMatrix[j,Ks[j]] <- 0
    
    for(n in representants){
      
      similarityMatrix[j,n] <- my.icp.2d.v2(nodes[[Ks[j]]][["representant"]],
                                            nodes[[n]][["representant"]], minIter=2, pSample=0.33)
    }
  }
  #computes the rank of similarity for the whole matrix
  similarityMatrix <- t(apply(similarityMatrix, 1, rank, ties.method = "random"))
  similarityMatrix <- similarityMatrix - 1
  
  groups <- rep(0, C)
  
  #puts the first representant into the first group
  groups[Ks[1]] <- 1
  groupIndex <- 2
  
  #exchanges all zeros by the greatest possible value
  similarityMatrix[which(similarityMatrix == 0)] <- C + 1
  
  for(j in 1:nGroups){
    
    #gets the closest representant
    closest <- which.min(similarityMatrix[j,])
    #checks whether this closest is one of the groups center
    sameGroup <- Position(function(x){ return(x == closest)}, Ks)
    #if it is, ...
    if(!is.na(sameGroup)){
      
      #picks the possible candidates for new group center (any which hasn't been picked yet)
      possibleCandidates <- c(1:C)[-Ks]
      #randomly chooses the replacement
      k <- sample(possibleCandidates, 1)
      #makes the replacement
      Ks[sameGroup] <- k
      
      representants <- c(1:C)[-Ks[sameGroup]]
      #sets the similiarity measure as 0 for the identity case
      similarityMatrix[sameGroup,Ks[sameGroup]] <- 0
      
      for(n in representants){
        
        similarityMatrix[sameGroup,n] <- my.icp.2d.v2(nodes[[Ks[sameGroup]]][["representant"]],
                                              nodes[[n]][["representant"]], minIter=2, pSample=0.33)
      }
    }
    
    #gets the similarity index for the closest
    value <- similarityMatrix[j,closest]
    
    if(groups[j] == 0){
      
      #if this value is smaller or equal than the threshold, ...
      if(value <= threshold){
        #and if closest already has a group, ...
        if(groups[closest] != 0){
          
          #assigns the group of the closest to it
          groups[j] <- groups[closest]
        }
        else{
          
          #creates a new group and assigns it to this representant
          groups[j] <- groupIndex
          #and to its closest
          groups[closest] <- groupIndex
          #updates group index
          groupIndex <- groupIndex + 1
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
      }
    }
  }
  
  return(groups)
}

print.hierarchicalModel <- function(model, file=""){
  
  N <- length(model)
  
  for(i in 1:N){
    
    cat("(", i, file=file, append=TRUE)
    
    print.hierarchicalModelAux(model[[i]], 1, file)
    
    cat("\n)\n", file=file, append=TRUE)
  }
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

computeErrorRanges <- function(dists, add=1){
  
  N <- length(dists)
  
  errorRanges <- matrix(c(0,0,0), ncol=3)
  
  errorRanges <- Reduce(function(errorRanges, d){
    
    n <- length(d[,1])
    
    for(i in 1:n){
      
      x <- Position(function(z){
        return(z == d[i,1])
      }, errorRanges[,1])
      
      if(!is.na(x)){
        
        if(d[i,2] > errorRanges[x,3])
          errorRanges[x,3] <- d[i,2]
        
        if(d[i,2] < errorRanges[x,2])
          errorRanges[x,2] <- d[i,2]
      }
      else{
        
        errorRanges <- rbind(errorRanges, matrix(c(d[i,1], d[i,2], d[i,2]), ncol=3))
      }
    }
    
    return(errorRanges)
    
  }, dists, errorRanges)
  
  errorRanges <- errorRanges[-1,]
  
  errorRanges[,2] <- errorRanges[,2] - add
  errorRanges[,3] <- errorRanges[,3] + add
  
  return(errorRanges)
}

checkErrorRange <- function(ranges, r){
  
  return(mapply(function(reference, target){
    
    n <- length(target[,1])
    
    for(i in 1:n){
      
      x <- target[i,1]
      refX <- Position(function(z){ return(z == x)}, reference[,1], nomatch = 0)
      
      if(refX != 0 && (target[i,2] > reference[refX,3] || target[i,2] < reference[refX,2]))
        return(FALSE)
    }
    return(TRUE)
    
  }, ranges, r))
}