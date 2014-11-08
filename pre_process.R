#loads the needed packages
library(jpeg) #open jpg images
library(stinepack) #interpolate curves
library(abind) #combine multi-dimensional arrays

# Writes the weights into the specified file. ----
# input:
#   weights = a list containing for each class, a weight for each descriptor
#   file = a file path
# output:
#   none, the file will be created or overwritten with the given weights
write.weights <- function(weights, file){
  
  #gets the number of classes' weights
  n <- length(weights)
  #gets the classes' names
  classes <- names(weights)
  
  #for each class...
  for(i in 1:n){
    #gets the number of descriptors
    col <- length(weights[[i]])
    #writes the name of the class into the file
    write(classes[i], file, 1, append=TRUE)
    #writes the 'col' weights into the file just after the class name
    write(weights[[i]], file, col, TRUE)
  }
  cat("Done!\n")
}

# Reads the weights from the given file ----
# input:
#   file = a file path where is the weights
# ouput:
#   a list containing for each class, each descriptor's weights
read.weights <- function(file){
  
  lines <- readLines(file)
  n <- length(lines)/2
  weights <- list()
  
  for(i in 1:n){
    
    weights[[lines[i*2-1]]] <- as.numeric(strsplit(lines[i*2], " ")[[1]])
  }
  (weights)
}

write.3dArray <- function(data, file){
  
  d <- dim(data)
  for(i in 1:d[1]){
    for(j in 1:d[2]){
      
      write(data[i,j,], file, d[3], TRUE)
    }
  }
}

read.3dArray <- function(file, d2, d3){
  
  lines <- readLines(file)
  n <- length(lines)
  
  data <- array(rep(0, (n/d2)*d2*d3), c((n/d2), d2, d3))
  
  for(i in 1:(n/d2)){
    for(j in 1:d2){
      
      data[i,j,] <- as.numeric(strsplit(lines[(i-1)*d2 + j], "[ ]")[[1]])
    }
  }
  (data)
}

# Performs the outlier correction using a mean image and a threshold.
# the difference between the image and the mean image greater than the
# threshold are considered an outlier and eliminated.
# input:
#   img = either a file path (string) or a 2D matrix containing the range image
#   meanImg= either a file path (string) or a 2D matrix containing the mean range image
#   threshold = a number. By default, the threshold is 3 times the mean of the 
#               differnece between the image and the mean image
#   file = the file path where the resulting image shall be saved. By default, nothing will be saved.
#   progress = whether the progress shall be displayed. By default the progress won't be displayed.
# output:
#   a 2D matrix containing the range iamge without the detected outliers.
outlierCorrection4 <- function(img, meanImg, threshold=0, file="", progress=FALSE){
  
  #if 'img' is a filepath, reads it, otherwise, do nothing
  if(is.character(img))
    img <- readImageData(img)
  
  #if 'meanImg' is a filepath, reads it, otherwise, do nothing
  if(is.character(meanImg))
    meanImg <- readImageData(meanImg)
  
  #retrieves the numbers of column and row of the image
  col <- length(img[1,])
  row <- length(img[,1])
  
  #computes the absolute difference between the image and the mean image
  diffImg <- abs(img - meanImg)
  
  #if the threshold equals zero, then...
  if(threshold == 0)
    #defines it as 3 times the mean of the difference between the image and the mean image
    threshold = 2*mean(diffImg)
  
  #initializes the unoutlied image as the image itself
  unoutlied <- img
  #finds the points considered outlier, the ones whose difference 
  #between the image and the mean image are greater than the threshold
  occlusionMask <- which(diffImg > threshold)
  #all such points becomes zero (black)
  unoutlied[occlusionMask] <- 0
  #diffImg <- matrix(rep(0, col*row), ncol=col)
  #diffImg[occlusionMask] <- 1
  
  #if the result must be saved
  if(file != ""){
    
    #saves the JPEG difference image
    my.writeJPEG(diffImg, concatenate(c(file, ".diff.jpg")), quality=1)
    #saves the JPEG image without outliers
    my.writeJPEG(unoutlied, concatenate(c(file, ".jpg")), quality=1)
    #saves the data of the difference image
    write(t(diffImg), concatenate(c(file, "diff.dat")), ncolumns=col)
    #saves the data of the image without outliers
    write(t(unoutlied), concatenate(c(file, ".dat")), ncolumns=col)
  }
  
  #returns the image without outliers
  (unoutlied)
}

# Generates a mean face image with the training set ----
# input:
#   trainingDir = a string containing the path to the training directory
#   file = a string with the file path where the resulting image shall be saved.
#          By default, the result won't be saved.
#   progress = whether the progress should be displayed. By default, the
#              progress will be displayed.
# output:
#   a 2D matrix containing the mean image.
meanFace <- function(trainingDir, file="", isDir=TRUE, progress=TRUE){
  
  #retrieves all files withing the training directory
  training <- trainingDir
  if(isDir){
    training <- dir(trainingDir)
    #takes out all JPG files
    training <- training[-which(my.strsplit(training, "[.]", 0) == "jpg")]
  }
  else
    trainingDir <- ""
  
  #'imgs' will contain all images from the tranining directory
  imgs <- list()  
  #gets the number of training images
  m <- length(training)
  
  #for each of the training image ...
  for(i in 1:m){
    
    #reads the image from the file
    imgs[[i]] <- readImageData(concatenate(c(trainingDir, training[i])))
    
    #if the progress is to be displayed, show the percentage of images read
    if(progress)
      cat("reading imgs ", i * 100/m, "%\n")
  }
  
  #gets the number of columns and rows of the images
  col <- length(imgs[[1]][1,])
  row <- length(imgs[[1]][,1])
  
  #initializes the mean image with all pixels black (value 0)
  meanImg <- matrix(rep(0, col*row), ncol=col)
  
  #for each pixel of the mean image, ...
  for(i in 1:row){
    for(j in 1:col){
      
      #computes the mean value for the coordinates 'i', 'j'
      #initializes the mean value with zero
      meanDepth <- 0
      for(k in 1:m){
        #sums all values of that coordinates
        meanDepth <- meanDepth + imgs[[k]][i,j]
      }
      #divides the sum by the number of images, obtaining the mean
      meanDepth <- meanDepth/k
      
      #attributes the mean to the corresponding point of the mean image
      meanImg[i,j] <- meanDepth
      
      #if the progress is to be displayed, show the percentage of pixels filled.
      if(progress)
        cat("computing mean ", (i*col + j) * 100/(col*row), "%\n")
    }
  }
  
  #if the result is to be saved ...
  if(file != ""){
    
    #saves the JPEG mean iamge
    my.writeJPEG(meanImg, file, quality=1)
    #saves the data mean image
    write(t(meanImg), concatenate(c(file, ".dat")), ncolumns=col)
  }
  
  #returns the mean image matrix
  (meanImg)
}

# Uses the transformation matrix given by PhD. Clement Creusot ----
# by ICP with the first face to align all faces.
# input:
#   file = a string containing the file path to the ABS file
#   transf_data = a list containing all transformation data indexed by the face ID,
#                 obtained by the 'getTransfData' function.
#   ldmk = a string containing the filepath to the landmarks file. By default
#          no landmark information will be used and transformed.
#   progress = whether the progress should be displayed. By default the progress
#              will be displayed.
# output:
#   a list containing 2 elements:
#     'points' = a 2D matrix continaing the points transformed (3 columns, one for each dimensions)
#     'ldmk' = the landmark coordinates transformed or an empty string
registerFaces <- function(file, transf_data, ldmk="", progress=TRUE){
  
  #gets the face ID from the file name
  name <- getFaceID(file)
  
  #reads the ABS file into a 2D matrix with 3 columns, one for each dimension
  face <- abs2matrix(file)
  if(progress)
    cat("abs data read\n")
  
  #cuts the face using the center coordinate and radius given with the 'transf_data'
  #face <- cutTheFace(face, nose=transf_data[[name]]$cropCenter, 0)
  if(is.null(ldmk[["nose"]]))
    ldmk[["nose"]] <- matrix(transf_data[[name]]$cropCenter, nrow=1, dimnames=list(c(), c("x", "y", "z")))
  
  face <- cutTheFace(face, nose=ldmk[["nose"]], 0, 100)
  
  if(progress)
    cat("face cut\n")
  
  #gets the number of 3D points
  m <- length(face[,1])
  #adds the W coordinate
  face <- cbind(face, matrix(rep(1,m), ncol=1))
  
  #applies the transformation by multiplying the points with the given transformation matrix
  face <- t(transf_data[[name]]$transfMatrix %*% t(face))
  if(progress)
    cat("face rotated\n")
  #removes the W coordinate
  face <- face[,1:3] / face[,4]
  if(progress)
    cat("coordinates back to normal\n")
  
  #if the landmarks were given ...
  if(length(ldmk) > 1){
    
    #adds the center point into the landmarks list
    #ldmk[["center"]] <- matrix(transf_data[[name]]$cropCenter, nrow=1, dimnames=list(c(), c("x", "y", "z")))
    
    for(i in 1:length(ldmk)){
      
      #adds the W coordinate
      l <- cbind(ldmk[[i]], matrix(1, ncol=1))
      #applies the transformation matrix as well
      l <- t(transf_data[[name]]$transfMatrix %*% t(l))
      #removes the W coordinate
      ldmk[[i]] <- l[,1:3]/l[,4]
    }
    
    if(progress)
      cat("landmarks rotated\n")
  }
  #returns the list with the transformed points of the image and of the landmakrs
  (list(points=face, ldmk=ldmk))
}

# Reads the data of the registration file that contains ----
# the crop and transformations information
# input:
#   file = a string containing the file path to the TXT file which contains all transformation informations
# output:
#   a list indexed by the face IDs, containing:
#     'cropCenter' = a vector containing the 3 coordinates of the face center point
#     'cropRadius' = a number containing the radius value to be used to cut the face
#     'transfMatrix' = a 4x4 matrix containing the coefficients to transform the faces aligning them
getTransfData <- function(file){
  
  #reads all lines from the file
  transf_data <- readLines(file)
  #removes the first 7 lines (header)
  transf_data <- transf_data[-(1:7)]
  #splits the lines by white spaces
  transf_data <- strsplit(transf_data, "[ ]")
  
  #gets the number of lines read, each corresponding to 1 face image
  n <- length(transf_data)
  
  data <- list()
  
  #for each face image information, ...
  for(i in 1:n){
    #reads the transformation data
    data[[transf_data[[i]][1]]] <- list(cropCenter = as.numeric(strsplit(transf_data[[i]][4], ",")[[1]]),
                                        cropRadius = as.numeric(transf_data[[i]][5]),
                                        transfMatrix = matrix(as.numeric(strsplit(transf_data[[i]][6], "[;,]")[[1]]), ncol=4, byrow=TRUE))
  }
  
  #returns the list
  (data)
}

# Computes the errors for each training sample and their closest ----
computeErrorsByCurves <- function(trainingDir, closest, nClosest, nDescriptors, toFile="", meanImgs="", ldmk="", method="error", maxIter=10, 
                                  minIter=5, pSample=0.33, range=0, progress=TRUE){
  #gets the training file names
  training <- dir(trainingDir)
  
  #if the range is given, gets only the samples within the range
  if(range[1] != 0)
    training <- training[range[1]:range[2]]
  
  N <- length(training) #gets the number of training samples
  M <- nDescriptors #gets the number of descriptors
  
  #initializes the 3D array of errors with zeros
  #(the 1st dimension is the training samples, the 2nd dimension is the closest, the 3rd dimension is the descriptors)
  errors <- array(rep(0, N*nClosest*M), c(N, nClosest, M))
  
  #for each training sample
  for(i in 1:N){
    #reads the ith training sample
    img1 <- readMainLines(paste(trainingDir, training[i], sep=""), "list")
    #gets the number of descriptors of the ith sample
    M <- length(img1)
    
    #gets all closest samples of the ith training sample
    closestImg <- readLines(concatenate(c(closest, "cl_", getFaceID(training[i]), ".txt")))
    #gets only the 'nClosest'th closest samples
    closestImg <- closestImg[2:(nClosest+1)]
    #turns them into the face sample name
    closestImg <- getFaceID(closestImg)
    #gets the number of closest samples
    n <- length(closestImg)
    #intializes a cronometer
    start <- getTime()
    
    #for each closest sample...
    for(j in 1:n){
      #reads the jth closest sample
      img2 <- readMainLines(concatenate(c(trainingDir, closestImg[j], ".lines")), "list")
      #reads the corresponding mean face sample
      meanImg <- getMainLines(readImageData(concatenate(c(meanImgs, getPersonID(closestImg[j]), "d000.mean.jpg.dat"))),
                              readLandmark(concatenate(c(ldmk, closestImg[j], ".ldmk"))))
      
      #computes the distance between the ith training sample and the jth closest sample
      if(method == "error"){
        
        for(k in 1:M){
	  
          errors[i,j,k] <- my.icp.2d.v2(curveCorrection3(img1[[k]], meanImg[[k]]), curveCorrection3(img2[[k]], meanImg[[k]]), maxIter, minIter, pSample)$error
        }
        
      }
      else if(method == "energy"){
        
        for(k in 1:M)
          errors[i,j,k] <- my.icp.2d.v2(img1[[k]], img2[[k]], maxIter, minIter, pSample)$energyMean
        
      }
    }
    
    if(progress)
      cat("computing errors: ", i*100/N, "%", crono.end(start), "\n")
  }
  
  if(toFile != "")
    write.3dArray(errors, toFile)

  (errors)
}

# Analyzes the descriptors and weights them for each person ----
# according to its discrimination power
familiarityWeights <- function(trainingDir, closest, nClosest, nDescriptors, errors=0, toFile="", meanImgs="", ldmk="", method="error", maxIter=10, range=0, progress=TRUE){
  
  #retrieves all files from 'trainingDir' directory
  training <- dir(trainingDir)    
  
  #if the range was specified, disconsider all others outside the given range
  if(range[1] != 0)
    training <- training[range[1]:range[2]]

  N <- length(training) #the number of training samples
  M <- nDescriptors #the number of descriptors
  
  #if the errors aren't given, computes them
  if(length(errors) == 1)
    errors <- computeErrorsByCurves(trainingDir, closest, nClosest, nDescriptors, toFile, meanImgs, ldmk, method, maxIter, range, progress)
  
  imgsClasses <- getPersonID(training) #the classes of each training sample
  classes <- unique(imgsClasses) #the classes
  Nclasses <- length(classes) #the number of classes
  classMeans <- list() #the list containing the means for each descriptor for each class
  classDev <- list() #the list containing the standard deviation for each descriptor for each class
  notClassMeans <- list() #the list containing the means for the matches classXnotClass for each class
  notClassDev <- list() #the list containing the standard deviation for each descriptor for each mismatch
  classWeights <- list()
  singleton <- list()
  
  #maxError <- max(errors)
  
  #initializes the lists with zeros
  for(l in 1:Nclasses){
    classMeans[[classes[l]]] <- rep(0, M)
    classDev[[classes[l]]] <- rep(0, M)
    notClassMeans[[classes[l]]] <- rep(0, M)
    notClassDev[[classes[l]]] <- rep(0, M)
    classWeights[[classes[l]]] <- rep(0, M)
    singleton[[classes[l]]] <- FALSE
  }
  
  for(i in 1:N){
    #retrieves the closest samples from each
    closestImg <- readLines(concatenate(c(closest, "cl_", getFaceID(training[i]), ".txt")))
    #removes the first one and keeps only the next 'nClosest'th 
    closestImg <- closestImg[2:(nClosest+1)]
    #gets the class for each of the closest
    closestImg <- getPersonID(closestImg)
    
    #retrieves the indexes of the closest of the same class as the ith sample
    matches <- which(imgsClasses[i] == closestImg)
    
    notClassErrors <- 0
    
    if(length(matches) != 0){
      
      #initializes the matrix that will contain the errors for each descriptor for the ith class match
      classErrors <- errors[i, matches,]
      
      if(is.null(dim(classErrors))){
        classErrors <- matrix(classErrors, nrow=1)
        classErrors <- rbind(classErrors, classErrors)
      }
      
      #computes the mean for each descriptor
      classMeans[[imgsClasses[i]]] <- colMeans(classErrors)
      
      #computes the standard deviation for each descriptor
      classDev[[imgsClasses[i]]] <- apply(classErrors, 2, sd)
      
      #initializes the matrix that will contain the errors for each descriptor for the ith class mismatch
      notClassErrors <- errors[i, -matches,]
    }
    else{
      singleton[[imgsClasses[i]]] <- TRUE
      
      #initializes the matrix that will contain the errors for each descriptor for the ith class mismatch
      notClassErrors <- errors[i, ,]
    }
    
    #computes the mean for each descriptor
    notClassMeans[[imgsClasses[i]]] <- colMeans(notClassErrors)
    
    #computes the standard deviation for each descriptor
    notClassDev[[imgsClasses[i]]] <- apply(notClassErrors, 2, sd)
    
    if(progress)
      cat("\ncomputing classes weights: ", i*100/N, "%\n")
  }
  
  for(cl in classes){
    
    cM <- classMeans[[cl]]
    cD <- classDev[[cl]]
    ncM <- notClassMeans[[cl]]
    ncD <- notClassDev[[cl]]
    
    #cat("class errors:\t", cW/nC, "\nnot-Class errors:\t", ncW/tnC, "\n")
    #if the number of matches classXclass isn't zero...
    if(!singleton[[cl]]){
      
      for(i in 1:M){ #for each descriptor
        sep <- computeSeparability(list(mean=cM[i], sd=cD[i]), list(mean=ncM[i], sd=ncD[i]), 1)
        #if the mean of errors of the matches is smaller than the mean of errors of the mismatches...
        if(cM[i] < ncM[i]){
          #computes the weight normally
          classWeights[[cl]][i] <- 3^(cM[i]/(7*sep + 0.5) - (sep)^2)
          
          if(cl == "04222")
            cat("line ", i, " normal\n")
        }
        #if it isn't...
        else{
          
          classWeights[[cl]][i] <- 3^(cM[i]/(5) + (sep)^2)
          if(cl == "04222")
            cat("line ", i, " mean error\n")
        }
        
        if(cl == "04222" && classWeights[[cl]][i] > 10)
          cat("too big!\n")
      }
    }
    #if it is zero...
    else{
      classWeights[[cl]] <- 4/(ncM+1/3) + 8
      cat("Singleton\t", cl, "\n")
    }
    
    #classWeights[[cl]] <- classWeights[[cl]] - min(classWeights[[cl]])
    #cat("class number:", numberOfClasses[[cl]], " not class number:", totalNumberOfClasses[[cl]], "\n")
  }
  
  (classWeights)
}

#' Computes a separability measure between two data sets.
#' The data sets are represented by two statiscal measures: mean and standard deviation.
#' ---------------------------------------------------------------------------------------
#' @aliases computeSeparability
#' @param d1 = A list containing the mean and the standard deviation, representing the 1st data set.
#' @param d2 = A list representing the 2nd data set.
#' @param reference = 0 if the reference data size should be picked automatically, 1 if the
#'                    shall be 'd1', 2 if the reference shall be 'd2'
#' @author Paulo Hack, Msc student at UNICAMP.
computeSeparability <- function(d1, d2, reference=0){
  
  #computes the dimensions ranges of each data set
  range1 <- list(min = d1$mean - d1$sd, max = d1$mean + d1$sd)
  range2 <- list(min = d2$mean - d2$sd, max = d2$mean + d2$sd)
  
  #computes the size of each range
  range1$size <- range1$max - range1$min
  range2$size <- range2$max - range2$min
  
  #cat("size1=", range1$size, " size2=", range2$size, "\n")
  
  ranges <- list(range1, range2)
  
  commonRange <- list(min = 0, max = 0, size = 0)
  commonRange$min = ranges[[which.max(c(range1$min, range2$min))]]$min
  commonRange$max = ranges[[which.min(c(range1$max, range2$max))]]$max
  commonRange$size = commonRange$max - commonRange$min
  
  separability <- 0
  
  #checks which range is smaller, that one will be used as size reference
  if(commonRange$size > 0){
    
    if(reference == 0){
      if(range1$size < range2$size)
        separability <- 1 - commonRange$size/range1$size
      else
        separability <- 1 - commonRange$size/range2$size
    }
    else if(reference == 1)
      separability <- 1 - commonRange$size/range1$size
    else
      separability <- 1 - commonRange$size/range2$size
  }
  else
    separability <- 1 - commonRange$size
  
  (separability)
}

graphizeErrorsOfClass <- function(errors, the_class, training, closest, nClosest=100, descriptor=0, maxError=0, progress=TRUE){
  
  dims <- dim(errors)
  error_classes <- getPersonID(training)
  samplesOfClass <- which(error_classes == the_class)
  N <- length(samplesOfClass)
  error_classes <- error_classes[samplesOfClass]
  training <- training[samplesOfClass]
  
  errors <- errors[samplesOfClass,,]
  
  matchError <- rep(0, dims[3])
  nMatches <- 0
  misMatchError <- rep(0, dims[3])
  nMisMatches <- 0
  
  if(maxError == 0)
    maxError <- max(errors)
  
  if(descriptor == 0)
    plot(c(0,dims[3]), c(0, maxError), col="white")
  else
    plot(c(0,N*nClosest), c(0, maxError), col="white")
  
  for(i in 1:N){
    
    #retrieves the closest samples from each
    closestImg <- readLines(concatenate(c(closest, "cl_", getFaceID(training[i]), ".txt")))
    #removes the first one and keeps only the next 'nClosest'th 
    closestImg <- closestImg[2:(nClosest+1)]
    closestIDs <- closestImg
    #gets the class for each of the closest
    closestImg <- getPersonID(closestImg)
    
    if(descriptor == 0){
      for(j in 1:nClosest){
        #cat(imgsClasses[i], " = ", closestImg[j], "\n")
        
        #if we have a match classXclass...
        if(error_classes[i] == closestImg[j]){
          for(k in 1:dims[3]){
            #sums the errors
            matchError[k] <- matchError[k] + errors[i,j,k]
            points(k-1, errors[i,j,k], col="blue", pch=7)
          }
          #counts the occurence
          nMatches <- nMatches + 1
        }
        else{ #if we have a match classXnotClass...
          for(k in 1:dims[3]){
            #sums the errors
            misMatchError[k] <- misMatchError[k] + errors[i,j,k]
            points(k-1, errors[i,j,k], col="red", pch=9)
          }
          #counts the occurence
          nMisMatches <- nMisMatches + 1
        }
      }
    }
    else{
      
      k <- descriptor
      
      for(j in 1:nClosest){
        #cat(imgsClasses[i], " = ", closestImg[j], "\n")
        
        #if we have a match classXclass...
        if(error_classes[i] == closestImg[j]){
          #sums the errors
          matchError[k] <- matchError[k] + errors[i,j,k]
          points((i-1)*nClosest + j, errors[i,j,k], col="blue", pch=7)
          #counts the occurence
          nMatches <- nMatches + 1
        }
        else{ #if we have a match classXnotClass...
          
          #sums the errors
          misMatchError[k] <- misMatchError[k] + errors[i,j,k]
          points((i-1)*nClosest + j, errors[i,j,k], col="red", pch=9)
          #counts the occurence
          nMisMatches <- nMisMatches + 1
        }
      }
    }
    
    cat(closestIDs[which.min(errors[i,,1])], "\n")
  
    if(progress)
      cat("\ncomputing classes weights: ", i*100/N, "%\n")
  }
  
  if(descriptor == 0)
    for(k in 1:dims[3]){
      
      points(k-1, matchError[k]/nMatches, col="green", pch=6)
      points(k-1, misMatchError[k]/nMisMatches, col="orange", pch=2)
    }
  else{
    
    points(-20, matchError[descriptor]/nMatches, col="green", pch=6)
    points(-20, misMatchError[descriptor]/nMisMatches, col="orange", pch=2)      
  }
}

# Retrieves the file names of all samples from the given class(es)-----
getAllSamplesFromClass <- function(folder, classes, toFile=""){
  
  files = dir(folder)
  exFiles <- list()
  k <- 1
  
  for(file in files){
    
    if(length(which(classes == getPersonID(file)))){
      exFiles[[k]] <- concatenate(c("rm ", folder, file))
      k <- k + 1
    }
  }
  
  if(toFile != "")
    write(list2vector(exFiles), toFile, 1)
  
  (exFiles)
}

my.icp.3d <- function(reference, target, matchingMethod="", sampling=TRUE, maxIter=10, threshold=0, progress=FALSE){
  
   N <- length(target[1,]) #gets the number of dimensions
   M <- length(target[,1]) #gets the number of points
   nSampling <- M/(round(M^(1/3))) #defines the number of samples used
   
   #samples some points
   tarSamples <- 0
   refSamples <- 0
   if(sampling){
     #samples <- sample(1:M, nSampling)
     tarSamples <- my.sample(target, nSampling)
     refSamples <- my.sample(reference, nSampling)
   }
   else{
     
     tarSamples <- 1:M
     refSamples <- 1:M
   }
   
   #performs the first transformation (initial guess) 
   #a neutral rotation
   Rotation <- matrix(rep(0, N*N), ncol=N)
   for(i in 1:N)
     Rotation[i,i] <- 1
   
   #a mean translation
   Translation <- matrix(rep(0,N), ncol=N)
#    for(i in 1:N)
#     Translation[1,i] <- ((max(reference[,i]) - max(target[,i])) + (min(reference[,i]) - min(target[,i])))/2
   
   #a neutral scaling
   #Scaling <- matrix(rep(1,N), ncol=N)
   
   #initializes the set of transformations
   #transformations <- list(list(R=Rotation, T=Translation, S=Scaling))
   transformations <- list()
   
   #computes the first error
   error <- threshold + 1
   
   #starts the iterations
   i <- 1
   while(error > threshold && i <= maxIter){
     
     transformations[[i]] <- list(R=Rotation, T=Translation)
     
     #performs the former transformations
     target <- (target %*% t(Rotation)) + matrix(rep(Translation[1,], M), nrow=M, byrow=TRUE)
     
     pairs <- list()
     j <- 1
     #computes the closest point for each sample
     for(p in tarSamples){
       
       start <- getTime()
       
       if(matchingMethod == ""){
         
         pairs[[as.character(j)]] <- findMatch(reference[refSamples,], target[p,], 0, M, N)
         j <- j + 1
       }
       else{
         
         distances <- as.matrix(dist(rbind(target[p,], reference[refSamples,])))[1,-1]
         pairs[[as.character(j)]] <- list(point=which.min(distances), distance=min(distances))
         j <- j + 1
       }
       
       cat(j, "th sample matched!\t", crono.end(start), "\n")
     }
     
     #removes the outliers
     distances <- getAllFieldFromList(pairs, "distance", 2)
     distances <- list2vector(distances)
     
     pairs <- pairs[-which(distances >= mean(distances) + sd(distances))]
     
     #builds the matrices with the collected points
     m <- length(pairs)
     tar <- matrix(rep(0, m*N), ncol=N)
     ref <- tar
     for(j in 1:m){
       
       tar[j,] <- target[as.numeric(names(pairs[j])),]
       ref[j,] <- reference[pairs[[names(pairs[j])]]$point,]
     }
     
     #finds the best rotation----
     muRef <- colMeans(ref) #gets the mean point of the reference
     muTar <- colMeans(tar) #gets the mean point of the target
     
     #computes the cross covariance matrix
     tar <- tar - matrix(rep(muTar, m), nrow=m, byrow=TRUE)
     ref <- ref - matrix(rep(muRef, m), nrow=m, byrow=TRUE)
     CrossCovariance <- (t(tar) %*% (ref))/M
     A <- CrossCovariance - t(CrossCovariance)
     Delta <- c(A[2,3], A[3,1], A[1,2])
     
     Q <- matrix(rep(0, (N + 1)^2), ncol=(N + 1))
     I <- matrix.identity(N)
     tr <- matrix.trace(CrossCovariance)
     Q[1,1] <- tr
     Q[1,-1] <- Delta
     Q[-1,1] <- Delta
     Q[-1,-1] <- CrossCovariance + t(CrossCovariance) - tr * I
     
     E <- eigen(Q, symmetric=TRUE)
     
     q <- E$vectors[,1]
     
     Rotation <- matrix(c((q[1]^2 + q[2]^2 - q[3]^2 - q[4]^2), (2*(q[2]*q[3] + q[1]*q[4])), (2*(q[2]*q[4] - q[1]*q[3])),
                   (2*(q[2]*q[3] - q[1]*q[4])), (q[1]^2 + q[3]^2 - q[2]^2 - q[4]^2), (2*(q[3]*q[4] + q[1]*q[2])),
                   (2*(q[2]*q[4] + q[1]*q[3])), (2*(q[3]*q[4] - q[1]*q[2])), (q[1]^2 + q[4]^2 - q[2]^2 - q[3]^2)), ncol=3)
     
     #finds the best scaling ----
#      sizesOfRef <- rep(0, N)
#      sizesOfTar <- rep(0, N)
#      for(j in 1:N){
#        
#        sizesOfRef[j] <- max(ref[,j]) - min(ref[,j])
#        sizeOfTar[j] <- max(tar[,j]) - min(tar[,j])
#        Scaling[1,j] <- sizesOfRef[j]/sizesOfTar[j]
#      }
     
     #finds the best translation ----
     Translation <- muRef - (muTar %*% Rotation)
     #Translation <- muRef - (Scaling * muTar) %*% Rotation
     Translation <- matrix(Translation, ncol=N)

     #computes the error
     tar <- (tar %*% t(Rotation)) + matrix(rep(Translation, m), ncol=N, byrow=TRUE)
     error <- 0
     for(j in 1:m){
       
       for(n in 1:N){
         
         error <- error + (ref[j,n] - tar[j,n])^2
       }
     }
     error <- error/m

     if(progress)
       cat("iteration ", i, "\n")

     i <- i + 1
   }

   transformations[[i]] <- list(R=Rotation, T=Translation)
    
   #performs the former transformations
   target <- (target %*% t(Rotation)) + matrix(rep(Translation[1,], M), nrow=M, byrow=TRUE)

   (list(transformations=transformations, target=target, error=error))
}

# Returns a identity matrix size x size----
matrix.identity <- function(size){
  
  I <- matrix(rep(0, size^2), ncol=size)
  for(i in 1:size)
    I[i,i] <- 1
  
  (I)
}

getCircledValue <- function(size, value, range=0){
  
  if(length(range) == 1)
    range <- c(1,size)
  
  return((value - range[1]) %% size + range[1])
}

matrix.trace <- function(m){
  
  n <- length(m[,1])
  tr <- 0
  for(i in 1:n)
    tr <- tr + m[i,i]
  
  (tr)
}

# Finds the target's closest point in the reference data ####
# input:
#   reference = a set of N-dimensional points;
#   target = a N-dimensional point;
#   guess = the initial guess range, a N-dimensional point.
#   M = the number of points in reference
#   N = the number of dimensions
# output:
#   a list with the index of the closest point in reference
#   and the distance
findMatch <- function(reference, target, guess=0, M=0, N=0){
  
  N <- length(reference[1,]) #gets the number of dimensions
  M <- length(reference[,1]) #gets the number of points
  
  #cat("N =", N, "M =", M, "\n")
  #print(reference)
  #print(target)
  
  #if the guess is not passed, chooses it automaticly
  if(guess == 0){
    
    guess <- rep(0, N)
    
    for(i in 1:N){
      guess[i] <- (max(reference[,i]) - min(reference[,i]))/(M^(1/N))
      #cat("guess ", i, ": ", guess[i], "\n")
    }
    
    guess <- mean(guess[which(guess != 0)])
    if(is.na(guess))
      guess <- 1
  }
  
  #initializes the candidates of reference to closest points with no point
  candidates <- list()
  cIndex <- 1
  
  while(length(candidates) == 0){
    
    for(i in 1:M){
      
      passed <- 0
      
      for(j in 1:N)
        if(reference[i,j] >= target[j] - guess && reference[i,j] <= target[j] + guess)
          passed <- passed + 1
      
      if(passed == N){
        
        candidates[[cIndex]] <- i
        cIndex <- cIndex + 1
        #cat("candidate ", cIndex-1, ": ", i, "\n")
      }
    }
    
    #cat("iteration ", cIndex, "\n")
    
    guess <- 2 * guess
  }
  candidates <- list2vector(candidates)
  
  #computes the distances
  distances <- as.matrix(dist(rbind(target, reference[candidates,])))[1,-1]

  #and returns the closest point with the correspoding distance  
  return(list(point=candidates[which.min(distances)], distance=min(distances)))
}

# Computes the curvature measure in each point of a given vector
# input:
#   data = a vector of numbers
#   factor = the proportions to consider each number of 'data'.
#            By default, 1.
# output:
#   a vector of numbers containing the curvature measures
curvatureVector <- function(data, factor=1){
  
  #gets the number of elements of 'data'
  n <- length(data)
  
  #initializes the result with 0s
  curvature <- rep(0, n)
  
  #for each element...
  for(i in 1:n){
    
    #computes the angle between the 'i'th point and its successor
    a <- 0
    #checks whether they are both different than zero and there is a successor
    if(i+1 <= n && data[i] != 0 && data[i+1] != 0)
      a <- abs(atan((factor*data[i] - factor*data[i+1])/(-1)))
    
    #computes the angle between the 'i'th point and its predecessor
    b <- 0
    #checks whether they are both differente than zero and there is a predecessor
    if(i-1 > 0 && data[i] != 0 && data[i-1] != 0)
      b <- abs(atan((factor*data[i] - factor*data[i-1])/(1)))
    
    #if both angles are different than zero...
    if(a != 0 && b != 0)
      #computes the curvature measures considering PI range
      curvature[i] <- pi/(pi - (a + b)) - 1
    else
      #otherwise, computes the curvature measure considereing PI/2 range
      curvature[i] <- (pi/2)/((pi/2) - (a + b)) - 1
  }
  
  #returns the result
  (curvature)
}

####### NOT USED FUNCTION #########
# Computes the rotation factors given the distance values
# and 2 lines/curves in order to approach the 2nd curve to 
# the 1st.
# input:
rotation.factors <- function(distances, reference, target){
  
  m <- length(distances)
  
  x1 <- min(reference[1,1], target[1,1])
  x2 <- min(reference[m,1], target[m,1])
  
  xm <- round((x2 - x1)/2)
  x1 = xm; x2 = xm;
  
  stopFlag <- FALSE
  
  charge <- 1
  if(distances[xm] < 0)
    charge <- -1
  
  while(x1 >= 1 && x2 <= m && !stopFlag){
    
    x1 <- x1-1
    if(charge == 1 && distances[x1] < 0){
      stopFlag <- TRUE
      xm <- x1+1
    }else if(charge == -1 && distances[x1] >= 0){
      stopFlag <- TRUE
      xm <- x1+1
    }
    
    x2 <- x2+1
    if(charge == 1 && distances[x2] < 0){
      stopFlag <- TRUE
      xm <- x2-1
    }else if(charge == -1 && distances[x2] >= 0){
      stopFlag <- TRUE
      xm <- x2-1
    }
  }
  
  if(stopFlag){
   
    xMax <- max(distances)
    x1 <- which(distances == xMax)[1]
    if(x1 >= xm)
      x1 <- atan(xMax/(x1 -xm))
    else
      x1 <- atan(xMax/(x1 -xm))
    
    xMin <- min(distances)
    x2 <- which(distances == xMin)[1]
    if(x2 >= xm)
      x2 <- atan(xMin/(x2 -xm))
    else
      x2 <- atan(xMin/(x2 -xm))
    
    if(abs(x1) > abs(x2))
      (list(p = xm, angle=x2))
    else
      (list(p=xm, angle=x1))
  }
  else
    (list(p=0, angle=0))
}

# Computes the scale factor to approach the target curve
# to the reference curve by performing a translation.
# input:
#   reference = a curve represented by a 2D matrix with 2 columns, 1 for each dimension
#   target = a curve as 'reference'
# output:
#   a number representing the translation factor
getXScaleFactor <- function(reference, target){
  
  #gets the number of points of the 'reference' curve
  m <- length(reference[,1])
  #finds the point with the greatest value in 'reference'
  h <- which.max(reference[,2])
  
  #finds the 
  refXvar1 <- getXvariation(reference[1:h,])
  tarXvar1 <- getXvariation(target[which(target[,1] <= h),])
  
  refXvar2 <- getXvariation(reference[h:m,])
  tarXvar2 <- getXvariation(target[which(target[,1] >= h),])
  
  scale1 <- 1
  if(tarXvar1$max - tarXvar1$min != 0)
    scale1 <- (refXvar1$max - refXvar1$min)/(tarXvar1$max - tarXvar1$min)
  
  scale2 <- 1
  if(tarXvar2$max - tarXvar2$min != 0)
    scale2 <- (refXvar2$max - refXvar2$min)/(tarXvar2$max - tarXvar2$min)
  
  if(abs(scale1 - 1) > abs(scale2 - 1))
    (scale1)
  else
    (scale2)
}

# Gets the max and min of the x coordinate (1st column) disconsidering
# the extreme regions, where the points forms a line parallel with the
# domain axis
getXvariation <- function(m){
  
  if(!is.null(dim(m))){
    n <- length(m[,1])
    xmin <- 0
    xmax <- n+1
    
    for(i in 1:(n-1)){
      
      if(m[i,2] == 0 && m[i,2] == m[i+1,2]){
        
        xmin <- i
      }
      else{
        break
      }
    }
    
    for(i in n:2){
      
      if(m[i,2] == 0 && m[i,2] == m[i-1,2]){
        
        xmax <- i
      }
      else{
        break
      }
    }
    
#     if(xmin == 0)
#       xmin <- m[xmin+1,1]
#     else
#       xmin <- m[xmin,1]
#     
#     if(xmax == n+1)
#       xmax <- m[xmax-1,1]
#     else
#       xmax <- m[xmax,1]-1
    
    xmax <- m[xmax-1,1]
    xmin <- m[xmin+1,1]
    
    (list(max = xmax, min = xmin))
  }
  else{
    
    (list(max=m[1], min=m[1]))
  }
  
}

# computes the factor of growing in each point by the target's coordinates
factor.p2p <- function(reference, target, distances){
  
  n <- length(target[,1])
  
  factors <- rep(0, n)
  
  #translate them so all 'y' will be positive
  minRefY <- min(reference[,2])
  minTarY <- min(target[,2])
  
  if(minTarY < minRefY && minTarY < 0){
    
    reference[,2] <- reference[,2] - minTarY
    target[,2] <- target[,2] - minTarY
  }
  
  j <- 1
  
  for(i in 1:n){
    
    x <- target[i,1]
    #j <- 1
    
    if(length(which(reference[,1] == x)) > 0){
      if(reference[which(reference[,1] == x),2] != 0){
        factors[i] <- distances[j] / max(target[i,2], reference[which(reference[,1] == x),2])
        #cat("OK\ti:", i, "j:", j, "f:", factors[i], "\n")
      }
      else{
        factors[i] <- distances[j]
        #cat("Bad\ti:", i, "j:", j, "f:", factors[i], "\n")
      }
      
      j <- j + 1
    }
  }
  
  (factors)
}

# Reads the main Lines from file returning either a vector or a list
# input:
#   file = the path to the file containing the lines
#   type = a string to tell whether a vector or a list must be returned
#          if it equals "vector", a vector is returned, otherwise, a list is returned
# output:
#   a vector with all lines in sequence or a list with each line separeted
readMainLines <- function(file, type="vector"){
  
  #reads the file, separing all each line and each value into a list
  resp <- strsplit(readLines(file), "[ ]")
  
  if(type == "vector"){ #if it is to return a vector
    
    #initiates a vector containing only the 1st line
    v <- resp[[1]]
    for(i in 2:length(resp)){
      
      v <- c(v, resp[[i]]) #adds all other lines into the vector
    }
    
    #parses the string values into numeric
    resp <- as.numeric(v)
  }
  else{ #if a it is to return a list
    
    for(i in 1:length(resp)){
      
      #parses the string values into numeric
      resp[[i]] <- as.numeric(resp[[i]])
    }
  }
  
  (resp)
}

# Gets 4 lines to describe the face: horizontal nose, perfil, horizontal eyes and forehead
# input:
#   img = a 2D matrix containing the image
#   ldmk = a list containing the landmarks nose, subnasale, nasion, righteye and lefteye
#   file = the path where the lines will be saved as a DAT file, containing the image values
#          of each line as strings in sequence at each line of the file.
#          if file = "", then no file will be created
#   halfFace = a boolean value indicating whether only half face will be used
#   smoothness = a numeric value indicating the intensity the image value must be smoothed
#                if smoothness = 0, then no change is done
# output:
#   a list containing the values of each line
getMainLines <- function(img, ldmk, file="", halfFace=FALSE, smoothness=0){
  
  col <- length(img[1,])
  row <- length(img[,1])
  
  lines <- list()
  
  #measures the distance between the nasion and the nose tip
  h <- dist(rbind(ldmk[["nasion"]][-3], ldmk[["nose"]][-3]))[1]
  
  #the first line passes through the nasion and the subnasale
  if(is.null(ldmk[["subnasale"]]))
    perfil <- findLine(ldmk[["nasion"]][-3], ldmk[["nose"]][-3], TRUE)
  else
    perfil <- findLine(ldmk[["nasion"]][-3], ldmk[["subnasale"]][-3], TRUE)
  
  lefteye <- list(a = perfil$a, b = perfil$b - (appLinear(perfil, ldmk[["lefteye"]][-3])[1] - ldmk[["lefteye"]][1] + (h/4)), inverted = perfil$inverted)
  righteye <- list(a = perfil$a, b = perfil$b + (-appLinear(perfil, ldmk[["righteye"]][-3])[1] + ldmk[["righteye"]][1] + (h/4)), inverted = perfil$inverted)
  
  #if the perfil line isn't perfectly vertical, we need the same line in function of x
  perfil2 <- perfil
  if(perfil2$a != 0){
    
    if(is.null(ldmk[["subnasale"]]))
      perfil2 <- findLine(ldmk[["nasion"]][-3], ldmk[["nose"]][-3])
    else
      perfil2 <- findLine(ldmk[["nasion"]][-3], ldmk[["subnasale"]][-3])
  }
  
  #the nose tip line is perpendicular to the perfil line, and it passes through the nose tip point
  nose <- findPerpendicularLine(perfil2, ldmk[["nose"]][-3])
  nose1 <- list(a = nose$a, b = nose$b - (h/4), inverted = nose$inverted)
  nose2 <- list(a = nose$a, b = nose$b + (h/4), inverted = nose$inverted)
  
  #the eyes line passes through the eyes
  eyes <- findLine(ldmk[["lefteye"]][-3], ldmk[["righteye"]][-3])
  eyes1 <- list(a = eyes$a, b = eyes$b - (h/5), inverted = eyes$inverted)
  eyes2 <- list(a = eyes$a, b = eyes$b + (h/5), inverted = eyes$inverted)
  #the forehead is perpendicular to the perfil line and is 1.8*h above the nose line
  forehead <- list(a = nose$a, b = nose$b - 1.8*h, inverted = nose$inverted)
  
  #the chin is perpendicular to the perfil line
  chin <- list(a = nose$a, b = nose$b + 1.2*h, inverted = nose$inverted)
  
  #writeImage(t(lines2image(list(perfil, nose, eyes, forehead), 150, 210, TRUE)), "lines.jpg", "jpg")
  
  #the perfil, right and left eye lines will end after the nose, so we stablished it as y=130
  lines[["perfil"]] <- rep(0, 130)
  lines[["righteye"]] <- rep(0,130)
  lines[["lefteye"]] <- rep(0,130)
  #for each y from 1 to 130, 
  for(i in 1:130){
    # finds the corresponding x of the perfil line, and the image value at that point
    x <- round(appLinear(perfil, i))
    if(x > 0 && x <= col)
      lines[["perfil"]][i] <- img[i, x]
    
    x <- round(appLinear(righteye, i))
    if(x > 0 && x <= col)
      lines[["righteye"]][i] <- img[i, x]
    
    x <- round(appLinear(lefteye, i))
    if(x > 0 && x <= col)
      lines[["lefteye"]][i] <- img[i, x]
  }
  
  #all other lines crosses all the image horizontally
  lines[["nose"]] <- rep(0, col)
  lines[["nose1"]] <- rep(0, col)
  lines[["nose2"]] <- rep(0, col)
  lines[["eyes"]] <- rep(0, col)
  lines[["eyes1"]] <- rep(0, col)
  lines[["eyes2"]] <- rep(0, col)
  lines[["forehead"]] <- rep(0, col)
  lines[["chin"]] <- rep(0, col)
  
  #for each x from 1 to col
  for(i in 1:col){
    
    #finds the corresponding y of the lines, and the image value at that point
    y <- round(appLinear(nose, i))
    if(y > 0 && y <= row)
      lines[["nose"]][i] <- img[y, i]
    
    y <- round(appLinear(nose1, i))
    if(y > 0 && y <= row)
      lines[["nose1"]][i] <- img[y, i]
    
    y <- round(appLinear(nose2, i))
    if(y > 0 && y <= row)
      lines[["nose2"]][i] <- img[y, i]
    
    y <- round(appLinear(eyes, i))
    if(y > 0 && y <= row)
      lines[["eyes"]][i] <- img[y, i]
    
    y <- round(appLinear(eyes1, i))
    if(y > 0 && y <= row)
      lines[["eyes1"]][i] <- img[y, i]
    
    y <- round(appLinear(eyes2, i))
    if(y > 0 && y <= row)
      lines[["eyes2"]][i] <- img[y, i]
    
    y <- round(appLinear(forehead, i))
    if(y > 0 && y <= row)
      lines[["forehead"]][i] <- img[y, i]
    
    y <- round(appLinear(chin, i))
    if(y > 0 && y <= row)
      lines[["chin"]][i] <- img[y, i]
  }
  
  #writeImage(lines2image(list(perfil, righteye, lefteye, forehead, chin, eyes, eyes1, eyes2, nose, nose1, nose2), 150, 210, FALSE), "04217dlinetest.jpg", "jpg")
  
  #if only half face is to be used
  if(halfFace){
    
    rightEyeEnergy <- measureEnergy(lines[["righteye"]])$meanEnergy
    leftEyeEnergy <- measureEnergy(lines[["lefteye"]])$meanEnergy
    
    #if the energy in the right eye is smaller
    if(rightEyeEnergy < leftEyeEnergy){
      #than, we take the right side
      for(i in 4:(length(lines)))
        lines[[i]] <- lines[[i]][-(1:round(col/2))]
    }
    else{
      #else, we take the left side
      for(i in 4:(length(lines)))
        lines[[i]] <- lines[[i]][-(round(col/2):col)]
    }
  }
  
  #if a file is to be created
  if(file != ""){
    
    #if a smooth operation is to be made
    if(smoothness > 0){
      #writes the smoothed lines' values into the file
      write(gaussianSmooth(lines[["perfil"]], c(smoothness)), file, 130)
      write(gaussianSmooth(lines[["righteye"]], c(smoothness)), file, 130, TRUE)
      write(gaussianSmooth(lines[["lefteye"]], c(smoothness)), file, 130, TRUE)
      
      for(i in 4:length(lines)){
        write(gaussianSmooth(lines[[i]], c(smoothness)), file, col, TRUE)
      }
    }
    else{
      #writes the lines' values into the file without changes
      write(lines[["perfil"]], file, 130)
      write(lines[["righteye"]], file, 130, TRUE)
      write(lines[["lefteye"]], file, 130, TRUE)
      
      for(i in 4:length(lines)){
        write(lines[[i]], file, col, TRUE)
      }
    }
  }
  #returns the list with the lines
  (lines)
}

# Applies the weight image into the range image
# input:
#   img = either a path to the image or a 2D matrix containing the image
#   weightImg = either a path to the weight image or a 2D matrix containing it
#   factor = a number to multiply the weights
#   file = the path where the new weighted image shall be saved
#          if file = "", then no file will be created
# output:
#   a 2D matrix containing the new weighted image
ponderateImage <- function(img, weightImg, factor=1, file=""){
  
  #if a path is given,
  if(is.character(img))
    #reads the image from the file
    img <- readImageData(img)
  
  #if a path is given,
  if(is.character(weightImg))
    #reads the image from the file
    weightImg <- readJPEG(weightImg)
  
  col <- length(img[1,])
  row <- length(img[,1])
  
  #initiates the weighted image
  weightedImg <- matrix(rep(0, col * row), ncol=col)
  
  #computes the weighted image, that is the image ponderated by the factor and the weights
  weightedImg <- img * factor * weightImg
  
  #if the weighted image is to be saved into a file
  if(file != ""){
    
    my.writeJPEG(weightedImg, file, quality=1)
    write(t(weightedImg), concatenate(c(file, ".dat")), ncolumns=col)
  }
  
  #returns the matrix containing the weighted image
  (weightedImg)
}

# Generates the curvature image.
# input:
#   rangeImage = either a path to the image file or a 2D matrix containing the image
#   weightImage = either a path to the image file or a 2D matrix containing the image
#   file = the path where the new image shall be saved, if it equals "", then the image
#          won't be saved
#   factor = the factor that will multiply the weights
#   printProgress = a boolean value to tell whether show the progress or not
# output:
#   a 2D matrix containing the curvature image
curvatureImage <- function(rangeImage, weightImage="", file="", factor=1, printProgress=FALSE){
  
  #if a path is given to the range image
  if(is.character(rangeImage))
    #reads it from the file
    rangeImage <- readImageData(rangeImage)
  
  #if a path is given to the weight image
  if(is.character(weightImage) && weightImage != "")
    #reads it from the file 
    weightImage <- readJPEG(weightImage)
  
  col <- length(rangeImage[1,])
  row <- length(rangeImage[,1])
  
  #initializes the curvature image
  curves <- matrix(rep(0, col * row), ncol = col)
  
  if(length(weightImage) > 1){
    #for each inner point of the image, computes the curvature
    for(y in 2:(row - 1)){
      for(x in 2:(col - 1)){
        
        #computes the horizontal angle of the point with its predecessor
        a <- 0
        if(rangeImage[y, x-1] != 0 && rangeImage[y, x] != 0)
          a <- abs(atan((factor*rangeImage[y,x] - factor*rangeImage[y, x-1])/(1)))
        
        #computes the horizontal angle of the point with its successor
        b <- 0
        if(rangeImage[y, x+1] != 0 && rangeImage[y, x] != 0)
          b <- abs(atan((factor*rangeImage[y,x] - factor*rangeImage[y, x+1])/(-1)))
        
        #computes the horizontal curvature
        curvesH <- 0
        if(a != 0 && b != 0)
          curvesH <- pi/(pi - (a + b)) - 1
        else
          curvesH <- (pi/2)/((pi/2) - (a + b)) - 1
        
        #computes the vertical angle of the point with its predecessor
        a <- 0
        if(rangeImage[y-1, x] != 0 && rangeImage[y, x] != 0)
          a <- abs(atan((factor*rangeImage[y,x] - factor*rangeImage[y-1, x])/(1)))
        
        #computes the vertical angle of the point with its successor
        b <- 0
        if(rangeImage[y+1, x] != 0 && rangeImage[y, x] != 0)
          b <- abs(atan((factor*rangeImage[y,x] - factor*rangeImage[y+1, x])/(-1)))
        
        #computes the vertical curvature
        curvesV <- 0
        if(a != 0 && b != 0)
          curvesV <- pi/(pi - (a + b)) - 1
        else
          curvesV <- (pi/2)/((pi/2) - (a + b)) - 1
        
        #computes the total curvature, that is a mean between the horizontal and vertical curvatures
        #ponderated by a factor and the weights
        curves[y,x] <- (curvesH + curvesV)/2 * factor * weightImage[y,x]
        
        #if it is to show the progress, do it
        if(printProgress)
          cat((y * col + x) / (col * row) * 100, "%\n")
      }
    }
  }
  else{
    
    #for each inner point of the image, computes the curvature
    for(y in 2:(row - 1)){
      for(x in 2:(col - 1)){
        
        #computes the horizontal angle of the point with its predecessor
        a <- 0
        if(rangeImage[y, x-1] != 0 && rangeImage[y, x] != 0)
          a <- abs(atan((factor*rangeImage[y,x] - factor*rangeImage[y, x-1])/(1)))
        
        #computes the horizontal angle of the point with its successor
        b <- 0
        if(rangeImage[y, x+1] != 0 && rangeImage[y, x] != 0)
          b <- abs(atan((factor*rangeImage[y,x] - factor*rangeImage[y, x+1])/(-1)))
        
        #computes the horizontal curvature
        curvesH <- 0
        if(a != 0 && b != 0)
          curvesH <- pi/(pi - (a + b)) - 1
        else
          curvesH <- (pi/2)/((pi/2) - (a + b)) - 1
        
        #computes the vertical angle of the point with its predecessor
        a <- 0
        if(rangeImage[y-1, x] != 0 && rangeImage[y, x] != 0)
          a <- abs(atan((factor*rangeImage[y,x] - factor*rangeImage[y-1, x])/(1)))
        
        #computes the vertical angle of the point with its successor
        b <- 0
        if(rangeImage[y+1, x] != 0 && rangeImage[y, x] != 0)
          b <- abs(atan((factor*rangeImage[y,x] - factor*rangeImage[y+1, x])/(-1)))
        
        #computes the vertical curvature
        curvesV <- 0
        if(a != 0 && b != 0)
          curvesV <- pi/(pi - (a + b)) - 1
        else
          curvesV <- (pi/2)/((pi/2) - (a + b)) - 1
        
        #computes the total curvature, that is a mean between the horizontal and vertical curvatures
        #ponderated by a factor and the weights
        curves[y,x] <- (curvesH + curvesV)/2 * factor
        
        #if it is to show the progress, do it
        if(printProgress)
          cat((y * col + x) / (col * row) * 100, "%\n")
      }
    }
  }
  
  #if the curvature is to be saved into a file
  if(file != ""){
    
    #gets the biggest curvature value
    #maxC <- max(curves)
    #curves <- curves/maxC
    
    #save the curvature image
    my.writeJPEG(curves, file, quality=1)
    write(t(curves), concatenate(c(file, ".dat")), ncolumns=col)
  }
  
  (curves)
}

# Generates an expression weight image.
# input:
#   rangeImage = either a path to the image file or a 2D matrix containing the image
#   ldmk = a list with the 9 landmarks
#   file = the path where the weight image shall be saved
#   printProgress = a boolean value that tells whether the progress must be showed
# output:
#   a 2D matrix containing the weight image
weightedRegions <- function(rangeImage, ldmk, file="", printProgress=FALSE){
  
  #if a path is given,
  if(is.character(rangeImage))
    #reads the image from the file
    rangeImage <- readImageData(rangeImage)
  
  #initiates the weight image
  weights <- matrix(rep(0, length(rangeImage[,1]) * length(rangeImage[1,])), ncol=length(rangeImage[1,]))
  
  #forms the region 1 (nose)
  line1 <- findLine(ldmk[["lefteye"]][-3], ldmk[["nasion"]][-3])
  line2 <- findLine(ldmk[["nasion"]][-3], ldmk[["righteye"]][-3])
  line3 <- findLine(ldmk[["righteye"]][-3], ldmk[["rightalare"]][-3], TRUE)
  line4 <- findLine(ldmk[["rightalare"]][-3], ldmk[["subnasale"]][-3])
  line5 <- findLine(ldmk[["subnasale"]][-3], ldmk[["leftalare"]][-3])
  line6 <- findLine(ldmk[["leftalare"]][-3], ldmk[["lefteye"]][-3], TRUE)
  
  #forms the region 2 (left eye)
  h <- dist(rbind(ldmk[["nasion"]][-3], ldmk[["nose"]][-3]))[1]
  line8 <- findLine(ldmk[["nasion"]][-3], ldmk[["subnasale"]][-3], TRUE)
  
  line85 <- line8
  if(line8$a != 0)
    line85 <- findLine(ldmk[["nasion"]][-3], ldmk[["subnasale"]][-3])
  
  line10 <- list(a = line8$a, b = line8$b - (3 * h/2), inverted = TRUE)
  line7 <- findPerpendicularLine(line85, pointAtDistance(line85, ldmk[["nasion"]][-3], h, "y-"))
  line9 <- list(a = line7$a, b = line7$b + (3 * h/2), inverted = line7$inverted)
  
  #gets the lines of region 3 (right eye)
  line11 <- list(a = line8$a, b = line8$b + (3* h/2), inverted = TRUE)
  
  #gets the lines of region 4 (right cheeck)
  y <- 0
  if(ldmk[["leftmouth"]][2] > ldmk[["rightmouth"]][2])
    y <- abs(ldmk[["leftmouth"]][2] - ldmk[["rightmouth"]][2])/2 + ldmk[["rightmouth"]][2]
  else
    y <- abs(ldmk[["leftmouth"]][2] - ldmk[["rightmouth"]][2])/2 + ldmk[["leftmouth"]][2]
    
  d <- dist(rbind(ldmk[["subnasale"]][-3], matrix(c(ldmk[["subnasale"]][1], y), nrow=1)))[1]
  
  line12 <- findLine(ldmk[["leftalare"]][-3], c(ldmk[["leftmouth"]][1] - (d/4), ldmk[["leftmouth"]][2]))
  line13 <- findPerpendicularLine(line85, ldmk[["subnasale"]][-3])
  
  #gets the lines of the region 5 (right cheeck)
  line14 <- findLine(ldmk[["rightalare"]][-3], c(ldmk[["rightmouth"]][1] + (d/4), ldmk[["rightmouth"]][2]))
  
  #gets the lines of the region 6 (mouth)
  line15 <- findLine(c(ldmk[["rightmouth"]][1] + (d/4), ldmk[["rightmouth"]][2]), c(ldmk[["rightmouth"]][1], ldmk[["rightmouth"]][2] + (d/2)))
  line16 <- list(a = line7$a, b = line13$b + 2*d, inverted = line7$inverted)
  line17 <- findLine(c(ldmk[["leftmouth"]][1], ldmk[["leftmouth"]][2] + (d/2)), c(ldmk[["leftmouth"]][1] - (d/4), ldmk[["leftmouth"]][2]))
  line20 <- findLine(c(ldmk[["subnasale"]][2], ldmk[["subnasale"]][1]), c(ldmk[["subnasale"]][2] + d, ldmk[["subnasale"]][1]))
  
  #gets the lines of the region 9 (chin)
  line18 <- findLine(c(ldmk[["rightmouth"]][2], ldmk[["rightmouth"]][1]), c(ldmk[["rightmouth"]][2] + (d/2), ldmk[["rightmouth"]][1]))
  line19 <- findLine(c(ldmk[["leftmouth"]][2], ldmk[["leftmouth"]][1]), c(ldmk[["leftmouth"]][2] + (d/2), ldmk[["leftmouth"]][1]))
  
  #creates an image with the lines
  #linesImage <- lines2image(list(line7, line8, line1),
  #                          progressFlag=TRUE)
  #writeImage(t(linesImage), "Images/lines.jpg", "jpg")
  
  #classifies each pixel of the range image
  for(y in 1:length(rangeImage[,1])){
    for(x in 1:length(rangeImage[1,])){
      
      #cat("x = ", x, "; y = ", y, "\n")
      
      #checks if the point belongs to the 1st region
      if((y >= appLinear(line1, x)) && (y >= appLinear(line2, x)) &&
           (x <= appLinear(line3, y)) && (y <= appLinear(line4, x)) &&
           (y <= appLinear(line5, x)) && (x >= appLinear(line6, y))){
        
        weights[y, x] = 1
      }
      else{
        
        #cat("x = ", x, ", y = ", y, "\n")
        
        #checks if the point belongs to the 2nd region
        if(((y >= appLinear(line7, x)) && (x <= appLinear(line8, y)) &&
              (y <= appLinear(line1, x)) && (x >= appLinear(line6, y))) ||
             ((y >= appLinear(line7, x)) && (x <= appLinear(line6, y)) &&
                (y <= appLinear(line9, x)) && (x >= appLinear(line10, y)))){
          
          weights[y, x] = 1/10
        }
        else{
          #checks if the point belongs to the 3rd region
          if(((y >= appLinear(line7, x)) && (x <= appLinear(line11, y)) &&
                (y <= appLinear(line9, x)) && (x >= appLinear(line3, y))) ||
               ((y >= appLinear(line7, x)) && (x <= appLinear(line3, y)) &&
                  (y <= appLinear(line2, x)) && (x >= appLinear(line8, y)))){
            
            weights[y, x] <- 1/10
          }
          else{
            #checks if the point belongs to the 4th region
            if((y >= appLinear(line9, x)) && (x <= appLinear(line6, y)) &&
                 (y <= appLinear(line12, x)) && (y <= appLinear(line13, x)) &&
                 (x >= appLinear(line10, y))){
              
              weights[y, x] = 1/2
            }
            else{
              #checks if the point belongs to the 5th region
              if((y >= appLinear(line9, x)) && (x <= appLinear(line11, y)) &&
                   (y <= appLinear(line13, x)) && (y <= appLinear(line14, x)) &&
                   (x >= appLinear(line3, y))){
                
                weights[y, x] = 1/2
              }
              else{
                #checks if the point belongs to the 6th region
                if(((y >= appLinear(line5, x)) && (x <= appLinear(line20, y)) &&
                      (y <= appLinear(line16, x)) && (y <= appLinear(line17, x)) &&
                      (y >= appLinear(line12, x))) ||
                     ((y >= appLinear(line4, x)) && (y >= appLinear(line14, x)) &&
                        (y <= appLinear(line15, x)) && (y <= appLinear(line16, x)) &&
                        (x >= appLinear(line20, y)))){
                  
                  weights[y, x] = 1/10
                }
                else{
                  #checks if the point belongs to the 7th region
                  if(((y >= appLinear(line13, x)) && (y <= appLinear(line12, x)) &&
                        (x >= appLinear(line10, y))) ||
                       ((y >= appLinear(line17, x)) && (y <= appLinear(line16, x)) &&
                          (x >= appLinear(line10, y)))){
                    
                    weights[y, x] = 1/10
                  }
                  else{
                    #checks if the point belongs to the 8th region
                    if(((y >= appLinear(line13, x)) && (x <= appLinear(line11, y)) &&
                          (y <= appLinear(line14, x))) ||
                         ((y >= appLinear(line15, x)) && (x <= appLinear(line11, y)) &&
                            (y <= appLinear(line16, x)))){
                      
                      weights[y, x] = 1/10
                    }
                    else{
                      #checks if the point belongs to the 9th region
                      if((x >= appYLinear(line12, y)) && (x <= appYLinear(line14, y)) &&
                           (y >= appLinear(line16, x))){
                        
                        weights[y, x] = 1/2
                      }
                      else{
                        #then it belongs to the 10th region
                        weights[y, x] = 1
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      
      if(printProgress)
        cat((y * length(rangeImage[1,]) + x) / (length(rangeImage[,1]) * length(rangeImage[1,])) * 100, "%\n")
    }
  }
  
  if(file != ""){
    
    writeJPEG(weights, file, quality=1)
    write(t(weights), concatenate(c(file, ".dat")), ncolumns=length(rangeImage[1,]))
  }
  
  (weights)
}

# Finds the perpendicular line of a given line that passes
# through a given point.
# intput:
#   line = a list with the paramenters 'a' and 'b' that describes a linear function
#          a linear function is given by f(x) = y = ax + b
#   point = a vector with 2 elements, 'x' and 'y' coordinate of a 2D point
# output:
#   a line represented as list, just like the input 'line'
findPerpendicularLine <- function(line, point){
  
  #initializes the angular coefficient 'a'
  a <- 0
  
  #if the line isn't paralell to the domain axis,...
  if(line$a != 0){
    
    #computes the angle of the line
    alpha <- atan(line$a)
    #computes the angle of the perpendicular line
    beta <- pi - (pi/2 + abs(alpha))
    
    #computes the angular coefficient of the perpendicular line
    #if alpha is positive, then beta should be negative
    if(alpha >= 0)
      a <- tan(-beta)
    #else, beta should be negative
    else
      a <- tan(beta)
    
    #computes the parameter 'b'
    b <- point[2] - a * point[1]
    
    #returns the line
    (list(a = a, b = b, inverted = line$inverted))
  }
  #if the line is parallel to the domain axis,...
  else{
    
    #if the line is inverted, ...
    if(line$inverted)
      #then, the parameter 'b' equals 'y' and it's not inverted
      (list(a = 0, b = point[2], inverted = FALSE))
    else
      #else, the parameter 'b' equals 'x' and it's inverted
      (list(a = 0, b = point[1], inverted = TRUE))
  }
}

# Given a line, a point in this line, a distance and a direction on this line,
# it finds the point in this line which is that distant from the 1st point.
# input:
#   line = a list containing the 'a' and 'b'parameters of a linear function
#          a linear function is given by f(x) = y = ax + b
#   point = a vector with 2 elements, the 2 coordinates, 'x' and 'y'
#   distance = the distance between the reference point and the one to be found
#   direction = a string specifying if the point has to be on the right, left,
#               upper, lower of the reference point. e.g. "y+" to upper,
#               "y-" to lower, "x+" to the right, "x-" to the left
# output:
#   a vector with 2 elements, the 2 coordiantes, 'x' and 'y' of the wanted point
pointAtDistance <- function(line, point, distance, direction){
  
  #initializes the point's coordinates
  x <- 0
  y <- 0
  
  #if the line isn't parallel to the domain axis
  if(line$a != 0){
    #computes the point's coordinates
    if(direction == "x+"){
      x <- point[1] + cos(atan(abs(line$a))) * distance
      y <- appLinear(line, x)
    }
    else
      if(direction == "x-"){
        x <- point[1] - cos(atan(abs(line$a))) * distance
        y <- appLinear(line, x)
      }
    else
      if(direction == "y+"){
        y <- point[2] + sin(atan(abs(line$a))) * distance
        x <- appYLinear(line, y)
      }
    else
      if(direction == "y-"){
        y <- point[2] - sin(atan(abs(line$a))) * distance
        x <- appYLinear(line, y)
      }
  }
  else{# if the line is parallel to the domain axis
    #and it is not inverted
    if(!line$inverted){
      #the y is the same
      y <- point[2]
      #the x is addition or subtraction of the distance,
      #depending on the direction
      if(direction == "x+" || direction == "y+")
        x <- point[1] + distance
      else
        x <- point[1] - distance
    }
    else{ #if the line is inverted
      #the x is the same
      x <- point[1]
      #the y is the addition or subtraction of the distance,
      #depending on the direction
      if(direction == "x+" || direction == "y+")
        y <- point[2] + distance
      else
        y <- point[2] - distance
    }
  }
  #returns the point  
  c(x,y)
}

# Applies an image value 'y' into a linear function 'line'
# and returns the corresponding domain 'x'.
# consider a linear function as f(x) = y = ax + b
# input:
#   line = a list containing 'a' and 'b' parameters of a linear function
#   y = the linear function's image
# output:
#   a value 'x', the domain of the given image
#   if 'a' = 0, then it will return Infinite
appYLinear <- function(line, y){
  
  #if the line isn't parallel to the domain axis
  if(line$a != 0)
    x <- (y - line$b)/line$a
  else #if it is,
    x <- Inf #return infinite
  (x)
}

# Applies a input x into a given linear function
# and returns the corresponding y.
# consider a linear function as f(x) = y = ax + b
# input:
#   line = a list containing 'a' and 'b' parameters of a linear function
#   x = the linear function's argument
# output:
#   a value 'y', the resulting of applying 'x' to the function 'line'
appLinear <- function(line, x){
  
  y <- line$a * x + line$b
  (y)
}

# Finds the line given by 2 points
# If the points have the same 'x' (2nd) coordinate, then the
# axis are inverted.
# input:
#   p1 = the 1st 2D point
#   p2 = the 2nd 2D point
#   inverted = whether the axis should be inverted despite the line
# output:
#   a list with the 3 parameters: 'a' and 'b', discribing the linear function
#   and a 3rd paramater telling whether the axis are inverted
#   a linear function is given by f(x) = y = ax + b
findLine <- function(p1, p2, inverted = FALSE){
  
  a <- 0
  b <- 0
  
  #If the x coordinate aren't equal for both points
  if(p1[1] - p2[1] != 0 && inverted == FALSE){
    #compute the parameters normally
    a <- (p1[2] - p2[2])/(p1[1] - p2[1])
    b <- p2[2] - p2[1] * a
  }
  else{ #if the are the same, invert the image and domain
    
    a <- (p1[1] - p2[1])/(p1[2] - p2[2])
    b <- p2[1] - p2[2] * a
    inverted <- TRUE
  }
  #returns the line
  (list(a = a, b = b, inverted = inverted))
}

# Plots in black the landmarks over the range image.
# input:
#   img = either a path to the range image or a 2D matrix containing the image
#   landmarks = either a path to the landmarks or a list with the landmarks coordinates
#   toFile = a path where the resulting image shall be saved
# output:
#   the range image with the landmarks "holes" (in black)
seeLandmakrs <- function(img, landmarks, toFile=""){
  
  #if a path is given for the image
  if(is.character(img))
    #reads it from the file
    img <- readImageData(img)
  
  #if a path is given for the landmarks
  if(is.character(landmarks))
    #reads it from the file
    landmarks <- readLandmark(landmarks)
  
  #for each landmark
  for(i in 1:length(landmarks)){
    #put a 'hole' (black point) into the image
    img[landmarks[[i]][2], landmarks[[i]][1]] <- 0
  }
  
  #if the image shall be saved, do it
  if(toFile != "")
    my.writeJPEG(img, toFile, quality=1)
  
  (img)
}

# Checks whether a given point belongs to a face from a range image
# input:
#   image = either a path to the image or a 2D matrix containing the image
#   row = the row coordinate of the point
#   col = the col coordinate of the point
# output:
#   a boolean value; TRUE if the point belongs to the face, FALSE, otherwise
isFace <- function(image, row, col){
  
  #if there is no point different from 0 to the right of the point, returns FALSE
  if(length(which(image[row, (col + 1):length(image[1,])] != 0)) == 0)
    (FALSE)
  #if there is no point different from 0 to the left of the point, returns FALSE
  else if(length(which(image[row, 1:(col-1)] != 0)) == 0)
    (FALSE)
  #if there is no point different from 0 above the point, returns FALSE
  else if(length(which(image[1:(row-1), col] != 0)) == 0)
    (FALSE)
  #if there is no point different from 0 below the point, returns FALSE
  else if(length(which(image[(row+1):length(image[,1]), col] != 0)) == 0)
    (FALSE)
  
  #otherwise, returns TRUE
  else (TRUE)
}

# Checks whether a given point belongs to a face from a range image
# input:
#   image = either a path to the image or a 2D matrix containing the image
#   row = the row coordinate of the point
#   col = the col coordinate of the point
# output:
#   a boolean value; TRUE if the point belongs to the face, FALSE, otherwise
isFace2 <- function(image, row, col){
  
  #if there is any point different from 0 to the right and left of the point, returns TRUE
  if((length(which(image[row, (col + 1):length(image[1,])] != 0)) != 0) &&
       (length(which(image[row, 1:(col-1)] != 0)) != 0))
    (TRUE)
  #if there is any point different from 0 above and below the point, returns TRUE
  else if((length(which(image[1:(row-1), col] != 0)) != 0) &&
            (length(which(image[(row+1):length(image[,1]), col] != 0)) != 0))
    (TRUE)
  
  #otherwise, returns TRUE
  else (FALSE)
}

# Fills the holes of a surface by analysing the 3x3 neighbourhood
# of the hole and replaces it with the mode or the mean of its neighbourhood
# input:
#   image = either the path to the jpg range image or a 2D matrix containing the image
#   toFile = (optional) the path to the new image
#   maxIter = the maximum number of iterations
#   printProgress = whether the progress must be shown
# output:
#   a 2D matrix with the new range image without holes
holeCorrection <- function(image, fileToSave = "", maxIter=1, printProgress=FALSE){
  
  #if a path for the image is given,
  if(is.character(image))
    #reads it from the file
    image <- readImageData(image)
  
  #initializes the variable that controls the number of holes that still persists
  counter <- 1
  #initializes the variable that controls the number of iterations completed
  times <- 0
  
  #while still there is left at least a hole and the maximum number of iterations
  #hasn't been reached, does...
  while(counter > 0 && times <= maxIter){
  
    counter <- 0 #the number of holes found starts with zero
    times <- times + 1 #increases the number of iterations completed
    
    #searches for the image, but one-pixel border
    for(row in 2:(length(image[,1])-1)){
      for(col in 2:(length(image[1,])-1)){
        
        #if this point is black (hole | 0) and it belongs to the face, ...
        if(image[row, col] == 0 && isFace(image, row, col)){
          
          #retrieves the 8 neighbours
          neighbours <- rep(0, 8)
          n <- 0
          for(i in -1:1){
            for(j in -1:1){
              
              if(i != 0 || j != 0){
                n <- n + 1
                neighbours[n] <- image[row + i, col + j]
              }
            }
          }
          
          #computes the mode of the neighbours, or it doesn't exist, the mean
          med <- getMode(neighbours[which(neighbours != 0)])
          
          #if all neighbours were 0s, then med = 0
          if(is.null(med) || is.na(med)){
            med = 0
            #cat("\nnull\n")
          }
          
          #if the mode or mean is different than zero, we have a hole filled
          if(med != 0){
            counter <- counter + 1
          }
          
          #cat("row = ", row, "col = ", col, " med = ", med, "\n")
          
#           if(col == 23 && row == 17){
#             cat("med = ", med, "\nneighbours:\n")
#             for(i in which(neighbours != 0))
#               cat(neighbours[i], ", ")
#             cat("\n")
#           }
          
          #the image receives the mode or mean
          image[row, col] <- med
        }
        
        #cat((row * col)/((length(image[,1]) - 2) * (length(image[1,]) - 2)) * 100, "%\n")
      }
    }
    
    if(printProgress)
      cat("iteration ", times, "\n")
  }
   
  if(fileToSave != ""){
    my.writeJPEG(image, fileToSave, quality=1)
    write(t(image), paste(fileToSave, ".dat", sep=""), ncolumns=length(image[1,]))
  }
  
  (image)
}

outlierCorrection <- function(file, fileToSave = "", maxIter=1){
  
  if(is.character(file))
    image <- readImageData(file)
  else
    image <- file
  
  counter <- 1
  times <- 0
  
  while(counter > 0 && times <= maxIter){
    
    counter <- 0
    times <- times + 1
    d <- 0
    
    for(row in 2:(length(image[,1])-1)){
      for(col in 2:(length(image[1,])-1)){
        
        if(isFace(image, row, col)){
          
          neighbours <- rep(0, 8)
          n <- 0
          maxDiff <- 0
          for(i in -1:1){
            for(j in -1:1){
              
              if(i != 0 || j != 0){
                
                n <- n + 1
                neighbours[n] <- image[row + i, col + j]
                
                if(is.logical(maxDiff > abs(image[row, col] - image[row + i, col + j])))
                  if(maxDiff > abs(image[row, col] - image[row + i, col + j]))
                    maxDiff <- abs(image[row, col] - image[row + i, col + j])
              }
            }
          }
          
          if(d < maxDiff){
            
            d <- maxDiff
          }
          
          if(maxDiff >= 0.1){
            counter <- counter + 1
            #replaces the value of the pixel by the mode of the neighbours
            image[row, col] <- getMode(neighbours[which(neighbours != 0)])
          }
        }
      }
    }
    #cat(times, ", difference = ", d, "\n")
  }
  
  if(fileToSave != ""){
    
    my.writeJPEG(image, fileToSave, quality=1)
    write(t(image), concatenate(c(fileToSave, ".dat")), ncolumns=col)
  }
  
  (image)
}

setConnectedComponentLabel <- function(img, label, point, labels){
  
  col <- length(img[1,])
  row <- length(img[,1])
  
  line <- list()
  line[[1]] <- point
  tail <- 2
  
  while(length(line) > 0){
    
    point <- line[[1]]
    
    if(point[1] >= 1 && point[1] <= row && point[2] >= 1 && point[2] <= col &&
         labels[point[1], point[2]] == 0 && img[point[1], point[2]] != 0){
      
      labels[point[1], point[2]] <- label
      
      line[[tail]] <- c(point[1] + 1, point[2])
      line[[tail + 1]] <- c(point[1], point[2] + 1)
      line[[tail + 2]] <- c(point[1] - 1, point[2])
      line[[tail + 3]] <- c(point[1], point[2] - 1)
      tail <- tail + 4
    }
    
    line <- line[-1]
    tail <- tail - 1
  }
  
  (labels)
}

# Removes all connected component from an image
# other than the biggest one.
# The connected component consider 4-size neighbourhood
connectedComponentAnalysis <- function(img, file="", progress=FALSE){
  
  if(is.character(img))
    img <- readImageData(img)
  
  col <- length(img[1,])
  row <- length(img[,1])
  
  labels <- matrix(rep(0, col*row), ncol=col)
  nComponent <- 0
  
  for(i in 2:(row-1)){
    for(j in 2:(col-1)){
      
      if(img[i,j] != 0 && labels[i,j] == 0){
        
        nComponent <- nComponent + 1
          
        labels <- setConnectedComponentLabel(img, nComponent, c(i,j), labels)
      }
      
      if(progress)
        cat((i*col + j) * 100/(col*row), "%\n", i, ", ", j, "\n")
    }
  }
  
  faceComponent <- getMode(labels[which(labels != 0)])
  
  img[which(labels != faceComponent)] <- 0
  
  if(file != ""){
    
    my.writeJPEG(img, file, quality=1)
    #labels <- labels/nComponent
    #writeImage(t(labels), paste(file, ".labels", sep=""), "jpg")
    write(t(img), concatenate(c(file, ".dat")), ncolumns=col)
  }
  
  (img)
}

outlierCorrection2 <- function(img, curvature, file="", correction=0, progress=FALSE){
  
  if(is.character(img))
    img <- readImageData(img)
  
  if(is.character(curvature))
    curvature <- readImageData(curvature)
  
  col <- length(img[1,])
  row <- length(img[,1])
  
  for(i in 1:col){
    imgVarX <- getXvariation(matrix(c(1:row, img[,i]), ncol=2))
    
    threshold <- 2*mean(curvature[which(curvature[,i] != 0),i])
    if(!is.na(threshold) && imgVarX$max > imgVarX$min){
      
      above <- which(curvature[(imgVarX$min):imgVarX$max,i] > threshold)
      below <- which(curvature[(imgVarX$min):imgVarX$max,i] <= threshold)
      
      curve <- 0
      
      if(length(above) > 1)
        for(j in 2:length(above)){
          
          diff <- above[j] - above[j-1]
          if(diff > 1 && diff <= 5){
            
            for(k in 1:(diff-1)){
              
              above <- c(above, above[j-1] + k)
              below <- below[-which(below == (above[j-1] + k))]
            }
          }
        }
      
      if(length(above) > 0){
        
        xs <- below + imgVarX$min-1
        ys <- img[below + imgVarX$min-1,i]
        
        if(above[1] == 1){
        
           belowX <- imgVarX$min - 1        
#           line <- linearInterpolation(matrix(c(below[1:7] + imgVarX$min - 1, img[below[1:7] + imgVarX$min - 1,i]), ncol=2))
#           belowY <- appLinear(line, belowX)
#           
#           xs <- c(belowX, xs)
#           ys <- c(belowY, ys)
          
          xs <- c(belowX, xs)
          ys <- c(0, ys)
        }
        
        if(max(above) == imgVarX$max - imgVarX$min +1){
          
          belowX <- imgVarX$max + 1
#           line <- linearInterpolation(matrix(c(below[(length(below)-5):length(below)] + imgVarX$min - 1, img[below[(length(below)-5):length(below)] + imgVarX$min - 1,i]), ncol=2))
#           belowY <- appLinear(line, belowX)
#           
#           xs <- c(xs, belowX)
#           ys <- c(ys, belowY)
          
          xs <- c(xs, belowX)
          ys <- c(ys, 0)
        }
        
        curve <- stinterp(xs, ys, above + imgVarX$min-1, method="sti")
        img[curve$x,i] <- curve$y
      }
      #else
      #  curve <- stinterp(below + imgVarX$min-1, img[below + imgVarX$min-1,i], above + imgVarX$min-1, method="sti")
      
      #img[curve$x,i] <- curve$y
    }
    
    if(progress)
      cat(i * 100/col, "% ", i, "\n")
  }
  
  if(correction != 0)
    img <- horizontalCorrection(img, correction)
  
  if(file != ""){
    
    my.writeJPEG(img, file, quality=1)
    write(t(img), paste(file, ".dat", sep=""), ncolumns=col)
  }
    
  (img)
}

horizontalCorrection <- function(img, factor=0.4, progress=FALSE){
  
  col <- length(img[1,])
  row <- length(img[,1])
  
  for(i in 1:row){
    for(j in 1:col){
      
      if(img[i,j] != 0){
        
        if(j == 1){
          
          if((img[i,j+1] - img[i,j]) / 1 >= factor)
            img[i,j] <- getMode(c(img[i,j], img[i,j+1], img[i,j+2]))
        }
        else if(j == col){
          
          if((img[i,j-1] - img[i,j]) / 1 >= factor)
            img[i,j] <- getMode(c(img[i,j], img[i,j-1], img[i,j-2]))
        }
        else{
          
          if((img[i,j+1] - img[i,j]) / 1 >= factor || (img[i,j-1] - img[i,j]) / 1 >= factor)
            img[i,j] <- getMode(c(img[i,j-1], img[i,j], img[i,j+1]))
        }
      }
      
      if(progress)
        cat((i*col + j)/(col*row)*100, "%\n")
    }
  }
  
  (img)
}

outlierCorrection3 <- function(img, curvature, file="", progress=FALSE, maxORmin="max"){
  
  if(is.character(img))
    img <- readImageData(img)
  
  if(is.character(curvature))
    curvature <- readImageData(curvature)
  
  col <- length(img[1,])
  row <- length(img[,1])
  
  colInterpolation <- img
  rowInterpolation <- img
  
  for(i in 1:col){
      
    imgVarX <- getXvariation(matrix(c(1:row, img[,i]), ncol=2))
    
    threshold <- 2*mean(curvature[which(curvature[,i] != 0),i])
    if(!is.na(threshold)){
      
      above <- which(curvature[(imgVarX$min):imgVarX$max,i] > threshold)
      below <- which(curvature[(imgVarX$min):imgVarX$max,i] <= threshold)
      
      curve <- 0
      
      if(length(above) > 1)
        for(j in 2:length(above)){
          
          diff <- above[j] - above[j-1]
          if(diff > 1 && diff <= 5){
            
            for(k in 1:(diff-1)){
              
              above <- c(above, above[j-1] + k)
              below <- below[-which(below == (above[j-1] + k))]
            }
          }
        }
      
      if(length(above) > 0){
        
        xs <- below + imgVarX$min-1
        ys <- img[below + imgVarX$min-1,i]
        
        if(above[1] == 1){
          
          belowX <- imgVarX$min - 1        
          line <- linearInterpolation(matrix(c(below[1:7] + imgVarX$min - 1, img[below[1:7] + imgVarX$min - 1,i]), ncol=2))
          belowY <- appLinear(line, belowX)
          
          xs <- c(belowX, xs)
          ys <- c(belowY, ys)
        }
        
        if(max(above) == imgVarX$max - imgVarX$min +1){
          
          belowX <- imgVarX$max + 1
          line <- linearInterpolation(matrix(c(below[(length(below)-5):length(below)] + imgVarX$min - 1, img[below[(length(below)-5):length(below)] + imgVarX$min - 1,i]), ncol=2))
          belowY <- appLinear(line, belowX)
          
          xs <- c(xs, belowX)
          ys <- c(ys, belowY)
        }
        
        curve <- stinterp(xs, ys, above + imgVarX$min-1, method="sti")
      }
      else
        curve <- stinterp(below + imgVarX$min-1, img[below + imgVarX$min-1,i], above + imgVarX$min-1, method="sti")
      
      colInterpolation[curve$x,i] <- curve$y
    }
    
    if(progress)
      cat(i * 100/(col + row), "% ", i, "\n")
  }
  
  for(i in 1:row){
    
    imgVarX <- getXvariation(matrix(c(1:col, img[i,]), ncol=2))
    
    threshold <- 2*mean(curvature[i,which(curvature[i,] != 0)])
    if(!is.na(threshold)){
      
      above <- which(curvature[i,imgVarX$min:imgVarX$max] > threshold)
      below <- which(curvature[i,imgVarX$min:imgVarX$max] <= threshold)
      
      curve <- 0
      
      if(length(above) > 1)
        for(j in 2:length(above)){
          
          diff <- above[j] - above[j-1]
          if(diff > 1 && diff <= 5){
            
            for(k in 1:(diff-1)){
              
              above <- c(above, above[j-1] + k)
              below <- below[-which(below == (above[j-1] + k))]
            }
          }
        }
      
      if(length(above) > 0){
        
        xs <- below + imgVarX$min-1
        ys <- img[i,below + imgVarX$min-1]
        
        if(above[1] == 1){
          
          belowX <- imgVarX$min - 1        
          line <- linearInterpolation(matrix(c(below[1:7] + imgVarX$min - 1, img[i,below[1:7] + imgVarX$min - 1]), ncol=2))
          belowY <- appLinear(line, belowX)
          
          xs <- c(belowX, xs)
          ys <- c(belowY, ys)
        }
        
        if(max(above) == imgVarX$max - imgVarX$min +1){
          
          belowX <- imgVarX$max + 1
          line <- linearInterpolation(matrix(c(below[(length(below)-5):length(below)] + imgVarX$min - 1, img[i,below[(length(below)-5):length(below)] + imgVarX$min - 1]), ncol=2))
          belowY <- appLinear(line, belowX)
          
          xs <- c(xs, belowX)
          ys <- c(ys, belowY)
        }
        
        curve <- stinterp(xs, ys, above + imgVarX$min-1, method="sti")
      }
      else
        curve <- stinterp(below + imgVarX$min-1, img[i,below + imgVarX$min-1], above + imgVarX$min-1, method="sti")
      
      rowInterpolation[i,curve$x] <- curve$y
    }
    
    if(progress)
      cat((i + col) * 100/(col + row), "% ", i, "\n")
  }
  
  img <- my.maxORmin(rowInterpolation, colInterpolation, maxORmin)
  
  if(file != ""){
    
    my.writeJPEG(img, file, quality=1)
    write(t(img), paste(file, ".dat", sep=""), ncolumns=col)
  }
  
  (img)
}

#Resample the face points, so the final amount of
#points is the one determined. 
#The k-means algorithm is used in such a task.
#input:
#   points = the 3D face points
#   number = the final amount of points to be resampled
#output:
#   a matrix containing the 3D points resampled
resample <- function(points, number){
  
  sample <- kmeans(points, number, nstart=3, iter.max=20)
  
  (sample$centers)
}

#Rotates the face points, so the eyes will be aligned
#in xy and xz planes
#input:
#   points = the 3D points of the face in a matrix format rows = points
#   eyes = a list with the landmark coordinates
#output:
#   a list with "points", a matrix contaning all points rotated on both planes
#   and with "ldmk", a list with the landmark updated coordinates
alignFace <- function(points, eyes){
  
  #computes the xz angle
  a <- (eyes[["righteye"]][1,"z"] - eyes[["lefteye"]][1,"z"])/(eyes[["righteye"]][1,"x"] - eyes[["lefteye"]][1,"x"])
  xzAngle <- atan(a)
  
  #measures the translation factors
  originFactor <- matrix(c(-eyes[["nasion"]][1,"x"], -eyes[["nasion"]][1,"y"], -eyes[["nasion"]][1,"z"]), nrow=1, dimnames=list(c(), c("x", "y", "z")))
  
  #builds the rotation matrix
  yMatrix <- matrix(c(cos(xzAngle), 0, -sin(xzAngle), 0, 1, 0, sin(xzAngle), 0, cos(xzAngle)), nrow=3)
  
  #rotates all points on XY plane
  for(i in 1:length(points[,1])){
    
    newPoint <- yMatrix %*% matrix(points[i,] + originFactor[1,], ncol=1)
    points[i,"x"] <- newPoint[1,1] - originFactor[1,"x"]
    points[i,"y"] <- newPoint[2,1] - originFactor[1,"y"]
    points[i,"z"] <- newPoint[3,1] - originFactor[1,"z"]
    
    #prints the progress
    #cat(i/(2*length(points[,1])) * 100, "%\n")
  }
  
  #updates the reference points
  eyes[["lefteye"]] <- t(yMatrix %*% t(eyes[["lefteye"]] + originFactor)) - originFactor
  eyes[["righteye"]] <- t(yMatrix %*% t(eyes[["righteye"]] + originFactor)) - originFactor
  eyes[["nose"]] <- t(yMatrix %*% t(eyes[["nose"]] + originFactor)) - originFactor
  eyes[["nasion"]] <- t(yMatrix %*% t(eyes[["nasion"]] + originFactor)) - originFactor
  eyes[["leftalare"]] <- t(yMatrix %*% t(eyes[["leftalare"]] + originFactor)) - originFactor
  eyes[["rightalare"]] <- t(yMatrix %*% t(eyes[["rightalare"]] + originFactor)) - originFactor
  eyes[["subnasale"]] <- t(yMatrix %*% t(eyes[["subnasale"]] + originFactor)) - originFactor
  eyes[["leftmouth"]] <- t(yMatrix %*% t(eyes[["leftmouth"]] + originFactor)) - originFactor
  eyes[["rightmouth"]] <- t(yMatrix %*% t(eyes[["rightmouth"]] + originFactor)) - originFactor
  
  #computes xy angle
  a <- (eyes[["righteye"]][1,"y"] - eyes[["lefteye"]][1,"y"])/(eyes[["righteye"]][1,"x"] - eyes[["lefteye"]][1,"x"])
  xyAngle <- atan(a)
  
  #builds the rotation matrix
  zMatrix <- matrix(c(cos(-xyAngle), sin(-xyAngle), 0, -sin(-xyAngle), cos(-xyAngle), 0, 0, 0, 1), nrow=3)
  
  #rotates all points on XZ plane
  for(i in 1:length(points[,1])){
    
    newPoint <- zMatrix %*% matrix(points[i,] + originFactor[1,], ncol=1)
    points[i,"x"] <- newPoint[1,1] - originFactor[1,"x"]
    points[i,"y"] <- newPoint[2,1] - originFactor[1,"y"]
    points[i,"z"] <- newPoint[3,1] - originFactor[1,"z"]
    
    #prints the progress
    #cat((i + length(points[,1]))/(2*length(points[,1])) * 100, "%\n")
  }
  
  #updates the reference points
  eyes[["lefteye"]] <- t(zMatrix %*% t(eyes[["lefteye"]] + originFactor)) - originFactor
  eyes[["righteye"]] <- t(zMatrix %*% t(eyes[["righteye"]] + originFactor)) - originFactor
  eyes[["nose"]] <- t(zMatrix %*% t(eyes[["nose"]] + originFactor)) - originFactor
  eyes[["nasion"]] <- t(zMatrix %*% t(eyes[["nasion"]] + originFactor)) - originFactor
  eyes[["leftalare"]] <- t(zMatrix %*% t(eyes[["leftalare"]] + originFactor)) - originFactor
  eyes[["rightalare"]] <- t(zMatrix %*% t(eyes[["rightalare"]] + originFactor)) - originFactor
  eyes[["subnasale"]] <- t(zMatrix %*% t(eyes[["subnasale"]] + originFactor)) - originFactor
  eyes[["leftmouth"]] <- t(zMatrix %*% t(eyes[["leftmouth"]] + originFactor)) - originFactor
  eyes[["rightmouth"]] <- t(zMatrix %*% t(eyes[["rightmouth"]] + originFactor)) - originFactor
  
  (list(points = points, ldmk = eyes))
}

#Cuts the face removing neck, ears, shoulders, etc,
#input:
# points = The original 3D points
# nose = A matrix containing the coordinate of the nose
# factor = the factor to compute the threshold
#output:
# a matrix containing the cut points
cutTheFaceByNoseLength <- function(points, ldmk){
  
  #gets the number of points
  rows <- length(points[,1])
  
  #initates the vector with the distances between each point to the nose point
  distances <- rep(0, rows)
  
  for(i in 1:rows){
    
    m <- rbind(ldmk[["nose"]], points[i,])
    #computes the distance between the ith point to the nose point
    distances[i] <- dist(m)[1]
    #prints the progress
    #cat(i/rows*100, "%\n")
  }
  
  #computes the threshold distance
  threshold <- dist(rbind(ldmk[["subnasale"]], ldmk[["nasion"]]))[1] * 4/2
  
  #removes all points above the threshold
  newPoints <- points[which(distances < threshold),]
  
  #returns the points within the threshold
  (newPoints)
}

#Cuts the face removing neck, ears, shoulders, etc,
#input:
# points = The original 3D points
# nose = A matrix containing the coordinate of the nose
# factor = the factor to compute the threshold
#output:
# a matrix containing the cut points
cutTheFace <- function(points, nose, factor=0.5, threshold=100){
  
  #gets the number of points
  rows <- length(points[,1])
  
  #initates the vector with the distances between each point to the nose point
  distances <- rep(0, rows)
  
  for(i in 1:rows){
    
    m <- rbind(nose, points[i,])
    #computes the distance between the ith point to the nose point
    distances[i] <- dist(m)[1]
    #prints the progress
    #cat(i/rows*100, "%\n")
  }
  
  #computes the maximum and minimum distances
  maxDist <- max(distances)
  minDist <- min(distances)
  
  #computes the threshold distance
  if(factor != 0)
    threshold <- factor * (maxDist - minDist)
  
  #removes all points above the threshold
  newPoints <- points[which(distances < threshold),]
  
  #returns the points within the threshold
  (newPoints)
}

#Writes the landmarks coordinates into a text file with LDMK extension.
#input:
#   ldmk = a list with a set of named elements containing a matrix 1x3
#          with the landmark coordinate.
#   file = a string with the path to the file to be saved.
#output:
#   none.
writeLandmark <- function(ldmk, file){
  
  text <- paste("# ", file, sep="")
  text <- paste(text, "\n# nb_of_pair key1 value1 key2 value2 ...\n", sep="")
  for(i in 1:length(ldmk)){
    
    text <- paste(text, "2 position ", sep="")
    text <- paste(text, as.character(ldmk[[i]][1]), sep="")
    text <- paste(text, ",", sep="")
    text <- paste(text, as.character(ldmk[[i]][2]), sep="")
    text <- paste(text, ",", sep="")
    text <- paste(text, as.character(ldmk[[i]][3]), sep="")
    text <- paste(text, " label ", sep="")
    text <- paste(text, dimnames(as.matrix(ldmk))[[1]][i], sep="")
    
    if(i != length(ldmk)){
      text <- paste(text, "\n", sep="")
    }
  }
  
  writeLines(text, file)
}

#Reads the coordinates of a set of landmarks.
#input:
#   file = the landmark file path
#output:
#   A list with 9 matrices containing the coordinates of the nose tip,
#   left and right eyes inner point, nose corners and base, mouth corners
#   and the point between eyes (the start of the nose)
readLandmark <- function(file){
  
  landmarks <- list()
  file_content <- readLines(file)
  
  for(i in 1:length(file_content)){
    
    if(strsplit(file_content[i], "[ \t]")[[1]][1] != "#"){
      
      x <- as.numeric(strsplit(strsplit(file_content[i], "[ \t]")[[1]][3], ",")[[1]][1])
      y <- as.numeric(strsplit(strsplit(file_content[i], "[ \t]")[[1]][3], ",")[[1]][2])
      z <- as.numeric(strsplit(strsplit(file_content[i], "[ \t]")[[1]][3], ",")[[1]][3])
      
      if((strsplit(file_content[i], "[ \t]")[[1]][5] == "1") || (strsplit(file_content[i], "[ \t]")[[1]][5] == "lefteye")){
        
        landmarks[["lefteye"]] <- matrix(c(x, y, z), nrow=1, dimnames=list(c(), c("x", "y", "z")))
      }
      else if((strsplit(file_content[i], "[ \t]")[[1]][5] == "3") || (strsplit(file_content[i], "[ \t]")[[1]][5] == "righteye")){
        
        landmarks[["righteye"]] <- matrix(c(x, y, z), nrow=1, dimnames=list(c(), c("x", "y", "z")))
      }
      else if((strsplit(file_content[i], "[ \t]")[[1]][5] == "2") || (strsplit(file_content[i], "[ \t]")[[1]][5] == "nasion")){
        
        landmarks[["nasion"]] <- matrix(c(x, y, z), nrow=1, dimnames=list(c(), c("x", "y", "z")))
      }
      else if((strsplit(file_content[i], "[ \t]")[[1]][5] == "5") || (strsplit(file_content[i], "[ \t]")[[1]][5] == "nose")){
        
        landmarks[["nose"]] <- matrix(c(x, y, z), nrow=1, dimnames=list(c(), c("x", "y", "z")))
      }
      else if((strsplit(file_content[i], "[ \t]")[[1]][5] == "6") || (strsplit(file_content[i], "[ \t]")[[1]][5] == "leftalare")){
  
        landmarks[["leftalare"]] <- matrix(c(x, y, z), nrow=1, dimnames=list(c(), c("x", "y", "z")))
      }
      else if((strsplit(file_content[i], "[ \t]")[[1]][5] == "7") || (strsplit(file_content[i], "[ \t]")[[1]][5] == "rightalare")){
  
        landmarks[["rightalare"]] <- matrix(c(x, y, z), nrow=1, dimnames=list(c(), c("x", "y", "z")))
      }
      else if((strsplit(file_content[i], "[ \t]")[[1]][5] == "8") || (strsplit(file_content[i], "[ \t]")[[1]][5] == "subnasale")){
  
        landmarks[["subnasale"]] <- matrix(c(x, y, z), nrow=1, dimnames=list(c(), c("x", "y", "z")))
      }
      else if((strsplit(file_content[i], "[ \t]")[[1]][5] == "9") || (strsplit(file_content[i], "[ \t]")[[1]][5] == "leftmouth")){
  
        landmarks[["leftmouth"]] <- matrix(c(x, y, z), nrow=1, dimnames=list(c(), c("x", "y", "z")))
      }
      else if((strsplit(file_content[i], "[ \t]")[[1]][5] == "10") || (strsplit(file_content[i], "[ \t]")[[1]][5] == "rightmouth")){
  
        landmarks[["rightmouth"]] <- matrix(c(x, y, z), nrow=1, dimnames=list(c(), c("x", "y", "z")))
      }
    }
  }
  
  (landmarks)
}

readImageData <- function(file){
  
  lines <- strsplit(readLines(file), "[ ]")
  
  row <- length(lines)
  col <- length(lines[[1]])
  
  img <- matrix(rep(0, row*col), ncol=col)
  
  for(i in 1:row){
    
    img[i,] <- as.numeric(lines[[i]])
  }
  
  (img)
}

#Reads a 3D matrix from file.
#The file must have 3 columns (X, Y, Z, respectvelly)
#input:
#   file = the path to the matrix file
#output:
#   a n x 3 matrix, where n is the number of rows (records)
data2matrix <- function(file){
  
  text <- readLines(file)
  
  Xs <- as.numeric(strsplit(text[1], "[ \t]")[[1]])
  Ys <- as.numeric(strsplit(text[2], "[ \t]")[[1]])
  Zs <- as.numeric(strsplit(text[3], "[ \t]")[[1]])
  
  data <- matrix(c(Xs, Ys, Zs), ncol=3, dimnames=list(c(), c("x", "y", "z")))
  (data)
}

# Converts an image represented by a NxM matrix into a N*Mx3 matrix----
# ...
img2matrix <- function(file, samplingVariation=1){
  
  img <- file
  if(is.character(file))
    img <- readImageData(file)
  
  #img <- t(img)
  
  N <- length(img[1,])
  M <- length(img[,1])
  
  m <- matrix(rep(0, 3), ncol=3)
  
  k <- 0
  total <- length(img[which(img != 0)])
  i <- 1
  while(i <= M){
    j <- 1
    while(j <= N){
      
      if(img[i, j] != 0){
        
        m <- rbind(m, c(j, i, img[i,j]))
        k <- k + 1
        
        cat(k*100/(total), "%\n")
      }
      j <- j + samplingVariation
    }
    i <- i + samplingVariation
  }
  
  maxY <- max(m[,2])
  
  m[,2] <- maxY + 1 - m[,2]
  
  (m[-1,])
}

#Reads a ABS file and generates a matrix, where each cell contains the range information.----
# 'file' is a string with the path to the ABS file.
# A 2D matrix is returned.
abs2matrix <- function(file){
  #reads the content of the ABS file
  file_content <- readLines(file)
  
  #extracts the number of rows
  row <- as.numeric(strsplit(file_content[1], " rows")[[1]])
  
  #extracts the number of columns
  col <- as.numeric(strsplit(file_content[2], " columns")[[1]])
  
  #extracts the flags
  #flags <- as.numeric(strsplit(file_content[4], "[ \t]")[[1]])
  
  #extracts the Xs values
  Xs <- as.numeric(strsplit(file_content[5], "[ \t]")[[1]])
  
  #extracts the Ys values
  Ys <- as.numeric(strsplit(file_content[6], "[ \t]")[[1]])
  
  #extracts the Zs values
  Zs <- as.numeric(strsplit(file_content[7], "[ \t]")[[1]])
  
  #removes the invalid readings
  Xs <- Xs[-which(Xs == -999999)]
  Ys <- Ys[-which(Ys == -999999)]
  Zs <- Zs[-which(Zs == -999999)]
  
  m <- matrix(c(Xs, Ys, Zs), ncol=3, dimnames = list(c(), c("x", "y", "z")))
  (m)
}

#converts the matrix of 3d points into a range image.
#inputs:
#   m = the 3d point matrix
#   file = the path to the image file with no extensions. e.g. "/home/username/face0001"
#   xfactor = the proportional factor to the x axis on the range image
#   yfactor = the proportional factor to the y axis on the range image
#   col = the number of columns (horizontal pixels) of the image
#   row = the number of rows (vertical pixels) of the image
#   ldmk = the list with the landmark coordinates to be updated into the new system
#outputs:
#   a list with "img", that is the 2D matrix (image). The image file will be created or replaced in JPG format
#   and with "ldmk", that is a list with the updated landmark coordinates to the image coordinate system
matrix2image <- function(m, file, dx = 0, dy = 0, col=640, row=480, ldmk=0){
  
  #finds the min and max values of both axis (X and Y)
  minXd <- min(m[,1])
  maxXd <- max(m[,1])
  minYd <- min(m[,2])
  maxYd <- max(m[,2])
  minZd <- min(m[,3])
  maxZd <- max(m[,3])
  
  #finds the variation in both axis (X and Y), delta_X and delta_Y
  if(dx == 0)
    dx <- (maxXd - minXd)/(col - 1)
  if(dy == 0)
    dy <- (maxYd - minYd)/(row - 1)
  
#   meanX <- (maxXd - minXd)/2 + minXd
#   tx <- trunc((meanX - minXd)/dx + 1.0 - (col/2))
#   
#   meanY <- (maxYd - minYd)/2 + minYd
#   ty <- trunc(-(meanY - maxYd)/dy + 1.0 - (row/2))
tx <- 0
ty <- 0
  
  if(length(ldmk) > 1){
    
    tx <- trunc((ldmk[["nose"]][1] - minXd)/dx + 1.0) - trunc(col/2)
    ty <- trunc(-(ldmk[["nose"]][2] - maxYd)/dy + 1.0) - trunc(row/2)
  }
  
  #instantiates the img matrix and fills it with 0s
  img <- matrix(rep(0, col * row), ncol=row)
  
  for(i in 1:length(m[,1])){
    
    x <- trunc((m[i,1] - minXd)/dx + 1.0) - tx
    y <- trunc(-(m[i,2] - maxYd)/dy + 1.0) - ty
    
    z <- m[i,3] - minZd + 1
    
    if(x > 0 && x <= col && y > 0 && y <= row && z > img[x,y]){
      
      img[x,y] <- z
    }
    
    #cat(i/length(m[,1]) * 100, "%\n")
  }
  
  #if the landmarks are passed as input
  if(length(ldmk) > 1){
    
    #for each ldmk, 
    for(i in 1:length(ldmk)){
      
      x <- trunc((ldmk[[i]][1] - minXd)/dx + 1.0)
      y <- trunc(-(ldmk[[i]][2] - maxYd)/dy + 1.0)
      z <- m[i,3] - minZd + 1
      
      ldmk[[i]][1] <- x
      ldmk[[i]][2] <- y
      ldmk[[i]][3] <- z
    }
  }
  
  if(file != ""){
    write(img, paste(file, ".jpg.dat", sep=""), ncolumns=col)
    img <- t(img)
    my.writeJPEG(img, paste(file, ".jpg", sep=""), quality=1)
  }
  else
    img <- t(img)
  
  (list(img = img, ldmk = ldmk))
}

#Computes a sphere around the given points.
# The points input must be a matrix (col = coordinates, row = points).
# It returns a list with a central point and the radius.
getSphere <- function(points){
  
  col = length(points[1,])
  row = length(points[,1])
  
  meanPoint <- rep(0, col)
  
  for(i in 1:col){
    
    meanPoint[i] = mean(points[,i])
  }
  
  distances <- rep(0, row)
  
  for(i in 1:row){
    
    m <- rbind(meanPoint, points[i,])
    #computes the distance between the ith point to the mean point
    distances[i] <- dist(m)[1]
    #prints the progress
    #cat(i/row*100, "%\n")
  }
  
  radius <- max(distances)
  
  (list(radius = radius, meanPoint = meanPoint))
}

#Computes a set of circles around the given sphere.
# It needs the list returned by getSphere ('sphere');
# and the number of (horizontal) circles.
# It returns a matrix with the points in each circle gernerated. The mean point is the first row of the matrix.
# The number of points in each circle is equal the double of the number of circles.
getSamplePoints <- function(sphere, numberOfCircles){
  
  points <- matrix(sphere$meanPoint, nrow=1)
  
  step <- 2 * sphere$radius / (numberOfCircles + 1)
  numberOfPoints <- 2 * numberOfCircles
  
  for(i in 1:numberOfCircles){
    
    #computes the new radius
    r <- sqrt(sphere$radius^2 - (sphere$radius - (i  * step))^2)
    
    #computes the new Y
    y <- sphere$meanPoint[2] + (sphere$radius - (i * step))
    
    #computes the X-Z-plane (horizontal) circle
    circle <- getCirclePoints(sphere$meanPoint[1], sphere$meanPoint[3], r, numberOfPoints)
    #plot(circle)
    
    #adds the circle's points into the matrix
    for(j in 1:numberOfPoints){
      
      if(circle[j,2] < sphere$meanPoint[3]){
        points = rbind(points, matrix(c(circle[j,1], y, circle[j,2]), nrow=1))
      }
    }
  }
  
  (points)
}

#Generates 'n' points of a circle with it's radius = 'r' and center = 'x', 'y'.
# It returns a matrix with the generated points.
getCirclePoints <- function(x, y, r, n){
  
  points <- matrix(rep(0, 2), nrow=1)
  
  for(i in 1:n){
    
    points <- rbind(points, matrix(c((r * cos(i * 2 * pi/n)) + x, (r * sin(i * 2 * pi/n)) + y), nrow=1))
  }
  
  (points[-1,])
}

measureEnergy <- function(curve){
  
  n <- length(curve)
  
  energy <- 0
  
  for(i in 2:(n-1)){
    #computes the horizontal angle of the point with its predecessor
    a <- abs(atan((3*curve[i] - 3*curve[i+1])/(1)))
    #computes the horizontal angle of the point with its successor
    b <- abs(atan((3*curve[i] - 3*curve[i-1])/(-1)))
    #computes the energy
    energy <- energy + (pi/(pi - (a + b)) - 1)^2
  }
  
  (list(energy = energy, meanEnergy = (energy/(n-2))))
}

# x = matrix(c(3,5,6,1,4,8), nrow=2, byrow=TRUE)
# circles = getSamplePoints(getSphere(x), 50)

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
    #curve <- gaussianSmooth(curve, c(smooth))
    curve <- gaussianSmooth(c(rep(curve[1], 2*smooth), curve, rep(curve[n], 2*smooth)), c(smooth))[(1 + 2*smooth):(n + 2*smooth)]
  
  (curve)
}