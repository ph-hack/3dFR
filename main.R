#loads the needed packages
library(jpeg) #open jpg images
source("pre_process.R")
source("interpolation.R")
source("classf.R")
source("timing.R")
source("teste.R")

# scatterplot3d(circles, highlight.3d=TRUE, col.axis="blue", col.grid="lightblue", main="Sphere", pch=20, box=FALSE)

# loads all ABS files from a directory into matrices.
# input:
#   directory = the directory where the ABS files are;
# output:
#   a list containing 'faces', a list with the matrices and
#   'names', a vector containing the file names
loadFaces <- function(directory){
  
  names <- dir(directory)
  n <- length(names)
  
  faces <- list()
  
  for(i in 1:n){
    
    faces[[i]] <- abs2matrix(file= paste(directory, names[i], sep=""))
    
    cat(i/n * 100, "%\n")
  }
  
  (list(faces = faces, names = names))
}

pre_process <- function(faceDir, ldmkDir, toDir, transfFile, range=0, exception="",keepProportions=TRUE, col=150, row=210, logFile="", append=FALSE){
  
  begin <- getTime()
  start <- getTime()
  
  #gets all file names from the face directory
  names <- dir(faceDir)
  
  #gets only the specified range
  if(length(range) > 1)
    names <- names[range[1]:range[2]]
  
  #removes the exceptions
  if(exception[1] != "")
    for(i in 1:length(exception)){
      exc <- which(names == exception[i])
      if(length(exc) > 0)
        names <- names[-exc]
    }
  
  ldmkFiles <- dir(ldmkDir)
  
  #reads the transformation data
  transf_data <- getTransfData(transfFile)
  
  n <- length(names)
  
  cat("n = ", n, "\n")
  
  dx <- 0
  dy <- 0
  
  cat("Reading file names and transformation data: ", crono.end(start), "\n", file=logFile, append=append)
  
  for(i in 1:n){
    
    start <- getTime()
    
    cat("processing ", names[i], "-------------------------------------------------\n")
    cat("Face ", names[i], "\n", file=logFile, append=TRUE)
    
#     face <- abs2matrix(concatenate(c(faceDir, names[i])))
#     cat("3D points loaded\n")
    
    ldmk <- readLandmark(file=concatenate(c(ldmkDir, ldmkFiles[which(ldmkFiles == paste(getFaceID(names[i]), ".groundtruth.ldmk", sep=""))])))
    cat("Landmarks loaded\n")
    
#     face <- cutTheFaceByNoseLength(face, ldmk)
#     cat("Face cut\n")
#     
#     face <- alignFace(face, ldmk)
#     cat("Face rotated\n")
    
    face <- registerFaces(concatenate(c(faceDir, names[i])), transf_data, ldmk)
    
    cat("Registering a face: ", crono.end(start), "\n", file=logFile, append=TRUE)
    
    start <- getTime()
    
    if(keepProportions && dx == 0 && dy == 0){
      
      minXd <- min(face$points[,1])
      maxXd <- max(face$points[,1])
      minYd <- min(face$points[,2])
      maxYd <- max(face$points[,2])
      
      dx <- (maxXd - minXd)/(col - 1)
      dy <- (maxYd - minYd)/(row - 1)
    }
    
    #img <- matrix2image(face$points, concatenate(c(toDir, getFaceID(names[i]), ".original")), dx, dy, col, row, face$ldmk)
    img <- matrix2image(face$points, "", dx, dy, col, row, face$ldmk)
    cat("Face resampled as range image\n")
    
    writeLandmark(img$ldmk, concatenate(c(toDir, getFaceID(names[i]), ".ldmk")))
    cat("New landmarks saved\n")
    
    cat("Creating range image: ", crono.end(start), "\n", file=logFile, append=TRUE)
    
    start <- getTime()
    
    #img <- connectedComponentAnalysis(img$img, concatenate(c(toDir, getFaceID(names[i]), ".connected.jpg")))
    img <- connectedComponentAnalysis(img$img)
    cat("Connected Component Analysis performed\n")
    
    img <- holeCorrection(img, concatenate(c(toDir, getFaceID(names[i]), ".unholed.jpg")), 20)
    cat("Hole correction made\n")
    
    cat("Hole correction: ", crono.end(start), "\n", file=logFile, append=TRUE)
    
    start <- getTime()
    
    #curves <- curvatureImage(img, factor=3)
    #cat("Curvature for outlier correction computed\n")
    
    #img <- outlierCorrection2(img, curves, concatenate(c(toDir, getFaceID(names[i]), ".unoutlied.jpg")), progress=FALSE)
    #cat("Outlier correction made\n")
    
    #weights <- weightedRegions(img, readLandmark(concatenate(c(toDir, getFaceID(names[i]), ".ldmk"))), concatenate(c(toDir, getFaceID(names[i]), ".weights.jpg")))
    #cat("Weights computed\n")
    
    #curves <- curvatureImage(img, weights, concatenate(c(toDir, getFaceID(names[i]), ".curvesweighted.jpg")), 3)
    curves <- curvatureImage(img, "", concatenate(c(toDir, getFaceID(names[i]), ".curves.jpg")), 3)
    cat("Curvature measures made\n")
    
    cat("Computing curvature image: ", crono.end(start), "\n", file=logFile, append=TRUE)
    
    cat(i/n * 100, "%\n")
  }
  
  cat("Total: ", crono.end(begin), "\n", file=logFile, append=TRUE)
}

pre_processJustWeigths <- function(names, ldmkDir, weightDir, fromDir, toDir, factor=1){
  
  if(is.character(names))
    names <- dir(names)
  
  ldmkFiles <- dir(ldmkDir)
  
  n <- length(names)
  
  cat("n = ", n, "\n")
  
  for(i in 1:n){
    
    cat("processing ", names[i], "-------------------------------------------------\n")
    
    img <- readJPEG(paste(fromDir, names[i], sep=""))
    
    #weights <- weightedRegions(img, readLandmark(paste(ldmkDir, paste(getFaceID(names[i]), ".ldmk", sep=""), sep="")), paste(toDir, paste(getFaceID(names[i]), paste(".weights", ".jpg", sep=""), sep=""), sep=""))
    #cat("Weights computed\n")
    
    weight <- readJPEG(paste(weightDir, paste(getFaceID(names[i]), ".weights.jpg", sep=""), sep=""))
    
    curves <- curvatureImage(img, weight, paste(toDir, paste(getFaceID(names[i]), paste(".curves", ".jpg", sep=""), sep=""), sep=""), factor)
    cat("Curvature measures made\n")
    
    cat(i/n * 100, "%\n")
  }
}

getWeightedRangeImages <- function(fromDir, toDir, factor=1){
  
  names <- dir(fromDir)
  
  n <- length(names)
  
  for(i in 1:n){
    
    names[i] <- paste(strsplit(names[i], "[.]")[[1]][1], ".jpg", sep="")
      
    ponderateImage(paste("processed_data/", names[i], sep=""), 
                   paste("processed_data/", paste(getFaceID(names[i]), paste(".weights", ".jpg", sep=""), sep=""), sep=""), 
                   factor, paste(toDir, paste(getFaceID(names[i]), ".weighted.jpg", sep=""), sep=""))
    
    cat(i * 100 / n, "%\n")
  }
}

getRangeImages <- function(fromDir, toDir){
  
  names <- dir(fromDir)
  
  n <- length(names)
  
  for(i in 1:n){
    
    names[i] <- paste(strsplit(names[i], "[.]")[[1]][1], ".jpg", sep="")
    writeJPEG(readJPEG(concatenate(c("processed_data/", names[i]))), paste(toDir, names[i], sep=""))
    
    cat(i * 100/n, "%\n")
  }
}

getLinesFiles <- function(imgsDir, ldmkDir, toDir, smoothness=0, halfFace=FALSE){
  
  if(length(imgsDir) == 1){
    names <- dir(imgsDir)
    names <- names[-which(my.strsplit(names, "[.]", 0) == "jpg")]
  }
  else{
    names <- imgsDir
    imgsDir <- ""
  }
  
  n <- length(names)
  
  for(i in 1:n){
    
    cat(names[i], ": ")
    
    lines <- getMainLines(readImageData(concatenate(c(imgsDir, names[i]))),
                          readLandmark(concatenate(c(ldmkDir, getFaceID(names[i]), ".ldmk"))),
                          concatenate(c(toDir, getFaceID(names[i]), ".lines")), halfFace, smoothness)
    
    cat(i * 100/n, "%\n")
  }
}

getOutliedImgs <- function(imgDir, toDir, meanImg, range=0, threshold=0.5, progress=TRUE){
  
  imgs <- dir(imgDir)
  imgs <- imgs[-which(my.strsplit(imgs, "[.]", 0) == "jpg")]
  if(length(range) == 2)
    imgs <- imgs[range[1]:range[2]]
  
  m <- length(imgs)
  
  for(i in 1:m){
    
    img <- outlierCorrection4(concatenate(c(imgDir, imgs[i])), meanImg, threshold, concatenate(c(toDir, getFaceID(imgs[i]), ".unoutlied")))
    
    if(progress)
      cat(i * 100/m, "%\n")
  }
}

separeteDataSet <- function(fromDir, trainingReference, testReference, trainingDir, testDir, sulfix, progress=TRUE){
  
  training <- dir(trainingReference)
  N <- length(training)
  test <- dir(testReference)
  M <- length(test)
  
  for(i in 1:N){
    
    writeJPEG(readJPEG(concatenate(c(fromDir, getFaceID(training[i]), sulfix, ".jpg"))),
             concatenate(c(trainingDir, getFaceID(training[i]), sulfix, ".jpg")), quality=1)
    
    write(t(readImageData(concatenate(c(fromDir, getFaceID(training[i]), sulfix, ".jpg.dat")))),
             concatenate(c(trainingDir, getFaceID(training[i]), sulfix, ".jpg.dat")), ncolumns=150)
    
    if(progress)
      cat("training ", i*100/(N+M), "%\n")
  }
  
  for(i in 1:M){
    
    writeJPEG(readJPEG(concatenate(c(fromDir, getFaceID(test[i]), sulfix, ".jpg"))),
             concatenate(c(testDir, getFaceID(test[i]), sulfix, ".jpg")), quality=1)
    
    write(t(readImageData(concatenate(c(fromDir, getFaceID(test[i]), sulfix, ".jpg.dat")))),
             concatenate(c(testDir, getFaceID(test[i]), sulfix, ".jpg.dat")), ncolumns=150)
    
    if(progress)
      cat("test ", (i+N)*100/(N+M), "%\n")
  }
}

resizeImgs <- function(fromDir, toDir, col=25, row=35, reason=1, progress=TRUE){
  
  imgs <- dir(fromDir)
  imgs <- imgs[-which(my.strsplit(imgs, "[.]", 0) == "jpg")]
  
  n <- length(imgs)
  
  for(i in 1:n){
    
    img <- readImageData(concatenate(c(fromDir, imgs[i])))
    img <- resize(img, concatenate(c(toDir, getFaceID(imgs[i]), ".jpg")), col, row, reason)
    
    if(progress)
      cat(i*100/n, "%\n")
  }
}

getMeanFacesByPerson <- function(fromDir, toDir, progress=TRUE){

  training <- fromDir
  if(length(fromDir) == 1){
    
    training <- dir(fromDir)
    training <- training[-which(my.strsplit(training, "[.]", 0) == "jpg")]
  }
  else
    fromDir <- ""

  trainingClasses <- getPersonID(training)
  classes <- unique(trainingClasses)

  nC <- length(classes)
  i <- 1

  for(cl in classes){

    samples <- which(trainingClasses == cl)

    if(length(samples) == 1)
      cat("fuck!\n")
    
    face <- meanFace(concatenate(list(fromDir, training[samples])), concatenate(c(toDir, cl, "d000.mean.jpg")), FALSE, FALSE)

    if(progress){
      cat(i * 100/nC, "\n")
      i <- i + 1
    }
  }
}
