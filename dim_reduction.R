# generates and saves the eigenVectors in TXT ##########################################################################
# it needs the folder where the images are, if the R isn't working in the images' folder ('directory')
#, the folder where the eigenVectors will be saved in ('eigenDir')
#, the folder where the means will be saved in ('meanDir')
# and the amount of energy which must be kept with the dimensionality reduction
eigenVectors <- function(directory="", eigenDir="", meanDir="", energy=0.99){
  imgs <- dir(directory)
  imgs <- imgs[-which(my.strsplit(imgs, "[.]", 0) == "jpg")]
  m <- length(imgs)
  
  #reads the first training sample
  trainingImg <- t(readImageData(paste(directory, imgs[1], sep="")))
  width <- length(trainingImg[,1])
  height <- length(trainingImg[1,])
  dim(trainingImg)  <- width * height
  
  trainingMatrix <- matrix(trainingImg, nrow=1)
  for(j in 2:m){
    
    #reads the jth training sample, for all j > 1
    trainingImg <- t(readImageData(paste(directory, imgs[j], sep="")))
    dim(trainingImg)  <- width * height
    
    #appends it to the matrix
    trainingMatrix <- rbind(trainingMatrix, matrix(trainingImg, nrow=1))
    
    cat("building training matrix: ", (j*100/m), "%\n")
  }
  
  print("applying PCA...")
  pcaResult <- pca(trainingMatrix, energy)
  
  #saves the best eigenVectors as images
  for(i in 1:pcaResult$nDim){
    
    eigenVector <- pcaResult$eigenVectors[,i] #gets the ith eigenVector
    
    if(eigenDir == "") eigenDir = directory
    write(eigenVector, concatenate(c(eigenDir, "eigenFace", as.character(i), ".txt")), sep=" ")
    
    cat("saving eigenFaces: ", (i*100/pcaResult$nDim), "%\n")
  }
  
  #saves the mean
  write(pcaResult$u, paste(meanDir, "means.txt", sep=""), sep=" ")
  
  (pcaResult$eigenVectors)
}

# applies the PCA into the data matrix ##############################################################################
# it needs the data matrix where to apply the PCA ('trainingMatrix')
# and the amount of energy which must be kept
# it returns a list with the number of dimensions
# and the eigenVectors
pca <- function(trainingMatrix, energy){
  cMean <- colMeans(trainingMatrix) #finds the mean of all cols (features)
  print("col means calculated...")
  
  #normalizes the transposed matrix so the mean of the cols = 0
  trainingMatrix <- trainingMatrix - matrix(rep(cMean, dim(trainingMatrix)[1]), ncol=dim(trainingMatrix)[2], byrow=TRUE)
  print("matrix normalized...")
  
  transposedTrainingMatrix <- t(trainingMatrix)
  print("transposed calculated...")
  
  covarianceMatrix <- transposedTrainingMatrix %*% trainingMatrix #computes the covariance matrix
  print("covariance matrix calculated...")
  
  E <- eigen(covarianceMatrix) #decomposes the covariance matrix into UtDU matrices
  print("decomposition made...")
  
  #finds the number of dimensions with 99% of energy
  nDimensions <- findMostSignificantDimensions(E$values, energy)
  
  (list(nDim = nDimensions,
        eigenVectors = E$vectors,
        u = cMean))
}

# finds the dimensions that carries 'energy'(a percentage) part of ####################################################
# the total energy.
# it needs the vector with the values ('values')
# and the minumum percentage of energy the dimensions must carry ('energy')
# by default the 'energy' is 0.99 (99% of energy)
# it returns the number of dimensions that carry 'energy' part of all energy.
findMostSignificantDimensions <- function(values, energy=0.99){
  i = 1
  totalValues <- sum(values)
  
  while(sum(values[1:i])/totalValues < energy){
    i = i + 1
  }
  
  (i)
}

# maps and save the images with the eigenVectors ########################################################################
# it needs the folder of the images ('imgsDir')
#, the folder containing the 'means.txt' file with the dimensions means
#, the folder where the mapped images will be saved in ('mapDir')
#, the folder of the eigenVectors ('eigenDir')
# and the dimensions to be ignored ('ignDim'), which can be one number or a vector with two numbers
# when 'ignDim' is one number, the 'ignDim' 1st eigenvectors will be ignored
# when 'ignDim' is two numbers, the interval between 'ignDim'[1] and 'ignDim'[2] wont be ignored
# By default all folders are set to the current R directory
mapWithEigenVectors <- function(imgsDir="", meanDir="", mapDir="", eigenDir="", ignDim=0){
  
  imgs <- dir(imgsDir) #gets the images' filenames
  imgs <- imgs[-which(my.strsplit(imgs, "[.]", 0) == "jpg")]
  m <- length(imgs) #gets the number of images
  
  eigens <- dir(eigenDir) #gets the eigenVectors' filenames]
  if(length(ignDim) == 1){
    eigens <- eigens[which(getPersonId(eigens) > ignDim)]
  }
  else{
    eigens <- eigens[which(getPersonId(eigens) >= ignDim[1])]
    eigens <- eigens[which(getPersonId(eigens) <= ignDim[2])]
  }
  d <- length(eigens) #gets the number of dimensions
  
  #reads the features' means
  u <- as.numeric(scan(paste(meanDir, "means.txt", sep=""), quiet=TRUE))
  
  #reads all the eigenVectors
  eigenVectors <- matrix(as.numeric(scan(paste(eigenDir, eigens[1], sep=""), quiet=TRUE)), nrow=1)
  for(j in 2:d){
    #reads the jth eigenFace
    eigenVectors <- rbind(eigenVectors, matrix(as.numeric(scan(paste(eigenDir, eigens[j], sep=""), quiet=TRUE)), nrow=1))
  }
  print("EigenVectors read...\n")
  
  for(i in 1:m){ #for each image
    
    #reads the jth training image
    img <- t(readImageData(paste(imgsDir, imgs[i], sep="")))
    dim(img) <- c(length(img[,1]) * length(img[1,]), 1) #converts it to a matrix with 1 column
    img <- img - matrix(u, ncol=1) #normalize by the mean
    #cat("training image ", j, " read...\n")
    
    #maps the image to new space
    newImg <- rep(0, d)
    
    for(k in 1:d){
      #reads the jth eigenFace
      eigenFace <- eigenVectors[k,]
      dim(eigenFace) <- c(1, length(eigenFace)) #converts it to a matrix with 1 row
      
      newImg[k] <- (eigenFace %*% img)[1,1] #maps
    }
    #saves the new image
    write(newImg, paste(mapDir, paste("mapped", paste(as.character(strsplit(imgs[i], "[.]")[[1]][1]), ".txt", sep=""), sep=""), sep=""), sep=" ")
    
    print(i*100/m) #prints the progress
  }
}

# converts an eigenVector into image ####################################################################################
# it needs the eigenvector ('eigenVector')
#, the image dimensions ('w' and 'h')
#, the name of the image to be saved ('filename')
# and the directory where to save the eigenface ('directory')
eigenFace <- function(eigenVector, w=0, h=0, filename="eigenFace.jpg", directory=""){
  
  maxV <- max(eigenVector)
  minV <- min(eigenVector)
  
  eigenVector <- (eigenVector - minV)/(maxV - minV)
  if(w != 0 && h != 0) dim(eigenVector) <- c(w, h)
  
  writeImage(eigenVector, paste(directory, filename, sep=""), type="jpg")
}

# generates and saves the eigenVectors from LDA in TXT ##########################################################################
# it needs the folder where the images are, if the R isn't working in the images' folder ('imgsDir')
#, the folder where the eigenVectors will be saved in ('eigenDir')
# and the method to be used, wheather the MASS algorithm (default) or not
ldaVectors <- function(imgsDir="", eigenDir="", method="R"){
  imgs <- dir(imgsDir)
  #imgs <- imgs[-which(my.strsplit(imgs, "[.]", 0) == "jpg")]
  m <- length(imgs)
  
  trainingImg <- 0
  width <- 0
  height <- 0
  #reads the first training sample
  if(strsplit(imgs[1], "[.]")[[1]][2] == "jpg"){
    trainingImg <- t(readJPEG(paste(imgsDir, imgs[1], sep="")))
    width <- length(trainingImg[,1])
    height <- length(trainingImg[1,])
    dim(trainingImg)  <- width * height
  }
  else{
    trainingImg <- scan(paste(imgsDir, imgs[1], sep=""), quiet=TRUE)
  }
  
  trainingMatrix <- matrix(trainingImg, nrow=1)
  for(j in 2:m){
    
    #reads the jth training sample, for all j > 1
    if(strsplit(imgs[j], "[.]")[[1]][2] == "jpg"){
      trainingImg <- t(readJPEG(paste(imgsDir, imgs[j], sep="")))
      dim(trainingImg)  <- width * height
    }
    else{
      trainingImg <- scan(paste(imgsDir, imgs[j], sep=""), quiet=TRUE)
    }
    
    #appends it to the matrix
    trainingMatrix <- rbind(trainingMatrix, matrix(trainingImg, nrow=1))
    
    cat("building training matrix: ", (j*100/m), "%\n")
  }
  
  print("applying LDA...")
  if(method == "R"){
    ldaResult <- lda(trainingMatrix, factor(as.character(getPersonID(imgs)))) 
  }
  else{
    ldaResult <- myLda(trainingMatrix, imgsDir)
  }
  
  d <- 0
  if(strsplit(imgs[j], "[.]")[[1]][2] == "jpg"){
    if(method == "R"){
      d <- length(ldaResult$lev)-1
    }
    else{
      d <- ldaResult$mc - 1
    }
  }
  else{
    d <- length(trainingMatrix[1,])
  }
  
  #saves the best C - 1 ldaVectors
  for(i in 1:d){  
    ldaVector <- 0
    if(method == "R"){
      ldaVector <- ldaResult$scaling[,i] #gets the ith ldaVector
    }
    else{
      ldaVector <- ldaResult$vectors[,i]
    }
    
    if(eigenDir == "") eigenDir = imgsDir
    if(is.complex(ldaVector)){
      write(Re(ldaVector), paste(eigenDir, paste("ldaVector", paste(as.character(i), ".txt", sep=""), sep=""), sep=""), sep=" ")
      cat("saving complex ldaVector: ", (i*100/d), "%\n")
    }
    else{
      write(ldaVector, paste(eigenDir, paste("ldaVector", paste(as.character(i), ".txt", sep=""), sep=""), sep=""), sep=" ")
      cat("saving ldaVector: ", (i*100/d), "%\n")
    }
  }
}

#maps and save the images with the eigenVectors ########################################################################
#it needs the folder of the images ('imgsDir')
#, the folder where the mapped images will be saved in ('mapDir')
#, the folder of the eigenVectors ('eigenDir')
#By default all folders are set to the current R directory
mapWithLDAVectors <- function(imgsDir="", mapDir="", eigenDir=""){
  
  imgs <- dir(imgsDir) #gets the images' filenames
  m <- length(imgs) #gets the number of images
  
  eigens <- dir(eigenDir) #gets the eigenVectors' filenames
  d <- length(eigens) #gets the number of dimensions
  
  #reads the features' means
  #u <- as.numeric(scan(paste(meanDir, "means.txt", sep=""), quiet=TRUE))
  
  #reads all the eigenVectors
  eigenVectors <- matrix(scan(paste(eigenDir, eigens[1], sep=""), quiet=TRUE), nrow=1)
  for(j in 2:d){
    #reads the jth eigenFace
    eigenVectors <- rbind(eigenVectors, matrix(as.numeric(scan(paste(eigenDir, eigens[j], sep=""), what=complex(), quiet=TRUE)), nrow=1))
    cat("reading lda vectors: ", (j*100/d), "%\n")
  }
  print("EigenVectors read...\n")
  
  for(i in 1:m){ #for each image
    
    #reads the jth training image
    img <- 0
    if(strsplit(imgs[j], "[.]")[[1]][2] == "jpg"){
      img <- t(readJPEG(paste(imgsDir, imgs[j], sep="")))
      dim(img) <- c(length(img[,1]) * length(img[1,]), 1) #converts it to a matrix with 1 column
    }
    else{
      img <- scan(paste(imgsDir, imgs[j], sep=""), quiet=TRUE)
      dim(img) <- c(length(img), 1) #converts it to a matrix with 1 column
    }
    
    #maps the image to new space
    newImg <- rep(0, d)
    
    for(k in 1:d){
      #reads the jth eigenFace
      eigenFace <- eigenVectors[k,]
      dim(eigenFace) <- c(1, length(eigenFace)) #converts it to a matrix with 1 row
      
      newImg[k] <- (eigenFace %*% img)[1,1] #maps
    }
    #saves the new image
    write(newImg, paste(mapDir, paste("mapped", paste(as.character(strsplit(imgs[i], "[.]")[[1]][1]), ".txt", sep=""), sep=""), sep=""), sep=" ")
    
    print(i*100/m) #prints the progress
  }
}