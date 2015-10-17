concatenate <- function(words, sep=""){
  
  phrase <- ""
  n <- length(words)
  
  if(is.list(words)){
    
    for(i in 1:n)
      phrase <- paste(phrase, words[[i]], sep=sep)
  }
  else{
    
    for(i in 1:n)
      phrase <- paste(phrase, words[i], sep=sep)
  }
  
  (phrase)
}

discretize <- function(x, n, range=c(0,1)){
  
  minX <- min(x)
  maxX <- max(x)
  
  deltaX <- (range[2] - range[1])/n
  
  x <- ((x - minX)*(range[2]-range[1]))/(maxX - minX)
  
  x <- floor(x/deltaX) * deltaX
  
  (x)
}

is.file <- function(string){
  
  chars <- strsplit(string, "")[[1]]
  
  if(chars[length(chars)] == "/")
    (FALSE)
  else
    (TRUE)
}

list2vector <- function(list){
  
  n <- length(list)
  
  v <- rep(0, n)
  
  for(i in 1:n){
    
    v[i] <- list[[i]]
  }
  
  (v)
}

my.maxORmin <- function(data1, data2, maxORmin="max"){
  
  col <- length(data1[1,])
  row <- length(data1[,1])
  
  output <- matrix(rep(0, col*row), ncol=col)
  
  for(i in 1:row){
    for(j in 1:col){
      
      if(maxORmin == "max")
        output[i,j] <- max(data1[i,j], data2[i,j])
      else
        output[i,j] <- min(data1[i,j], data2[i,j])
    }
  }
  
  (output)
}

AND <- function(x){
  
  return (Reduce(function(p, x){
    return(p && x)
  }, x, TRUE))
}

# Computes the cosine distance between two vectors ----
cosineDist <- function(v1, v2){
  
  return(1 - sum(v1 * v2)/(sqrt(sum(v1^2)) * sqrt(sum(v2^2))))
}

euclidianDist <- function(v1, v2){
  
  return(sqrt(sum((v1-v2)^2)))
}

# Retrieves all values of a given field of a list in a given level ----
# input:
#   theList = the list which have the values;
#   index = the field
#   level = the level of 'theList' where is the field 'index'
# output:
#   a list containing all occurences of such field
# Example:
#   y = list(list(p=5, d=10), list(p=3, d=10), list(p=2, d=5))
#   the output of the call with theList=y, index="d", level=2)
#   will be 10,10,5
getAllFieldFromList <- function(theList, index=1, level=1){
  
  if(level > 1){
    
    n <- length(theList)
    field <- list()
    
    for(i in 1:n){
      
      field[[i]] <- getAllFieldFromList(theList[[i]], index, level-1)
    }
    
    return(concatenateList(field))
  }
  else{
    return (list(theList[[index]]))
  }
}

concatenateList <- function(listOfLists){
  
  n <- length(listOfLists)
  m <- lapply(listOfLists, length)
  l <- list()
  
  for(i in 1:n){
    
    for(j in 1:m[[i]]){
      
      if(i == 1)
        l[[j]] <- listOfLists[[i]][[j]]
      else
        l[[(i-1) * m[[i-1]] + j]] <- listOfLists[[i]][[j]]
    }
  }
  
  return(l)
}

# Samples a specific amount of points from a data set by using k-means algorithm----
# input:
#   x = a matrix where each row contains a point and each column is one dimension of the points
#   n = the number of samples to be collected
#   iter.max = the maximum number of iterations fot the k-mean execution
#   nstart = the number of restart for the k-mean execution
# output:
#   a vector containing the index of each point (row) selected
my.sample <- function(x, n, iter.max=100, nstart=3){
  
  begin <- getTime()
  start <- getTime()
  kmResult <- kmeans(x, n, iter.max, nstart)
  cat("kmeans executed", crono.end(start), "|")
  start <- getTime()
  
  clusters <- kmResult$cluster
  M <- length(clusters)
  
  samples <- c(0, n)
  
  for(i in 1:n){
    
    pointsOfClusterI <- which(kmResult$cluster == i)
    if(length(pointsOfClusterI) == 1)
      p <- findMatch(matrix(x[pointsOfClusterI,], nrow=1), kmResult$centers[i,])
    else
      p <- findMatch(x[pointsOfClusterI,], kmResult$centers[i,])
    samples[i] <- pointsOfClusterI[p$point]
    
    #cat(i * 100/n, "%\n")
  }
  cat("sample selected ", crono.end(start), "|")
  cat("total: ", crono.end(begin), "\n")
  
  (samples)
}

# Splits a given string 'str' by a pattern 'split' and returns ----
# the 'index'th part of the splitting.
# intput:
#   str = either a string or a vector of strings
#   split = a string representing the pattern to split the string
#   index = a integer specifying which part of the splitting shall be returned.
# output:
#   either a string or a vector of strings with the 'index'th part of the splitting.
my.strsplit <- function(str, split, index){
  
  #if 'str' is a vector ...
  if(length(str) > 1){
    #uses the standard function to split 'str'
    #'aux' will be a list will the splitting for each element of 'str'
    aux <- strsplit(str, split)
    #finds the length of the resulting splitting
    n <- length(aux)
    #initiates the result as 'n' empty strings
    res <- rep("", n)
    
    #for each of the 'n' splittings
    for(i in 1:n)
      #if index is less than 1, starts the search by the end of the splitting
      #in other words, gets the reverse position
      if(index < 1)
        res[i] <- aux[[i]][length(aux[[i]])-index]
    else
      #otherwise, gets the 'index'th position
      res[i] <- aux[[i]][index]
    
    #returns all 'n' selected parts of the splittings
    (res)
  }
  #otherwise ...
  else{
    
    #splits 'str' with the standard function
    res <- strsplit(str, split)[[1]]
    
    #if index is less than 1, gets the reverse position
    if(index < 1)
      res <- res[length(res) -index]
    else
      #otherwise, gets the 'index'th position
      res <- res[index]
    
    #returns the selected part of the splitting
    (res)
  }
}

# Prints a set of lines into an image in JPEG format.
# input:
#   lines = the set of lines, it has to be a list;
#   col = the number of columns of the image; (default 150)
#   row = the number of rows of the image; (default 210)
#   progressFlag = a boolean value, whether the progress must be printed; (default FALSE)
# output:
#   a 2D matrix containing the image
lines2image <- function(lines, col=150, row=210, progressFlag=FALSE){
  
  n <- length(lines)
  
  img <- matrix(rep(1, col * row), ncol=col)
  
  for(i in 1:n){
    
    stopFlag <- FALSE
    
    cat(lines[[i]]$inverted, "\n")
    
    if(lines[[i]]$inverted){
      y <- 1
      while(!stopFlag){
        
        x <- appLinear(lines[[i]], y)
        if(x <= col && x >= 1 && y <= row)
          img[y, x] <- 0
        else
          stopFlag <- TRUE
        
        y <- y + 1
      }
    }
    else{
      
      x <- 1
      
      while(!stopFlag){
        
        y <- appLinear(lines[[i]], x)
        if(y >= 1 && y <= row && x <= col)
          img[y, x] <- 0
        else
          stopFlag <- TRUE
        
        x <- x + 1
      }
    }
    
    if(progressFlag)
      cat(i * 100 / n, "%, ", i, "\n")
  }
  
  (img)
}

my.writeJPEG <- function(data, file, quality=1){
  
  #cat(dim(data), "\n")
  
  maxData <- max(data[which(data != 0)])
  minData <- min(data[which(data != 0)])
  
  data[which(data != 0)] <- (data[which(data != 0)] - minData)/(maxData - minData)
  
  #cat(dim(data), "\n")
  
  writeJPEG(data, file, quality)
}

resize <- function(img, file="", col=25, row=35, reason=1){
  
  if(is.character(img))
    img <- readImageData(img)
  
  colO <- length(img[1,])
  rowO <- length(img[,1])
  
  if(reason != 1){
    
    col <- colO * reason
    row <- rowO * reason
  }
  
  newImg <- matrix(rep(0, col*row), ncol=col)
  
  for(i in 1:row){
    
    y <- round(i*rowO/row)
    
    for(j in 1:col){
      
      x <- round(j*colO/col)
      
      newImg[i,j] <- img[y,x]
    }
  }
  
  if(file != ""){
    
    my.writeJPEG(newImg, file)
    write(t(newImg), concatenate(c(file, ".dat")), col)
  }
  
  (newImg)
}

# Computes the mode of a set.
# If there is no mode, it returns the mean
# input:
#   x = a vector containing the set
# output:
#   a value, either the mode or mean of x
getMode <- function(x){
  #gets all unique elements of x
  ux <- unique(x)
  #gets the frequencies of x's elements
  h <- tabulate(match(x, ux))
  #gets the values which maximizes h
  y <- ux[which(h == max(h))]
  
  #if there is more than one which maximizes h,...
  if(length(y) != 1){
    
    #computes the mean, instead
    y <- mean(x)
  }
  
  #returns the mode, or the mean
  (y)
}

# Computes the difference between each point of a line/curve
# intput:
#   data = a vector 'D' of numbers with 'n' elements
# output:
#   a vector 'V' of numbers with 'n'-1 elements containing the difference
#   between consecutives points of 'data'.
#   in other words, V[i] = |D[i] - D[i+1]|
differenceVector <- function(data){
  
  #gets the number of points of 'data'
  n <- length(data)
  
  #initializes the result with 0
  diff <- rep(0, n-1)
  
  #for 'n'-1 times...
  for(i in 1:(n-1)){
    
    #computes the absolute difference of each 2 consecutive elements of 'data'
    diff[i] <- abs(data[i] - data[i+1])
  }
  
  #returns the result
  (diff)
}

#' Converts a list containing some set of values to a matrix
list2matrix <- function(l){
  
  m <- Reduce(function(y, x){
    if(!is.matrix(y)){
      return (matrix(x, nrow=1))
    }
    return(rbind(y, matrix(x, nrow=1)))
  },  l, 0)
  
  (m)
}

#' Checks if a number has passed a limit, if it has
#' returns that limit. Otherwise, returns that number.
#' the type can be 'lower' (default) or 'upper'.
limit <- function(x, limit, type="lower"){
  
  if(type == "lower" && x < limit)
    return(limit)
  
  else if(type == "upper" && x > limit)
    return(limit)
  
  return(x)
}

getDirectory <- function(fileName){

  chars <- strsplit(fileName, "")[[1]]
  pos <- Position(function(x){return(x == "/")}, chars, right=TRUE)
  path <- substr(fileName, 1, pos)
  return(path)
}

computeNumberOfCombinations <- function(n, p){
  
#   if(p == 2)
#     return((n^2 - n)/2)
#   else
#     return(factorial(n)/(factorial(p) * factorial(n - p)))
  x <- 1
  y <- n
  for(i in 1:p){
    
    x <- x * y
    y <- y - 1
  }
  return(x/factorial(p))
}

cat.matrix <- function(m, file="", append=FALSE){
  
  if(!is.null(dim(m)))
  
    res <- apply(m, 1, function(x){
      
      cat(x, "\n", file=file, append=append)
    })
  else
    cat(paste(m, sep = " "), "\n", file=file, append=append)
}

my.as.matrix <- function(vector){
  
  n <- length(vector)
  
  return(matrix(c(1:n, vector), ncol=2))
}

strcount <- function(x, pattern, split){
  
  count <- 0;

  for(text in x){
    
    count <- count + unlist(lapply(
       strsplit(text, split),
       function(z) na.omit(length(grep(pattern, z)))
    ))
  }

  return(count)
}

addRow <- function(m, r){
  
  errorMsg <- list("The new row is neither a vector nor a matrix!\n",
                   "The number of columns are incorrect!\n")
  
  if(is.vector(r)){
    
    if(dim(m)[2] == length(r)){
      
      m <- rbind(m, matrix(r, nrow=1))
    }
    else{
      cat(errorMsg[[2]])
    }
  }
  else if(is.matrix(r)){
    
    if(dim(m)[2] == dim(r)[2]){
      
      m <- rbind(m, r)
    }
    else{
      cat(errorMsg[[2]])
    }
  }
  else{
    cat(errorMsg[[1]])
  }
  
  return(m)
}

merge.list <- function(x, y){
  
  n <- length(x)
  m <- length(y)
  
  if(m > 0){
  
    for(i in 1:m){
      
      if(is.null(names(y[i])) || names(y[i]) == ""){
        x[[i + n]] <- y[[i]]
      }
      else{
        x[[names(y[i])]] <- y[[i]]
      }
    }
  }
  
  return(x)
}

my.smooth <- function(X, d=1){
  
  N <- length(X)
  
  if(d > 0){
    X <- gaussianSmooth(c(rep(X[1], 2*d), X, rep(X[N], 2*d)), c(d))[(1 + 2*d):(N + 2*d)]
  }
  
  return(X)
}

removeNullFromList <- function(L, returnVector=FALSE){
  
  A <- list()
  N <- length(L)
  
  for(i in 1:N){
    
    if(!is.null(L[[i]]))
      A[[length(A)+1]] <- L[[i]]
  }
  
  if(returnVector)
    return(list2vector(A))
  
  return(A)
}

meanOfInterval <- function(X){
  
  return((max(X) + min(X))/2)
}

#Generates a data matrix randomically, with M samples of N variables ----
dataGenerator <- function(M, N, range=c(0:100)/100, replace=TRUE){
  
  return(matrix(sample(range, M*N, replace = replace), ncol=N))
}

#Computes the smallest angle between two line segments ----
angleBetweenSegments <- function(l1, l2){
  
  l1[[2]][1] <- l1[[2]][1] - l1[[1]][1]
  l1[[2]][2] <- l1[[2]][2] - l1[[1]][2]
  l1[[1]][1] <- 0
  l1[[1]][2] <- 0
  
  l2[[2]][1] <- l2[[2]][1] - l2[[1]][1]
  l2[[2]][2] <- l2[[2]][2] - l2[[1]][2]
  l2[[1]][1] <- 0
  l2[[1]][2] <- 0
  
  a1 <- 0
  a2 <- 0
  if(abs(l1[[2]][2]) != 0){
    a1 <- l1[[2]][1]/l1[[2]][2]
    l1[[2]][2] <- l1[[2]][2]/abs(l1[[2]][2]) * sqrt(1/(a1^2 + 1))
  }
  else{
    a1 <- l1[[2]][1]/0.0001
    l1[[2]][2] <- sqrt(1/(a1^2 + 1))
  }
  
  if(abs(l2[[2]][2]) != 0){
    a2 <- l2[[2]][1]/l2[[2]][2]
    l2[[2]][2] <- l2[[2]][2]/abs(l2[[2]][2]) * sqrt(1/(a2^2 + 1))
  }
  else{
    a2 <- l2[[2]][1]/0.0001
    l2[[2]][2] <- sqrt(1/(a2^2 + 1))
  }
  
  l1[[2]][1] <- a1 * l1[[2]][2]
  l2[[2]][1] <- a2 * l2[[2]][2]
  
  cossine <- round(l1[[2]][1] * l2[[2]][1] + l1[[2]][2] * l2[[2]][2], 8)
  res <- acos(cossine)
  
  p2 <- rotateCurve(matrix(l1[[2]], ncol=2), 0, res)
  p3 <- rotateCurve(matrix(l1[[2]], ncol=2), 0, -res)
  
  if(sqrt((p2[1] - l2[[2]][1])^2 + (p2[2] - l2[[2]][2])^2) > sqrt((p3[1] - l2[[2]][1])^2 + (p3[2] - l2[[2]][2])^2))
    res <- -res
  
  if(is.nan(res))
    cat("NAN!", c, "v1 =", l1[[2]], "v2 =", l2[[2]], "\n")
  
  return(res)
}

computeCombinations <- function(M, N){
  
  comb <- matrix(rep(0, 2), ncol=2)
  
  for(i in 1:M){
    
    if(i <= N){
      for(j in i:N){
        
        comb <- rbind(comb, matrix(c(i,j), ncol=2))
      }
    }
    else{
      for(j in 1:N){
        
        comb <- rbind(comb, matrix(c(i,j), ncol=2))
      }
    }
  }
  
  return(comb[-1,])
}

And <- function(bools){
  
  return(Reduce(function(y,x){
    
    return(y && x)
    
  }, bools, TRUE))
}

list.find <- function(l, x){
  
  M <- length(l)
  for(i in 1:M){
    
    if(identical(l[[i]],x))
      
      return(i)
  }
  
  return(0)
}