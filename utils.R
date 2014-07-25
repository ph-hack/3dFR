#' Checks for the integrety of the closest samples.
#' In other words, it checks if the 'nClosest' first closest selected
#' contains at least one sample of that class.
#' @example
#' The file 'closestDir'/cl_02463d452.txt contains the selected
#' closest for the individual 452 of the class 02463.
#' This function checks if the file contains at least one sample of the
#' class 02463.
checkClosestIntegrity <- function(closestDir, nClosest, isTraining=FALSE){
  
  closest <- dir(closestDir)
  N <- length(closest)
  count <- 0
  
  for(i in 1:N){
    minI <- 1
    maxI <- nClosest
    prefix <- "cl__"
    if(isTraining){
      minI <- 2
      maxI <- nClosest + 1
      prefix <- "cl_"
    }
    cl <- readLines(concatenate(c(closestDir, closest[i])))[minI:maxI]
    
    name <- strsplit(closest[i], prefix)[[1]][2]
    
    if(length(which(getPersonID(cl) == getPersonID(name))) == 0){
      cat(name, " has failed!\n")
      count <- count + 1
    }
  }
  
  cat("done!", count, "has failed!\n")
}

#' Goes through the log file and counts the amount
#' of time it finds the word "found" and "missed",
#' respectively, meaning that a sample was correctly
#' classified and wasn't.
count_found <- function(logFile){
  
  lines <- readLines(logFile)
  n <- length(lines)
  
  founds <- 0
  missed <- 0
  
  for(i in 1:n){
    
    line <- strsplit(lines[i], "[ ]")[[1]][1]
    
    if(line == "-"){
      
      line <- strsplit(lines[i], "[ ]")[[1]][2]
      
      if(line == "found")
        founds <- founds + 1
      else
        missed <- missed + 1
    }
  }
  
  (list(found=founds, missed=missed))
}