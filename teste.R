
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
