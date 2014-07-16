#' Computes the angular coefficient of a line
#' given two points.
#' ----------------------------------------------
angularCoeff <- function(x1, y1, x2, y2){
  
  if(x1 - x2 != 0)
    return ((y1 - y2)/(x1 - x2))
  else
    return (0)
}

#' Computes the line model that best fits the given data points.
#' -------------------------------------------------------------
linearInterpolation <- function(data, isOpt=FALSE){
  
  m <- length(data[,1])
  
  angularCoefficients <- rep(0, m-1)
  
  if(isOpt){
    
    angularCoefficients <- mapply(angularCoeff, data[1:(m-1),1], data[1:(m-1),2], data[2:m,1], data[2:m,2])
  }
  else
    for(i in 1:(m-1)){
      
      if(data[i+1,1] - data[i,1] != 0)
        angularCoefficients[i] <- (data[i+1,2] - data[i,2])/(data[i+1,1] - data[i,1])
      else
        angularCoefficents[i] <- 0
    }
  
  #angularMean <- mean(angularCoefficients)
  aux <- which(abs(angularCoefficients) <= abs(mean(angularCoefficients)))
  if(length(aux) > 0)
    angularMean <- mean(angularCoefficients[aux])
  else
    angularMean <- mean(angularCoefficients)
#   coefficientsMean <- mean(angularCoefficients)
#   dev <- sd(angularCoefficients)
#   angularMean <- mean(angularCoefficients[intersect(which(angularCoefficients <= coefficientsMean + dev), which(angularCoefficients >= coefficientsMean - dev))])
  
  xMean <- mean(data[,1])
  yMean <- mean(data[,2])
  
  b <- yMean - angularMean * xMean
  
  (list(a = angularMean, b = b, inverted = FALSE, coefficients=angularCoefficients))
}

# the polynomial function
polynomial <- function(x, param){
  
  m <- length(x)
  n <- length(param)
  y <- rep(0, m)
  
  for(j in 1:m){
    for(i in 1:n){
    
      y[j] <- y[j] + x[j]^(i-1) * param[[i]]
    }
  }
  
  (y)
}

# this function applies a polynomial regression
# and returns the parameters that describe the polynomial function
polynomialRegression2D <- function(data, alpha, degree=3, maxIter=1000, minError=1, printProgress=FALSE){
  
  M <- length(data[,1])
  
  #data[,2] <- data[,2] * M
  
  #initializes the parameters
  parameters <- list()
  
  for(i in 1:(degree + 1)){
    
    parameters[[i]] <- sample(1:10, 1)/5
  }
  
  #initializes the maximum-iteration control variable
  k <- 1
  
  #initializes the cost function
  jCost <- c(1000)
  
  while(jCost[k] > minError && k <= maxIter){
    #initializes the cost value and the derivatives
    jCost[k] = 0.0
    
    d <- rep(0, degree + 1)
    
    #computes the cost function and the derivatives
    for(i in 1:M){
      f <- polynomial(data[i,1], parameters)
      jCost[k] = jCost[k] + (f - data[i,2])^2
      
      for(j in 1:(degree + 1))
        d[j] <- d[j] + data[i,1]^(j-1) * (f - data[i,2])
    }
    jCost[k] = jCost[k]/(2*M)
    
    #updates the parameters
    for(i in 1:(degree + 1)){
      
      d[i] <- d[i]/M
      parameters[[i]] <- parameters[[i]] - alpha * d[i]
    }
    
    #updates k e prepare the next cost value
    
    if(printProgress)
      cat("iteration ", k, "\n")
    
    if(k > 1 && jCost[k] > jCost[k-1])
      k = maxIter + 1
    
    jCost = c(jCost, jCost[k])
    k = k + 1
  }
  
  #returns the cost function for each iteration and the parameters
  (list(cost = jCost,
        param = parameters))
}