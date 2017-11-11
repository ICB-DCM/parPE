library(rhdf5)

getStarts <- function(file) {
  names(h5dump(file, load=F)$multistarts)
}

#######################################################
# Plot some optimization results
#######################################################

plotCostTrajectory <- function(cost) {
  #naCols <- colSums(apply(cost, 2, is.na)) == nrow(cost)
  #cost <- cost[, !naCols]
  cost[1] <- NA # skip first (high) cost value
  par(las = 1)
  plot(cost, xlab="Iteration", ylab="Cost", type="l")
  points(cost, pch=20)
}

plotParameterTrajectory <- function(trajectory) {
  matplot(trajectory, type="l", xlab="Iteration", ylab="Parameter value")
  
  legend("topright", legend=paste("Param ", 1:ncol(trajectory)), col = 1:6, lty=1:5)
}

getCostTrajectory <- function(file, starts=getStarts(file)) {
  cost <- NULL
  for(start in starts) {
    newCost <- tryCatch(h5read(file, paste0("/multistarts/", start, "/iterCostFunCost")), error=function(cond) {NULL})
    if(is.null(newCost)) {
      if(is.null(cost))
        newCost <- NA
      else
        newCost <- as.matrix(rep(NA, nrow(cost)))
    }
    maxIter <- max(nrow(cost), nrow(newCost))
    if(! is.null(cost)) {
      cost <- rbind(cost, matrix(NA, maxIter - nrow(cost), ncol(cost)))
      newCost <- c(newCost, rep(NA, maxIter - nrow(newCost)))
    }
    cost <- cbind(cost, newCost)
    colnames(cost)[ncol(cost)] <- start
  }
  return(cost)
}

getExitCodes <- function(file, starts=getStarts(file)) {
  status <- rep(NA, length(starts))
  for(i in 1:length(starts)) {
    start <- starts[i]
    newStatus <- tryCatch(h5read(file, paste0("/multistarts/", start, "/exitStatus")), error=function(cond) {NA})
    status[i] <- newStatus
    #colnames(cost)[ncol(status)] <- start
  }
  return(status)
}


removeIncrease <- function(x) {
  # replace values x_i by nan if x_i > min(x_0, ..., x_i) 
  xMin <- Inf
  for (i in 1:length(x)) {
    if(is.na(x[i]))
      next
    if(x[i] < xMin) {
      xMin <- x[i]
    } else {
      x[i] <- NA
    }
  }
  return(x)
}
#######################################################

getCorrelationByObservable <- function(simulationResultFile, starts) {
  numObservables <- dim(h5read(simulationResultFile, paste0("/multistarts/0/ySim")))[2]
  tmp <- rep(NA, numObservables * length(starts))
  corrDf <- data.frame(observable=tmp, start=tmp, corr=tmp, nonNA=tmp, cumdiff=tmp, meandiff=tmp)
  i <- 0;
  for(start in starts) {
    #cat(sprintf("Start %d:\n", start))
    for(observable in 1:numObservables) {
      ysim <- tryCatch(h5read(simulationResultFile, paste0("/multistarts/", start, "/ySim"))[,observable], 
                       error=function(cond) {NULL})
      if(is.null(ysim))
        next
      
      ymes <- h5read(simulationResultFile, paste0("/multistarts/", start, "/yMes"))[,observable]
      
      i <- i+1
      
      numNonNA <- sum(!is.na(ysim + ymes))
      
      squaredRes <- (ysim - ymes)^2
      corrDf$cumdiff[i] <- sum(squaredRes, na.rm = T)
      corrDf$meandiff[i] <- mean(squaredRes, na.rm = T)
      corrDf$observable[i] <- observable
      corrDf$start[i] <- start
      corrDf$nonNA[i] <- numNonNA
      
      hasUsableData <- numNonNA > 2
      if(!hasUsableData)
        next
      
      ct = cor.test(ymes, ysim)
      #cat(sprintf("%d: %f (%d non-NA)\n", observable, ct$estimate, numNonNA))
      
      corrDf$corr[i] <- ct$estimate
    }
  }
  corrDf$start <- as.factor(corrDf$start)
  corrDf <- corrDf[!is.na(corrDf$observable),]
  corrDf$observable <- as.factor(corrDf$observable)
  return(corrDf)
}
#######################################################
#######################################################

getParameterTrajectory <- function(file, start) {
    return(h5read(file, paste0("/multistarts/", start, "/iterCostFunParameters")))
}



getFinalCost <- function(cost) {
  final <- cost[nrow(cost), ]
  for(i in 1:length(final)) {
    final[i] <- getLastNonNA(cost[, i])
  }
  return(final)
}

getFirstNonNA <- function(x) {
  for(i in 1:length(x)) {
    if(!is.na(x[i]))
      return(x[i])
  }
  return(NA)
}

getLastNonNA <- function(x) {
  for(i in length(x):1) {
    if(!is.na(x[i]))
      return(x[i])
  }
  return(NA)
}

# test
# 1 : getLastNonNA(c(NA, 2, 1, NA))

