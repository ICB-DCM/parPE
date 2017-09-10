library(rhdf5)

getStarts <- function(file) {
  names(h5dump(file, load=F)$multistarts)
}

#######################################################
# Plot some optimization results
#######################################################

plotCostTrajectory <- function(cost) {
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
    newCost <- tryCatch(h5read(optimizationResultFile, paste0("/multistarts/", start, "/iterCostFunCost")), error=function(cond) {NULL})
    if(is.null(newCost))
      newCost <- as.matrix(rep(NA, nrow(cost)))
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
  status <- NULL
  for(start in starts) {
    newStatus <- tryCatch(h5read(optimizationResultFile, paste0("/multistarts/", start, "/exitStatus")), error=function(cond) {NULL})
    status <- c(status, newStatus)
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