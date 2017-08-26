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

#######################################################
resultFile <- "/home/dweindl/src/parPE-build/amici/examples/steadystate/bla_rank00000.h5"
resultFile <- "/home/dweindl/Desktop/_rank00000_noprot.h5"
numStarts <- 5
#######################################################

library(rhdf5)
# library(ggplot2)

#######################################################
# Plot cost over iterations
#######################################################
# combined plot
cost <- NULL
for(start in 0:(numStarts - 1)) {
  cost <- cbind(cost, h5read(resultFile, paste0("/multistarts/", start, "/iterCostFunCost")))
}
matplot(cost, xlab="Iteration", ylab="Cost", type="l", ylim=c(4.5e3, 1e4), lty=1)
legend("bottomleft", legend = 0:(numStarts - 1), col = palette()[1:numStarts], lty=1)
# individual plots
par(mfrow=c(2,3))
for(start in 0:(numStarts - 1)) {
  cost <- h5read(resultFile, paste0("/multistarts/", start, "/iterCostFunCost"))
  plotCostTrajectory(cost)
}
#######################################################
# Plot parameters over iterations
trajectory <- h5read(resultFile, "/multistarts/0/iterCostFunParameters")
plotParameterTrajectory(trajectory)

#######################################################
# Compare to true parameters
inputFile <- "/home/dweindl/src/parPE/amici/examples/steadystate/data.h5"
inputDataAttr <- h5readAttributes(inputFile, "/data")
pTrue <- inputDataAttr$ptrue
pOpt <- h5read(resultFile, "/multistarts/0/finalParameters")
rbind(pTrue, pOpt)

#######################################################
# Plot gradient over iterations
gradient <- h5read(resultFile, "/multistarts/0/iterCostFunGradient")
matplot(gradient, type="l", xlab="Iteration", ylab="Gradient")
# last 5 iterations only 
idx <- nrow(gradient)-0:4
matplot(idx, gradient[idx, ], type="l", xlab="Iteration", ylab="Gradient")
