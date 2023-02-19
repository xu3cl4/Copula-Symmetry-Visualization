if (!require(foreach)){
  install.packages("foreach")
}
library(foreach)

if (!require(doParallel)){
  install.packages("doParallel")
}
library(doParallel)

if (!require(copula)){
  install.packages("copula")
}
library(copula)

if (!require(VineCopula)){
  install.packages("VineCopula")
}
library(VineCopula)

if (!require(fdaoutlier)){
  install.packages("fdaoutlier")
}
library("fdaoutlier")

library(stats)

#### subroutines
# resampling method
simulate.jsym <- function(n, cpl){

  smpl <- rCopula(n, cpl)
  dcs <- runif(n)
  ind.1 <- (dcs > 0.25 & dcs <= 0.5)
  ind.2 <- (dcs > 0.5 & dcs <= 0.75)
  ind.3 <- (dcs > 0.75) 

  smpl[ind.1,] <- cbind(1 - smpl[ind.1, 1], smpl[ind.1, 2])
  smpl[ind.2,] <- cbind(smpl[ind.2, 1], 1 - smpl[ind.2, 2])
  smpl[ind.3,] <- 1 - smpl[ind.3, c(1,2)] 

  return(smpl)
}

# testing function
f.jsym.1 = function(U){
  function (t, u2){
    C.n(cbind(t, u2), U) + C.n(cbind(t, 1 - u2), U) - t
  }
}

f.jsym.2 = function(U){
  function (t, u2){
    C.n(cbind(t, u2), U) + C.n(cbind(1 - t, u2), U) - u2
  }
}

# get data matrix for functional boxplot
get.y.1 <- function(U, n.curves, t){
  p <- length(t)
  n <- dim(U)[1]
  if (n.curves != n){
     u2 <- runif(n.curves)
  } else {
     u2 <- U[,2]
     n.curves <- n
  }
  y <- matrix(f.jsym.1(U)(rep(t, n.curves), rep(u2, each=p)), nrow = p)
  ind <- sample(1:n.curves, floor(n.curves/2), replace=FALSE)
  y[,ind] <- (-1)*y[,ind]
  return(y)
}

get.y.2 <- function(U, n.curves, t){
  p <- length(t)
  n <- dim(U)[1]
  if (n.curves < n){
     u2 <- sample(U[,2], n.curves, replace = FALSE)
  } else {
     u2 <- U[,2]
     n.curves <- n
  }
  y <- matrix(f.jsym.2(U)(rep(t, n.curves), rep(u2, each=p)), nrow = p)
  ind <- sample(1:n.curves, floor(n.curves/2), replace=FALSE)
  y[,ind] <- (-1)*y[,ind]
  return(y)
}

# get the test statistic W
get.W <- function(y, y.s){
  n <- dim(y)[2]
  n.s <- dim(y.s)[2]

  data <- t(cbind(y, y.s))
  dps.mbd <- modified_band_depth(data)
  W.mbd <- sum(rank(dps.mbd, ties.method = "random")[1:n])

  return(W.mbd)
}

# exclude outliers in functional data
non.out <- function(data){
  n <- dim(data)[2]
  idx.out <- fbplot(data, plot=F)$outpoint
  idx <- setdiff(1:n, idx.out)
  return(data[,idx])
}

# get the empirical size/power of experiments
# N stands for the number of simulations to estimate the empirical size / power 
size <- function(cpl, n, sig=0.05, p=100, N=1000, N.b=250, n.curves=NULL){

  n.cores <- parallel::detectCores() - 5
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)

  # setup for time parameter
  e <- 10^(-2)
  t <- seq(e, 1 - e, len=p)
  n.s <- n
  if is.null(n.curves){
    n.curves <- ceiling(1.2*n)
  }

  # hypothesis testing procedure
  res = foreach(a=1:N, .combine="rbind", .inorder=FALSE,
                .packages=c("copula", "VineCopula", "deSolve", "locfit", "fdaoutlier", "foreach"),
                .export=c("f.jsym.1", "f.jsym.2", "simulate.jsym", "get.y.1", "get.y.2", "get.W", "non.out")) %dopar% {
                  
	  	  U <- rCopula(n, cpl)
          y.1 <- get.y.1(U, n.curves, t)
		  y.2 <- get.y.2(U, n.curves, t)

          # use U to construct simulated sample
          ec <- empCopula(U)
          U.s <- simulate.jsym(n.s, ec)
          y.s.1 <- get.y.1(U.s, n.curves, t)
		  y.s.2 <- get.y.2(U.s, n.curves, t)

		  W.1 <- get.W(y.1, y.s.1)
		  W.2 <- get.W(y.2, y.s.2)

          # bootstrap samples of the test statistic W
          Ws.b = foreach(b=1:(N.b), .combine="rbind", .inorder=FALSE,
                .packages=c("copula", "VineCopula", "deSolve", "locfit", "fdaoutlier"),
                .export=c("f.jsym.1", "f.jsym.2", "simulate.jsym", "get.y.1", "get.y.2", "get.W", "non.out")) %do% {
                                   
			  	   U.b <- simulate.jsym(n, ec)
                   y.b.1 <- get.y.1(U.b, n.curves, t)
				   y.b.2 <- get.y.2(U.b, n.curves, t)

                   # use U to construct simulated sample
                   ec.b <- empCopula(U.b)
                   U.bs <- simulate.jsym(n.s, ec.b)
                   y.bs.1 <- get.y.1(U.bs, n.curves, t)
				   y.bs.2 <- get.y.2(U.bs, n.curves, t)

				   W.b.1 <- get.W(y.b.1, y.bs.1)
				   W.b.2 <- get.W(y.b.2, y.bs.2)
				   return (c(W.b.1, W.b.2))	
                                 }

		  # calculate p value
          p.val.1 <- sum(Ws.b[,1] <= W.1)/N.b
		  p.val.2 <- sum(Ws.b[,2] <= W.2)/N.b
		  p.vals <- p.adjust(c(p.val.1, p.val.2), method="BY")
                  return (p.vals)
  		}

  stopCluster(cl)
  #print(res)
  
  res.min <- apply(res, 1, min)
  return (sum(res.min <= sig)/N)
}

### experiments

##### experiment setup
#n <- c(100, 250, 500, 1000)
n <- c(1000)
# paramters - test of radial/joint symmetry
tau <- c(1/4, 1/2, 3/4)


print("#### Joint Symmetry ####")
# independnet
print("independent copula")
#ic <- indepCopula()
#for (k in 1:length(n)){
#	print(sprintf("n = %f", n[k]))
#	system.time(print(size(ic, n[k])))
#}

# frank
print("frank copula")
#for (i in 1:length(tau)){
#  for (j in 1:length(n)){
#    print(sprintf("tau = %f, n = %f", tau[i], n[j]))
#    theta <- BiCopTau2Par(5, tau[i])
#    fc <- frankCopula(theta)
#    system.time(print(size(fc, n[j])))
#  }
#}

## gaussian
print("gaussian copula")
for (i in 1:length(tau)){
  for (j in 1:length(n)){
    print(sprintf("tau = %f, n = %f", tau[i], n[j]))
    theta <- BiCopTau2Par(1, tau[i])
    nc <- normalCopula(theta)
    system.time(print(size(nc, n[j])))
  }
}

# clayton
print("clayton copula")
for (i in 1:length(tau)){
  for (j in 1:length(n)){
    print(sprintf("tau = %f, n = %f", tau[i], n[j]))
    theta <- BiCopTau2Par(3, tau[i])
    cc <- claytonCopula(theta)
    system.time(print(size(cc, n[j])))
  }
}

## gumbel
print("gumbel copula")
for (i in 1:length(tau)){
  for (j in 1:length(n)){
    print(sprintf("tau = %f, n = %f", tau[i], n[j]))
    theta <- BiCopTau2Par(4, tau[i])
    gc <- gumbelCopula(theta)
    system.time(print(size(gc, n[j])))
  }
}
