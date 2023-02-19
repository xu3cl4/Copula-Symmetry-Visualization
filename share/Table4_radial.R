if (!require(fda)){
  install.packages("fda")
}
library(fda)

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

#### subroutines
# resampling method
simulate.rsym <- function(n, cpl){

  smpl <- rCopula(n, cpl)
  dcs <- runif(n)
  more <- (dcs > 0.5) 
  smpl[more,] <- 1 - smpl[more, c(1,2)] 

  return(smpl)
}

# testing function
f.rsym = function(U){
  function (t, u2){
    C.n(cbind(t, u2), U) - C.n(cbind(1 - t, 1 - u2), U) + 1 - t - u2
  }
}

# get data matrix for functional boxplot
get.y <- function(U, n.curves, t){
  p <- length(t)
  n <- dim(U)[1]
  if (n.curves < n){
     u2 <- sample(U[,2], n.curves, replace = FALSE)
  } else {
     u2 <- U[,2]
     n.curves <- n
  }
  y <- matrix(f.rsym(U)(rep(t, n.curves), rep(u2, each=p)), nrow = p)
  ind <- sample(1:n.curves, floor(n.curves/2), replace=FALSE)
  y[,ind] <- (-1)*y[,ind]
  return(y)
}

get.W <- function(y, y.s){
  n <- dim(y)[2]
  n.s <- dim(y.s)[2]

  data <- t(cbind(y, y.s))
  dps <- modified_band_depth(data) 
  W <- sum(rank(dps, ties.method = "random")[1:n])
  
  return(W)
}

# exclude outliers in functional data
non.out <- function(data){
  n <- dim(data)[2]
  idx.out <- fbplot(data, plot=F)$outpoint
  idx <- setdiff(1:n, idx.out)
  return(data[,idx])
}

# get the empirical size/power of experiments
size <- function(cpl, n, sig=0.05, p=100, N=1000, N.b=250, n.curves=NULL){

  n.cores <- parallel::detectCores() - 4
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
  res = foreach(a=1:N, .combine="c",
                .packages=c("fda", "copula", "VineCopula", "deSolve", "locfit", "fdaoutlier", "foreach"),
                .export=c("f.rsym", "simulate.rsym", "get.y", "get.W")) %dopar% {
                  
	  	  U <- rCopula(n, cpl)
                  y <- get.y(U, n.curves, t)

                  # use U to construct simulated sample
                  ec <- empCopula(U)
                  U.s <- simulate.rsym(n.s, ec)
                  y.s <- get.y(U.s, n.curves, t)

		  W <- get.W(y, y.s)

                  # bootstrap samples of the test statistic W
                  Ws.b = foreach(b=1:(N.b), .combine="c",
                                 .packages=c("fda", "copula", "VineCopula", "deSolve", "locfit"),
                                 .export=c("f.rsym", "simulate.rsym", "get.y", "get.W")) %do% {
                                   
			  	   U.b <- simulate.rsym(n, ec)
                                   y.b <- get.y(U.b, n.curves, t)

                                   # use U to construct simulated sample
                                   ec.b <- empCopula(U.b)
                                   U.bs <- simulate.rsym(n.s, ec.b)
                                   y.bs <- get.y(U.bs, n.curves, t)

				   W.b <- get.W(y.b, y.bs)
				   return (W.b)	
                                 }

		  # calculate p value
                  p.val <- sum(Ws.b <= W)/N.b
                  return (p.val)
  		}

  stopCluster(cl)
  print(res)
  return (sum(res <= sig)/N)
}

### experiments

##### experiment setup
n <- c(100, 250, 500, 1000)
# paramters - test of radial/joint.rsymmetry
tau <- c(1/4, 1/2, 3/4)

#### test simulated.rsymmetric copula
ic <- indepCopula()
for (k in 1:length(n)){
	print(sprintf("n = %f", n[k]))
	system.time(print(size(ic, n[k])))
}

print("#### Radial Symmetry ####")
# frank
print("frank copula")
for (i in 1:length(tau)){
  for (j in 1:length(n)){
    print(sprintf("tau = %f, n = %f", tau[i], n[j]))
    theta <- BiCopTau2Par(5, tau[i])
    fc <- frankCopula(theta)
    system.time(print(size(fc, n[j])))
  }
}

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
