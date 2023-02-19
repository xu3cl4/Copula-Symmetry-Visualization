if (!require(devtools)){
  install.packages("devtools", lib="~/R/x86_64-pc-linux-gnu-library/3.6")
}

if (!require(deSolve)){
  install_version("deSolve", version="1.28", lib="~/R/x86_64-pc-linux-gnu-library/3.6")
}

if (!require(locfit)){
  install_version("locfit", version="1.5-9.4", lib="~/R/x86_64-pc-linux-gnu-library/3.6")
}

if (!require(foreach)){
  install.packages("foreach", repos = "http://cran.us.r-project.org")
}
library(foreach)

if (!require(doParallel)){
  install.packages("doParallel", repos = "http://cran.us.r-project.org")
}
library(doParallel)

if (!require(copula)){
  install.packages("copula", repos = "http://cran.us.r-project.org")
}
library(copula)

if (!require(VineCopula)){
  install.packages("VineCopula",repos = "http://cran.us.r-project.org")
}
library(VineCopula)

if (!require(fdaoutlier)){
  install.packages("fdaoutlier", repos = "http://cran.us.r-project.org")
}
library(fdaoutlier)

#### subroutines
# resampling method
simulate.sym <- function(n, cpl){

  smpl <- rCopula(n, cpl)
  dcs <- runif(n)
  more <- (dcs > 0.5) 
  smpl[more,] <- smpl[more, c(2,1)] 

  return(smpl)
}

# testing function
f.sym <- function(U){
  function (t, u2){
    C.n(cbind(t, u2), U) - C.n(cbind(u2, t), U)
  }
}

# get data matrix for functional boxplot
get.y <- function(U, n.curves, t){
  p <- length(t)
  n <- dim(U)[1]
  if (n.curves != n){
     u1 <- runif(n.curves)
     u2 <- runif(n.curves)
  } else {
     u1 <- U[,1]  
     u2 <- U[,2]
     n.curves <- n
  }
  y <- matrix(f.sym(U)(rep(t, n.curves), rep(u2, each=p)), nrow = p)
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

# get the empirical size/power of experiments
size <- function(cpl, n, sig=0.05, p=100, N=1000, N.b=250, n.curves=NULL){

  av <- parallel::detectCores()	
  n.cores <- max(floor(av/2) + 6, av - 4)
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)

  # setup for time parameter
  e <- 10^(-2)
  t <- seq(e, 1 - e, len=p)
  n.s <- n
  
  # number of curves 
  if is.null(n.curves) {
    n.curves <- ceiling(1.2*n)
  }

  # hypothesis testing procedure
  res = foreach(a=1:N, .combine="c", .inorder=FALSE,
                .packages=c("copula", "VineCopula", "deSolve", "locfit", "fdaoutlier", "foreach"),
                .export=c("f.sym", "simulate.sym", "get.y", "get.W") ) %dopar% {
                  
	  	  U <- rCopula(n, cpl)
                  y <- get.y(U, n.curves, t)

                  # use U to construct simulated sample
                  ec <- empCopula(U)
                  U.s <- simulate.sym(n.s, ec)
                  y.s <- get.y(U.s, n.curves, t)

		  W <- get.W(y, y.s)

                  # bootstrap samples of the test statistic W
                  Ws.b = foreach(b=1:(N.b), .combine="c", .inorder=FALSE,
                                 .packages=c("copula", "VineCopula", "deSolve", "locfit", "fdaoutlier"),
                                 .export=c("f.sym", "simulate.sym", "get.y", "get.W") ) %do% {
                                   
			  	   U.b <- simulate.sym(n, ec)
                                   y.b <- get.y(U.b, n.curves, t)

                                   # use U to construct simulated sample
                                   ec.b <- empCopula(U.b)
                                   U.bs <- simulate.sym(n.s, ec.b)
                                   y.bs <- get.y(U.bs, n.curves, t)

                  		   W.b <- get.W(y.b, y.bs)               
		  }

		  p.val <- sum(Ws.b <= W)/N.b
		  return (p.val)
  		}

  stopCluster(cl)
  print(res)
  return (sum(res <= sig)/N)
}

### experiments

##### experiment setup
#n <- c(100, 250, 500, 1000)

# paramaters - test of symmetry
#delta <- c(0, 1/4, 1/2, 3/4)
delta <- c(3/4)
tau.asym <- c(0.5, 0.7, 0.9)
tau.sym <- c(1/4, 1/2, 3/4)
# paramters - test of radial/joint symmetry
tau <- c(1/4, 1/2, 3/4)

#### test simulated symmetric copula
print("#### Symmetry ####")

print("independent copula")
ic <- indepCopula()
#for (k in 1:length(n)){
#	print(sprintf("n = %f", n[k]))
#	system.time(print(size(ic, n[k])))
#}
#
### clayton
#print("clayton copula")
#for (i in 1:length(delta)){
#  for (j in 1:3){
#    for (k in 1:length(n)){
#      if (delta[i] == 0){
#        print(sprintf("delta = %f, tau = %f, n = %f", delta[i], tau.sym[j], n[k]))
#        theta <- BiCopTau2Par(3, tau.sym[j])
#        cc <- claytonCopula(theta)
#        system.time(print(size(cc, n[k])))
#	writeLines('\n')
#      } else {
#        print(sprintf("delta = %f, tau = %f, n = %f", delta[i], tau.asym[j], n[k]))
#        theta <- BiCopTau2Par(3, tau.asym[j])
#        cc <- claytonCopula(theta)
#        kcc <- khoudrajiCopula(copula1 = ic, copula2 = cc, shapes = c((1 - delta[i]), 1))
#        system.time(print(size(kcc, n[k])))
#        writeLines('\n')
#      }
#    }
#  }
#}

### gaussian
#print("gaussian copula")
#for (i in 1:length(delta)){
#  for (j in 1:3){
#    for (k in 1:length(n)){
#      if (delta[i] == 0){
#        print(sprintf("delta = %f, tau = %f, n = %f", delta[i], tau.sym[j], n[k]))
#        theta <- BiCopTau2Par(1, tau.sym[j])
#        nc <- normalCopula(theta)
#        system.time(print(size(nc, n[k])))
#        writeLines('\n')
#      } else {
#        print(sprintf("delta = %f, tau = %f, n = %f", delta[i], tau.asym[j], n[k]))
#        theta <- BiCopTau2Par(1, tau.asym[j])
#        nc <- normalCopula(theta)
#        knc <- khoudrajiCopula(copula1 = ic, copula2 = nc, shapes = c((1 - delta[i]), 1))
#        system.time(print(size(knc, n[k])))
#        writeLines('\n')
#      }
#    }
#  }
#}

## gumbel
print("gumbel copula")
for (i in 1:length(delta)){
  for (j in 1:3){
    for (k in 1:length(n)){
      if (delta[i] == 0){
      	print(sprintf("delta = %f, tau = %f, n = %f", delta[i], tau.sym[j], n[k]))
        theta <- BiCopTau2Par(4, tau.sym[j])
        gc <- gumbelCopula(theta)
        system.time(print(size(gc, n[k])))
        writeLines('\n')
      } else {
        print(sprintf("delta = %f, tau = %f, n = %f", delta[i], tau.asym[j], n[k]))
        theta <- BiCopTau2Par(4, tau.asym[j])
        gc <- gumbelCopula(theta)
        kgc <- khoudrajiCopula(copula1 = ic, copula2 = gc, shapes = c((1 - delta[i]), 1))
        system.time(print(size(kgc, n[k])))
        writeLines('\n')
      }
    }
  }
}
