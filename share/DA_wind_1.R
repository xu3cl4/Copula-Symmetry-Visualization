# set working directory 
setwd("/Users/ASUS/Desktop/projects/kaust/codes")

# source other files 
source("HT_sym_noref.R")
source("plot_curves.R")

# explore and clean the current data 
wd <- readRDS(file="Wind_Dataset1.rds")
head(wd)
dim(wd$ws)
# the dataset has no header 
ws <- t(wd$ws)
ds <- dim(ws)

# remove the effect of marginals   
n <- ds[1]
ws.nor <- pobs(as.matrix(ws))
head(ws.nor)

# run data analysis 
tests <- combn(1:3, 2)

# # scatter plots 
par(mfrow=c(2, 3))
for (i in 1:(dim(tests)[2])){
  pair <- tests[,i]
  v1 <- sprintf("Position %.0f", pair[1])
  v2 <- sprintf("Position %.0f", pair[2])
  title <-  paste(v1, " v.s ", v2)
  plot(ws.nor[, pair[1]], ws.nor[, pair[2]], pch = 16, main=title, xlab = v1, ylab = v2)
}

# visualization of curves
for (i in 1:(dim(tests)[2])){
  pair <- tests[,i]
  v1 <- sprintf("Position %0.f", pair[1])
  v2 <- sprintf("Position %0.f", pair[2])

  U <- cbind(ws.nor[, pair[1]], ws.nor[, pair[2]])

  ec <- empCopula(U)
  U.s <- simulate.sym(n, ec)

  plot.curves(f.sym, U, v1, U.s, v2)
}

# Hypothesis Testing 
for (i in 1:(dim(tests)[2])){
  pair <- tests[,i]
  print(sprintf("Positions %.0f and %.0f", pair[1], pair[2]))
        
  U <- cbind(ws.nor[, pair[1]], ws.nor[, pair[2]])
  print(size(U))
}
