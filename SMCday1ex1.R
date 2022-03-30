# PACKAGES AND DIRECORIES------------------------------------------------------
rm(list=ls())
thiswd <- "/Users/alice/Library/CloudStorage/OneDrive-UniversityofWarwick/COURSESandCONF/SMCmasterclass"
setwd(thiswd)
library(MASS)
library(mvtnorm)


# 1 - importance weights ------------------------------------------------------

# sample from the importance distribution
rq <- function(N, D){
  y <- rmvnorm(N, mean = rep(0, D), sigma = diag(D)*2)
  return(y)
}

# importance pdf
dq <- function(y){
  dq <- dmvnorm(y, mean = rep(0, D), sigma = diag(D)*2)
  return(dq)
}

# target pdf
dp <- function(x){
  dp <- dmvnorm(x, mean = rep(0, D), sigma = diag(D))
  return(dp)
}


N=10
D=1000
Y <- rq(N=N, D=D) 
W <- rep(NA, 10)
for (n in 1:N){
  W[n] <- dp(x=Y[n, ])/dq(y=Y[n, ])
}

W
# what is happening? 
dp(x=Y[1, ])
dq(y=Y[1, ])



# 2 - importance log-weights ---------------------------------------------------

# importance lpdf
ldq <- function(y){
  ldq <- dmvnorm(y, mean = rep(0, D), sigma = diag(D)*2, log=TRUE)
  return(ldq)
}

# target lpdf
ldp <- function(x){
  ldp <- dmvnorm(x, mean = rep(0, D), sigma = diag(D), log=TRUE)
  return(ldp)
}

lW <- rep(NA, 10)
for (n in 1:N){
  lW[n] <- ldp(x=Y[n, ])-ldq(y=Y[n, ])
}

lW
# what is happening? 
ldp(x=Y[1, ])
ldq(y=Y[1, ])


# 3 - self-normalised weights ---------------------------------------------
sum(exp(lW))
Wtilde <- exp(lW)/ sum(exp(lW))

# first dimension plot
plot(Y[,1], Wtilde, type="h", lwd=4)


# 4 - max trick -----------------------------------------------------------
m <- max(lW)

t <- lW-m
sum(exp(t))
ttilde <- exp(t)/sum(exp(t))

# first dimension plot
plot(Y[,1], ttilde, type="h", lwd=4)

