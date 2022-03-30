# PACKAGES AND DIRECORIES------------------------------------------------------
rm(list=ls())
thiswd <- "/Users/alice/Library/CloudStorage/OneDrive-UniversityofWarwick/COURSESandCONF/SMCmasterclass"
setwd(thiswd)
library(MASS)
library(mvtnorm)


# 0 - data ----------------------------------------------------------------

data <- read.csv("seOMXlogreturns2012to2014.csv", header=F)
y    <- as.vector(t(data))
plot(y, type="l", xlab = "time", ylab ="log returns")
T <- length(y)

# 1 - SSM functions -------------------------------------------------------

try_theta <- c(0.98,0.16,0.70) # phi sigma beta

rstate <- function(theta, xtm1){
  xt <- rnorm(1, mean=theta[1]*xtm1, sd=theta[2])
  return(xt)
}

dstate <- function(theta, xtm1, xt){
  dxt <- dnorm(xt, mean=theta[1]*xtm1, sd=theta[2])
  return(dxt)
}

robs   <- function(theta, xt){
  yt <- rnorm(1, mean=0, sd=theta[3]*exp(xt/2))
  return(yt)
}

dobs   <- function(theta, xt, yt){
  dyt <- dnorm(yt, mean=0, sd=theta[3]*exp(xt/2))
  return(dyt)
}

r0     <- function(theta){
  x0 <-rnorm(1, mean=0, sd=theta[2])
  return(x0)
}

x0   <- r0(theta=try_theta)
xsim <- rep(NA, T)
ysim <- rep(NA, T)

xsim[1] <- rstate(theta=try_theta, xtm1=x0)
ysim[1] <- robs(theta=try_theta, xt=xsim[1])
for (t in 2:T){
  xsim[t] <- rstate(theta=try_theta, xtm1=xsim[t-1])
  ysim[t] <- robs(theta=try_theta, xt=xsim[t])
}

# plot of the hidden states
plot(0:T, c(x0, xsim), type="l", lwd=3, col="blue", xlab="time", ylab="volatility", 
     ylim=max(abs(c(ysim, xsim)))*c(-1,1))
lines(1:T, c(ysim), type="l", lwd=1.5, col="black")
legend(0, max(abs(c(ysim, xsim))), c("X, volatility", "Y, log-returns"), 
       col=c("blue", "black"), lwd=c(3, 1.5))

# 2 - particles -----------------------------------------------------------
# Note : for R users it would be possible also to use the function Vectorize()
# which is always useful to speed up

rstate <- function(theta, xtm1, N){
  xt <- rnorm(N, mean = theta[1]*xtm1, sd = theta[2])
  return(xt)
}

dstate <- function(theta, xtm1, xt){
  dxt <- dnorm(xt, mean=theta[1]*xtm1, sd=theta[2])
  return(dxt)
}

robs   <- function(theta, xt, N){
  yt <- rnorm(N, mean=0, sd=theta[3]*exp(xt/2))
  return(yt)
}

ldobs   <- function(theta, xt, yt){
  ldyt <- dnorm(yt, mean=0, sd=theta[3]*exp(xt/2), log=T)
  return(ldyt)
}

r0     <- function(theta, N){
  x0 <-rnorm(N, mean=0, sd=theta[2])
  return(x0)
}



# 3 - BPF -----------------------------------------------------------------
# particle number
N     <- 100
Xprop   <- matrix(NA, N, T) 
X0prop  <- rep(NA, N) 
LWsamp  <- matrix(NA,N,T)
Xfilt   <- matrix(NA,N,T)

# sample
X0prop      <- r0(theta=try_theta, N=N)
Xprop[,1]   <- rstate(theta=try_theta, xtm1 = X0prop, N=N)
# weight
LWsamp[,1]  <- ldobs(theta=try_theta, xt=Xprop[,1], yt=y[1])
# max trick
m <- max(LWsamp[,1])
lwm <- LWsamp[,1]-m
# self normalised weights
snw  <-  lwm/sum(lwm)
# resample
idx <- sample(1:N, size=N, prob=snw, replace=T)
# filtered
Xfilt[, 1] <- Xprop[idx, 1] 

for (t in 2:T){
  # sample
  Xprop[,t]   <- rstate(theta=try_theta, xtm1 = Xfilt[,t-1], N=N)
  # weight
  LWsamp[,t]  <- ldobs(theta=try_theta, xt=Xprop[,t], yt=y[t])
  # max trick
  m <- max(LWsamp[,t])
  lwm <- LWsamp[,t]-m
  # self normalised weights
  snw  <-  exp(lwm)/sum(exp(lwm))
  # resample
  idx <- sample(1:N, size=N, prob=snw, replace=T)
  # filtered
  Xfilt[, t] <- Xprop[idx, t] 
}



pdf("pic2.pdf", height=5, width=9)
plot(y, type="l", xlab = "time", ylab ="volatility", col="white", xlim=c(0,30))
for (t in 1:30){
  points(rep(t-0.3, N), Xprop[,t], pch=19, col=rgb(0.5,0.5,0.5,0.1), cex=1)
  points(rep(t, N), Xprop[,t], pch=19, col=rgb(1,0,0,0.1), cex=50*(exp(LWsamp[, t])/sum(exp(LWsamp[, t]))))
  points(rep(t+0.3, N), Xfilt[,t], pch=19, col=rgb(0,0.8,0.2,0.1), cex=1)
}
lines(y)
dev.off()


# 4 - BPF likelihood approximation ----------------------------------------

BPF <- function(try_theta, N, y){
  LWsamp  <- matrix(NA,N,T)
  X0prop      <- r0(theta=try_theta, N=N)
  Xprop_t  <- rstate(theta=try_theta, xtm1 = X0prop, N=N)
  LWsamp[,1]  <- ldobs(theta=try_theta, xt=Xprop_t, yt=y[1])
  m <- max(LWsamp[,1])
  lwm <- LWsamp[,1]-m
  snw  <-  lwm/sum(lwm)
  idx <- sample(1:N, size=N, prob=snw, replace=T)
  Xfilt_t <- Xprop_t[idx] 
  
  for (t in 2:T){
    Xprop_t   <- rstate(theta=try_theta, xtm1 = Xfilt_t, N=N)
    LWsamp[,t]  <- ldobs(theta=try_theta, xt=Xprop_t, yt=y[t])
    m <- max(LWsamp[,t])
    lwm <- LWsamp[,t]-m
    snw  <-  exp(lwm)/sum(exp(lwm))
    idx <- sample(1:N, size=N, prob=snw, replace=T)
    # filtered
    Xfilt_t <- Xprop_t[idx] 
  }
  
  llike <- sum(log(apply(exp(LWsamp), 2, mean)))
  return(llike)
}

BPF(try_theta = c(0.98,0.16,0.70), N=100 , y=y)
BPF(try_theta = c(0.98,0.16,0.70), N=100 , y=y)
BPF(try_theta = c(0.98,0.16,0.70), N=100 , y=y)


