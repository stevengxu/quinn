##########################################################
#Code for reproducing simulation Design 1-4 in the paper.#
##########################################################
library(sn)

n <- 100

#Design 1
X <- runif(n,0,5)
y <- X + sin(2*X) + 3*rsn(n,0,1,4)


#Design 2
X <- runif(n,0,1)
fx <- 3*X
gx <- 0.5 + 2*X + sin(3*pi*X + 1)
eps <- rnorm(n)
y <- fx + gx*eps

#Deisgn 3
X <- matrix(runif(n*2),ncol=2)
fx <- sin(2*pi*X[,1]) + cos(2*pi*X[,2])
gx <- sqrt(2*(X[,1]^2+X[,2]^2))
eps <- rt(n,3)
y <- fx + gx*eps

#Design 4
X <- matrix(runif(size*10),ncol=10)
u <- runif(size)
q <- 3*(u-0.5)*(X[,1]+0.6)^3 + 15*(X[,2]+4*(X[,2]-0.5)^2)*exp(-X[,2]^2) +
  12*exp((X[,3]+0.5)^2*(X[,4]-0.5)^2) + 5*(u-1)*(X[,5]+0.4)*(X[,6]+0.5)^2 + qnorm(u,0,0.25)
