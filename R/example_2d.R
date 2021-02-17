library(tidyverse)
library(plotly)
source("quinn.R")
source("ALEplot_qr.R")

set.seed(1)
n <- 500
X <- matrix(runif(n*2),ncol=2)
fx <- sin(2*pi*X[,1]) + cos(2*pi*X[,2])
gx <- sqrt(2*(X[,1]^2+X[,2]^2))
eps <- rnorm(n,0,1)
y <- fx + gx*eps

min_y <- min(y)-0.001
max_y <- max(y)+0.001

z <- (y-min_y)/(max_y-min_y)


train.model <- list(n.hidden=8,n.knots=11,n.var=2)
iter=2000
warmup=500

mcmc <- quinn_samp(X=X,z=z,train.model=train.model,iter=iter,warmup=warmup)

post_id <- seq(iter-999,iter)
post <- mcmc$par[post_id,]

pred.model <- c(train.model,list(n.z=101,samp=post))

tau <- seq(0.05,0.95,0.05)

X1.grid <- seq(0,1,length.out = 51)
X2.grid <- seq(0,1,length.out = 51)
X.grid <- as.matrix(expand.grid(X1.grid,X2.grid))

q.pred <- quinn_pred(pred.model=pred.model,newX=X.grid,tau=tau)

q.pred <- q.pred*(max_y-min_y)+min_y

Q1 <- q.pred[1,]
Q2 <- q.pred[10,]
Q3 <- q.pred[19,]
dim(Q1) <- dim(Q2) <- dim(Q3) <- c(51,51)

fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(x = X2.grid, y = X1.grid, z = Q1) %>% 
  add_surface(x = X2.grid, y = X1.grid, z = Q2) %>% 
  add_surface(x = X2.grid, y = X1.grid, z = Q3)

fig

pred.fun <- function(pred.model,newX,tau){
  q.pred <- quinn_pred(pred.model=pred.model,newX=newX,tau=tau)
  q.pred <- q.pred*(max_y-min_y)+min_y
  return(q.pred)
}

q.ale_V2 <- ALEPlot_qr(X=X,pred.model=pred.model,pred.fun=pred.fun,tau=tau,J=2,K=50)
