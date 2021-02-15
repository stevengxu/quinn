library(tidyverse)
library(plotly)
source("quinn.R")

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

n.hidden=8
n.knots=11
iter=2000
warmup=500

mcmc <- quinn_samp(X=X,z=z,n.hidden=n.hidden,n.knots=n.knots,iter=iter,warmup=warmup)

post_id <- seq(iter-999,iter)
post <- mcmc$par[post_id,]

tau <- seq(0.05,0.95,0.05)

X1.grid <- seq(0,1,length.out = 51)
X2.grid <- seq(0,1,length.out = 51)
X.grid <- as.matrix(expand.grid(X1.grid,X2.grid))

q.pred <- quinn_pred(X=X.grid,param=post,tau=tau)

q.pred <- q.pred*(max_y-min_y)+min_y

V1 <- q.pred[1,]
V2 <- q.pred[10,]
V3 <- q.pred[19,]
dim(V1) <- dim(V2) <- dim(V3) <- c(51,51)

fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(x = X2.grid, y = X1.grid, z = V1) %>% 
  add_surface(x = X2.grid, y = X1.grid, z = V2) %>% 
  add_surface(x = X2.grid, y = X1.grid, z = V3)

fig