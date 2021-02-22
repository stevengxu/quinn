library(tidyverse)
library(plotly)
source("quinn.R")
source("ALEplot_qr.R")

#Simulate a bivariate data set
set.seed(1)
n <- 500
X <- matrix(runif(n*2),ncol=2)
fx <- sin(2*pi*X[,1]) + cos(2*pi*X[,2])
gx <- sqrt(2*(X[,1]^2+X[,2]^2))
eps <- rnorm(n,0,1)
y <- fx + gx*eps

#min-max normalize covariates and response
min_y <- min(y)-0.001 #small shifts from the boundary, can be arbitrarly small
max_y <- max(y)+0.001

z <- (y-min_y)/(max_y-min_y)

#Create model input for training by specifying number of hidden neurons, 
#number of spline knots, and dimension of the covariate vector
train.model <- list(n.hidden=8,n.knots=11,n.var=2)

#Specifying number of iterations and warmups, optionally can specify thinning
iter=2000
warmup=500

#Run MCMC
mcmc <- quinn_samp(X=X,z=z,train.model=train.model,iter=iter,warmup=warmup)

#Extract posterior samples, here they are the last 1000 iterations
post_id <- seq(iter-999,iter)
post <- mcmc$par[post_id,]

#Create model input for fitting/prediction
#n.z is the grid size for interpolating the CDF and QF
pred.model <- c(train.model,list(n.z=101,samp=post))

#Quantile levels at which predictions are sought
tau <- seq(0.05,0.95,0.05)

#X values at which predictions are sought
X1.grid <- seq(0,1,length.out = 51)
X2.grid <- seq(0,1,length.out = 51)
X.grid <- as.matrix(expand.grid(X1.grid,X2.grid))

#Predict the QFs
q.pred <- quinn_pred(pred.model=pred.model,newX=X.grid,tau=tau)

#Back-transform to original scale
q.pred <- q.pred*(max_y-min_y)+min_y

#Plot quantile surfaces, here tau = 0.05, 0.5, 0.95
Q1 <- q.pred[1,]
Q2 <- q.pred[10,]
Q3 <- q.pred[19,]
dim(Q1) <- dim(Q2) <- dim(Q3) <- c(51,51)

fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(x = X2.grid, y = X1.grid, z = Q1) %>% 
  add_surface(x = X2.grid, y = X1.grid, z = Q2) %>% 
  add_surface(x = X2.grid, y = X1.grid, z = Q3)

fig

#Wrapper prediction function to be feed into ALEPlot_qr
pred.fun <- function(pred.model,newX,tau){
  q.pred <- quinn_pred(pred.model=pred.model,newX=newX,tau=tau)
  q.pred <- q.pred*(max_y-min_y)+min_y
  return(q.pred)
}

#Estimate marginal quantile main effect of X_2 for tau = 0.05,...,0.95
q.ale_V2 <- ALEPlot_qr(X=X,pred.model=pred.model,pred.fun=pred.fun,tau=tau,J=2,K=50)

#Plot the quantile ALE estimate against X_2
plot(NULL,NULL,xlab="",ylab="",xlim=range(q.ale_V2[[1]]),ylim=range(q.ale_V2[[2]][,10]),axes=F)
axis(side = 1, at = round(seq(min(q.ale_V2[[1]]),max(q.ale_V2[[1]][,10]), length.out = 5)),cex.axis = 1.5)
axis(side = 2, at = round(seq(min(q.ale_V2[[2]][,10]),max(q.ale_V2[[2]][,10]), length.out = 5)),cex.axis = 1.5)
title(xlab=expression(italic(X)),line=2.5,cex.lab = 1.5)
title(ylab=latex2exp::TeX("$\\bar{Q}(\\tau,x)$"),line=2,cex.lab = 1.7)
grid (NULL,NULL, lty = 3)