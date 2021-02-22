library(tidyverse)
library(ggplot2)
library(latex2exp)
source("quinn.R")
data <- read.csv("../data/hurricane.csv")
X <- data$Year
y <- data$WmaxST

min_X <- min(X)
max_X <- max(X)
min_y <- min(y)-0.001
max_y <- max(y)+0.001

X <- (X-min_X)/(max_X-min_X)
z <- (y-min_y)/(max_y-min_y)

train.model <- list(n.hidden=5,n.knots=9,n.var=1)
iter=2000
warmup=500

mcmc <- quinn_samp(X=X,z=z,train.model=train.model,iter=iter,warmup=warmup)

post_id <- seq(iter-999,iter)
post <- mcmc$par[post_id,]

pred.model <- c(train.model,list(n.z=101,samp=post))

tau <- seq(0.05,0.95,0.05)

X.grid <- seq(0,1,length.out = 101)

q.pred <- quinn_pred(pred.model=pred.model,newX=X.grid,tau=tau)

q.pred <- q.pred*(max_y-min_y)+min_y

qplot.df = as.data.frame(cbind(X.grid*(max_X-min_X)+min_X,q.pred))

colnames(qplot.df) = c("X", paste("tau",1:length(tau),sep="."))

qplot.df = gather(qplot.df,tau,y,tau.1:tau.19,factor_key=TRUE)

ggplot()+geom_point(data = data, aes(Year,WmaxST))+labs(x = "Year", y="WmaxST")+
  geom_line(data = qplot.df, aes(X,y,col = tau))+scale_color_discrete(name = expression(tau), labels = seq(0.05,0.95,0.05))


qf.post <- matrix(nrow=length(post_id),ncol=length(tau))
for(i in 1:length(post_id))
{
  pred.model_i <- c(train.model,list(n.z=101,samp=post[i,]))
  qf.post[i,] <- quinn_pred(pred.model=pred.model_i,newX=0.8,tau=tau)
}
qf.post <- qf.post * (max_y-min_y) + min_y
qf.post.df <- as.data.frame(cbind(rep(tau,length(post_id)),c(t(qf.post))))
qf.post.df$V3 <- rep(paste0("id.",1:length(post_id)),each=length(tau))
colnames(qf.post.df) <- c("tau","qf","id")

ggplot()+geom_line(data = qf.post.df, aes(tau,qf,group=id,color=id),alpha=0.05)+theme_bw()+
  labs(x=TeX("$\\tau$"),y="",title = TeX("Posterior distribution of $Q(\\tau |Year=2001),\\tau\\in\\[0.05,0.95\\]$"))+
  theme(text = element_text(size = 18),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18))+guides(color=FALSE)
