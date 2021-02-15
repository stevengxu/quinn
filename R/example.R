source("quinn.R")
source("utils.R")
data <- read.csv("../data/hurricane.csv")
X <- data$Year
y <- data$WmaxST

min_X <- min(X)
max_X <- max(X)
min_y <- min(y)-0.001
max_y <- max(y)+0.001

X <- (X-min_X)/(max_X-min_X)
z <- (y-min_y)/(max_y-min_y)

n.hidden=5
n.knots=9
iter=2000
warmup=500

mcmc <- quinn_samp(X=X,z=z,n.hidden=n.hidden,n.knots=n.knots,iter=iter,warmup=warmup)

post_id <- seq(iter-999,iter)
post <- mcmc$par[post_id,]

tau <- seq(0.05,0.95,0.05)

X.grid <- seq(0,1,length.out = 101)
z.grid <- seq(0,1,length.out = 101)

q.pred <- quinn_pred(X=X.grid,z=z.grid,param=post,tau=tau)

q.pred <- q.pred*(max_y-min_y)+min_y

qplot.df = as.data.frame(cbind(X.grid*(max_X-min_X)+min_X,t(q.pred)))

colnames(qplot.df) = c("X", paste("tau",1:length(tau),sep="."))

qplot.df = gather(plot.df,tau,y,tau.1:tau.19,factor_key=TRUE)

ggplot()+geom_point(data = data, aes(Year,WmaxST))+labs(x = "Year", y="WmaxST")+
  geom_line(data = qplot.df, aes(X,y,col = tau))+scale_color_discrete(name = expression(tau), labels = seq(0.05,0.95,0.05))

