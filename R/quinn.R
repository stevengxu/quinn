library(splines2)
library(pracma)
library(brms)
library(rmutil)
library(fdrtool)
library(loo)
library(rlang)
source("utils.R")

quinn_samp <- function(X,
                       z,
                       train.model,
                       iter,
                       warmup=floor(iter/2),
                       chain=1,
                       thin=1,
                       hyper.a=30,
                       control=list(adapt_delta=0.999, max_td=6, metric="diag"),
                       seed=1,
                       diagnose=F)
{
  data.X <- X
  data.z <- z
  n.hidden <- train.model$n.hidden
  n.knots <- train.model$n.knots
  n.var <- train.model$n.var
  n.par <- (n.var+1)*n.hidden + (n.hidden+1)*n.knots + n.var+2
  hyper.a <- sqrt(pi/2)/hyper.a
  isp <- iSpline(seq(0,1,length.out = 101), knots = seq(0,1,length.out=n.knots)[-c(1,n.knots)], degree = 2, intercept = F)
  msp <- iSpline(seq(0,1,length.out = 101), knots = seq(0,1,length.out=n.knots)[-c(1,n.knots)], degree = 2, derivs = 1, intercept = F)
  
  ##A container of information passed to parent scope
  quinn.env <<- new_environment(list(data.X=data.X,data.z=data.z,n.hidden=n.hidden,
                                     n.knots=n.knots,n.var=n.var,n.par=n.par,
                                     hyper.a=hyper.a,isp=isp,msp=msp))
  
  nz <- length(data.z)
  init <- rep(0,n.par)
  fn <- .loglik
  gr <- .loglik.grad
  
  if(!is.null(seed)) set.seed(seed)
  control <- .update_control(control)
  eps <- control$stepsize
  M <- control$metric
  ## Initialized as unit matrix
  if(is_string(M)) M <- rep(1, len=n.par)
  if( !(is.vector(M) | is.matrix(M)) )
    stop("Metric must be vector or matrix")
  max_td <- control$max_treedepth
  adapt_delta <- control$adapt_delta
  adapt_mass <- control$adapt_mass
  ## Mass matrix adapatation algorithm arguments. Same as Stan defaults.
  w1 <- control$w1; w2 <- control$w2; w3 <- control$w3
  aws <- w2 # adapt window size
  anw <- w1+w2 # adapt next window
  if(warmup < (w1+w2+w3) & adapt_mass){
    warning("Too few warmup iterations to do mass matrix adaptation.. disabled")
    adapt_mass <- FALSE
  }
  ## Using a mass matrix means redefining what fn and gr do and
  ## backtransforming the initial value.
  rotation <- .rotate_space(fn=fn, gr=gr, M=M, y.cur=init)
  fn2 <- rotation$fn2; gr2 <- rotation$gr2
  theta.cur <- rotation$x.cur
  chd <- rotation$chd
  sampler_params <-
    matrix(numeric(0), nrow=iter, ncol=6, dimnames=list(NULL,
                                                        c("accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__",
                                                          "divergent__", "energy__")))
  lp__matrix <- matrix(0, nrow=iter, ncol=nz)
  ## This holds the rotated but untransformed variables ("y" space)
  theta.out <- matrix(NA, nrow=iter, ncol=n.par)
  ## how many steps were taken at each iteration, useful for tuning
  j.results <- rep(NA, len=iter)
  useDA <- is.null(eps)               # whether to use DA algorithm
  if(useDA){
    epsvec <- Hbar <- epsbar <- rep(NA, length=warmup+1)
    eps <- epsvec[1] <- epsbar[1] <-
      .find.epsilon(theta=theta.cur, fn=fn2, gr=gr2, eps=.01, verbose=FALSE)
    mu <- log(10*eps)
    Hbar[1] <- 0; gamma <- 0.05; t0 <- 10; kappa <- 0.75
  } else {
    ## dummy values to return
    epsvec <- epsbar <- Hbar <- NULL
  }
  ## Start of MCMC chain
  time.start <- Sys.time()
  message('')
  message(paste('Starting NUTS at', time.start))
  for(m in 1:iter){
    ## Initialize this iteration from previous in case divergence at first
    ## treebuilding. If successful trajectory they are overwritten
    theta.minus <- theta.plus <- theta.cur
    theta.out[m,] <-
      if(is.vector(M)) chd*theta.cur else t(chd %*% theta.cur)
    r.cur <- r.plus <- r.minus <-  rnorm(n.par,0,1)
    H0 <- .calculate.H(theta=theta.cur, r=r.cur, fn=fn2)
    
    ## Draw a slice variable u~U(0,exp(H0)) in log space
    logu <- log(runif(1)) + H0
    j <- 0; n <- 1; s <- 1; divergent <- 0
    ## Track steps and divergences; updated inside .buildtree
    info <- as.environment(list(n.calls=0, divergent=0))
    while(s==1) {
      #Uniformly decides direction
      v <- sample(x=c(1,-1), size=1)
      if(v==1){
        ## move in right direction
        res <- .buildtree(theta=theta.plus, r=r.plus, logu=logu, v=v,
                          j=j, eps=eps, H0=H0,
                          fn=fn2, gr=gr2, info=info)
        theta.plus <- res$theta.plus
        r.plus <- res$r.plus
      } else {
        ## move in left direction
        res <- .buildtree(theta=theta.minus, r=r.minus, logu=logu, v=v,
                          j=j, eps=eps, H0=H0,
                          fn=fn2, gr=gr2, info=info)
        theta.minus <- res$theta.minus
        r.minus <- res$r.minus
      }
      ## test whether to accept this state
      if(!is.finite(res$s)) res$s <- 0
      if(res$s==1) {
        if(runif(n=1, min=0,max=1) <= res$n/n){
          theta.cur <- res$theta.prime
          ## Rotate parameters
          theta.out[m,] <-
            if(is.vector(M)) chd*theta.cur else t(chd %*% theta.cur)
        }
      }
      n <- n+res$n
      s <- as.vector(res$s*.test.nuts(theta.plus, theta.minus, r.plus, r.minus))
      if(!is.finite(s)) s <- 0
      j <- j+1
      ## Stop doubling if too many or it's diverged enough
      if(j>=max_td) {
        ## warning("j larger than max_treedepth, skipping to next m")
        break
      }
    }
    j.results[m] <- j-1
    
    alpha2 <- res$alpha/res$nalpha
    if(!is.finite(alpha2)) alpha2 <- 0
    ## ---------------
    ## Step size adapation with the
    ## Do the adapting of eps.
    if(useDA){
      if(m <= warmup){
        ## Adaptation during warmup:
        Hbar[m+1] <- (1-1/(m+t0))*Hbar[m] +
          (adapt_delta-alpha2)/(m+t0)
        ## If logalpha not defined, skip this updating step and use
        ## the last one.
        ## if(is.nan(Hbar[m+1])) Hbar[m+1] <- abs(Hbar[m])
        logeps <- mu-sqrt(m)*Hbar[m+1]/gamma
        epsvec[m+1] <- exp(logeps)
        logepsbar <- m^(-kappa)*logeps + (1-m^(-kappa))*log(epsbar[m])
        epsbar[m+1] <- exp(logepsbar)
        eps <- epsvec[m+1]
      } else {
        ## Fix eps for sampling period
        eps <- epsbar[warmup]
      }
    }
    ## ---------------
    ## Do the adaptation of mass matrix. The algorithm is working in X
    ## space but I need to calculate the mass matrix in Y space. So need to
    ## do this coversion in the calcs below.
    if(adapt_mass & .slow_phase(m, warmup, w1, w3)){
      ## If in slow phase, update running estimate of variances
      if(m==w1){
        ## Initialize algorithm from end of first fast window
        m1 <- theta.out[m,]
        s1 <- 
          if(control$metric=="diag") rep(0, len=n.par) else diag(n.par)
        k <- 1
      } else if(m==anw){
        ## If at end of adaptation window, update the mass matrix to the estimated
        ## variances
        M <- s1/(k-1) # estimated variance
        ## Update density and gradient functions for new mass matrix
        if(any(!is.finite(M))){
          warning("Non-finite estimates in mass matrix adaptation -- reverting to unit")
          M <- rep(1, length(M))
        }
        rotation <- .rotate_space(fn=fn, gr=gr, M=M,  y.cur=theta.out[m,])
        fn2 <- rotation$fn2; gr2 <- rotation$gr2; chd <- rotation$chd;
        theta.cur <- rotation$x.cur
        ## Reset the running variance calculation
        k <- 1; m1 <- theta.out[m,]
        s1 <- 
          if(control$metric=="diag") rep(0, len=n.par) else diag(n.par)
        ## Calculate the next end window. If this overlaps into the final fast
        ## period, it will be stretched to that point (warmup-w3)
        aws <- 2*aws
        anw <- .compute_next_window(m, anw, warmup, w1, aws, w3)
        ## Find new reasonable eps since it can change dramatically when M
        ## updates
        eps <- .find.epsilon(theta=theta.cur, fn=fn2, gr=gr2, eps=.01, verbose=FALSE)
        #if(!is.null(control$verbose))
        #  print(paste(m, ": new range(M) is:",
        #              round(min(M),5), round(max(M),5), ", pars",
        #              which.min(M), which.max(M), ", eps=", eps))
      } else {
        ## The Welford online variance calculation, see
        ## https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm
        k <- k+1; m0 <- m1; s0 <- s1
        ## Update M and S
        m1 <- m0 + (theta.out[m,]-m0)/k
        s1 <- s0 + 
          if(control$metric=="diag") (theta.out[m,]-m1)*(theta.out[m,]-m0) else tcrossprod(theta.out[m,]-m1,theta.out[m,]-m0)
      }
    }
    ## End of mass matrix adaptation
    ##---------------
    
    ## Save adaptation info.
    lp__matrix[m,] <- .loglik.data(theta.out[m,])
    
    sampler_params[m,] <-
      c(alpha2, eps, j, info$n.calls, info$divergent, fn2(theta.cur))
    if(m==warmup) time.warmup <- difftime(Sys.time(), time.start, units='secs')
    .print.mcmc.progress(m, iter, warmup, chain)
  } ## end of MCMC loop
  
  ## Process the output for returning
  ndiv <- sum(sampler_params[-(1:warmup),5])
  if(ndiv>0)
    message(paste0("There were ", ndiv, " divergent transitions after warmup"))
  msg <- paste0("Final acceptance ratio=", sprintf("%.2f", mean(sampler_params[-(1:warmup),1])))
  
  theta.out <- theta.out[seq(warmup+1, iter, by=thin),]
  lp__matrix <- lp__matrix[seq(warmup+1, iter, by=thin),]
  lp__ <- rowMeans(lp__matrix)
  suppressWarnings(waic <- waic(lp__matrix)$estimates[3])
  sampler_params <- sampler_params[seq(warmup+1, iter, by=thin),]
  
  if(useDA) msg <- paste0(msg,", and target=", adapt_delta)
  message(msg)
  if(useDA) message(paste0("Final step size=", round(eps, 3),
                           "; after ", warmup, " warmup iterations"))
  time.total <- difftime(Sys.time(), time.start, units='secs')
  .print.mcmc.timing(time.warmup=time.warmup, time.total=time.total)
  if(!diagnose) return(list(par=theta.out, lp__=lp__, waic=waic)) else return(list(par=theta.out, lp__=lp__, waic=waic, diagnose=sampler_params))
}


quinn_pred <- function(pred.model,newX,tau,type="qf")
{
  n.hidden <- pred.model$n.hidden
  n.knots <- pred.model$n.knots
  n.var <- pred.model$n.var
  n.z <- pred.model$n.z
  post.samp <- pred.model$samp
  if(is.null(n.z)) n.z <- 101
  nsamp <- nrow(post.samp)
  if(is.null(nsamp))
  {
    nsamp <- 1
    dim(post.samp) <- c(1,length(post.samp))
  } 
  X.shape <- dim(newX)
  if(is.null(X.shape)) dim(newX) <- X.shape <- c(length(newX),1)
  z.grid <- seq(0,1,length.out = n.z)
  is_pdf <- type == "pdf"
  pred.sp <- iSpline(z.grid, knots = seq(0,1,length.out=n.knots)[-c(1,n.knots)], degree = 2, intercept = F, derivs = is_pdf)
  pred.sp <- pred.sp[rep(1:nrow(pred.sp),X.shape[1]),]
  newX <- newX[rep(1:X.shape[1],each = n.z),]
  if(X.shape[2]==1) dim(newX) <- c(length(newX),1)
  pred.df <- numeric(nrow(newX))
  for(i in 1:nsamp)
  {
    theta <- post.samp[i,]
    loc.par <- .location.par(n.hidden,n.knots,n.var)
    
    B1 <- matrix(theta[loc.par$B1],nrow=n.var+1,ncol=n.hidden)
    B2 <- matrix(theta[loc.par$B2],nrow=n.hidden+1,ncol=n.knots)
    logs1 <- theta[loc.par$logs1]
    logs2 <- theta[loc.par$logs2]
    
    enn <- exp(.nn(newX,exp(logs1)*B1,exp(logs2)*B2,.tanh))
    pred.df <- pred.df + 1/nsamp*rowSums(pred.sp*enn)/rowSums(enn)
  }
  dim(pred.df) <- c(n.z,X.shape[1])
  if(type!="qf") return(pred.df)
  pred.qf <- matrix(nrow = X.shape[1], ncol = length(tau))
  for(i in 1:X.shape[1])
  {
    pred.qf[i,] <- approx(pred.df[,i],z.grid,xout = tau)$y 
  }
  return(pred.qf)
}

##Hyperbolic tanh with gradient
.tanh = function(x, grad = FALSE){
  if(!grad) return(base::tanh(x))
  else return(1-base::tanh(x)^2)
}

##Feed-forward
.fwd <- function(input, weight) 
{
  output <- cbind(1,input)%*%weight
  return(output)
}

##Neural-network
.nn <- function(X, W1, W2, f)
{
  fout <- .fwd(f(.fwd(X,W1)),W2)
  return(fout)
}

##Full log-Likelihood
.loglik <- function(theta)
{
  invisible(list2env(as.list(quinn.env),environment()))
  loc.par <- .location.par(n.hidden,n.knots,n.var)
  
  B1 <- matrix(theta[loc.par$B1],nrow=n.var+1,ncol=n.hidden)
  B2 <- matrix(theta[loc.par$B2],nrow=n.hidden+1,ncol=n.knots)
  logs1 <- theta[loc.par$logs1]
  logs2 <- theta[loc.par$logs2]
  logs <- c(logs1,logs2)
  
  enn <- exp(.nn(data.X,exp(logs1)*B1,exp(logs2)*B2,.tanh))
  sp <- predict(msp,newx=data.z)
  loglik <- sum(log(rowSums(enn*sp)))-sum(log(rowSums(enn)))+sum(dnorm(B1,log=T))+
    sum(dnorm(B2,log=T))+sum(dhalfnorm(exp(logs),hyper.a,log=T))+sum(logs)
  return(loglik)
}

##Data log-Likelihood
.loglik.data <- function(theta)
{
  invisible(list2env(as.list(quinn.env),environment()))
  loc.par <- .location.par(n.hidden,n.knots,n.var)
  
  B1 <- matrix(theta[loc.par$B1],nrow=n.var+1,ncol=n.hidden)
  B2 <- matrix(theta[loc.par$B2],nrow=n.hidden+1,ncol=n.knots)
  logs1 <- theta[loc.par$logs1]
  logs2 <- theta[loc.par$logs2]
  
  enn <- exp(.nn(data.X,exp(logs1)*B1,exp(logs2)*B2,.tanh))
  sp <- predict(msp,newx=data.z)
  loglik_array <- log(rowSums(enn*sp))-log(rowSums(enn))
  return(loglik_array)
}

##Gradient of log-Likelihood
.loglik.grad = function(theta)
{
  invisible(list2env(as.list(quinn.env),environment()))
  grad_theta <- numeric(n.par)
  
  loc.par <- .location.par(n.hidden,n.knots,n.var)
  
  B1 <- matrix(theta[loc.par$B1],nrow=n.var+1,ncol=n.hidden)
  B2 <- matrix(theta[loc.par$B2],nrow=n.hidden+1,ncol=n.knots)
  logs1 <- theta[loc.par$logs1]
  logs2 <- theta[loc.par$logs2]
  
  hidden1 <- .fwd(data.X,exp(logs1)*B1)
  hidden2 <- .fwd(.tanh(hidden1),exp(logs2)*B2)
  Z1 <- cbind(1,data.X)
  Z2 <- cbind(1,.tanh(hidden1))
  sp <- predict(msp,newx = data.z)
  enn <- exp(hidden2)
  snn <- sp*enn
  
  TT1 <- (tcrossprod(snn,exp(logs2)*B2[-1,])*.tanh(hidden1,T))/rowSums(snn)
  TT2 <- (tcrossprod(enn,exp(logs2)*B2[-1,])*.tanh(hidden1,T))/rowSums(enn)
  TT3 <- snn/rowSums(snn)
  TT4 <- enn/rowSums(enn)
  
  grad_theta[loc.par$B1] <- crossprod(t(exp(logs1)*t(Z1)),TT1) - 
    crossprod(t(exp(logs1)*t(Z1)),TT2) - B1
  
  grad_theta[loc.par$B2] <- crossprod(exp(logs2)*Z2,TT3) - 
    crossprod(exp(logs2)*Z2,TT4) - B2
  
  grad_theta[loc.par$logs1] <- diag(tcrossprod(crossprod(Z1,TT1),exp(logs1)*B1)) - 
    diag(tcrossprod(crossprod(Z1,TT2),exp(logs1)*B1)) - 
    2*hyper.a^2/pi*exp(logs1)^2 + 1
  
  grad_theta[loc.par$logs2] <- sum(TT3*hidden2)- sum(TT4*hidden2) - 
    2*hyper.a^2/pi*exp(logs2)^2 + 1
  
  return(grad_theta)
}

##Return location of parameters in the sample
.location.par <- function(n.hidden,n.knots,n.var)
{
  B1 <- 1:((n.var+1)*n.hidden)
  B2 <- ((n.var+1)*n.hidden+1):((n.var+1)*n.hidden+(n.hidden+1)*n.knots)
  logs1 <- ((n.var+1)*n.hidden+(n.hidden+1)*n.knots+1):((n.var+1)*n.hidden+(n.hidden+1)*n.knots+n.var+1)
  logs2 <- (n.var+1)*n.hidden+(n.hidden+1)*n.knots+n.var+2
  return(list(B1=B1,B2=B2,logs1=logs1,logs2=logs2))
}
