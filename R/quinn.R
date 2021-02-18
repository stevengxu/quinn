library(splines2)
library(pracma)
library(brms)
library(rmutil)
library(fdrtool)
library(loo)
library(rlang)
source("utils.R")

quinn_samp <- function(X,z,train.model,iter,warmup=floor(iter/2),chain=1,thin=1,hyper.a=30,control=list(adapt_delta=0.999, max_td=6),seed=1)
{
  data.X <- X
  data.z <- z
  n.hidden=train.model$n.hidden
  n.knots=train.model$n.knots
  n.var=train.model$n.var
  hyper.a <- sqrt(pi/2)/hyper.a
  isp = iSpline(seq(0,1,length.out = 101), knots = seq(0,1,length.out=n.knots)[-c(1,n.knots)], degree = 2, intercept = F)
  msp = iSpline(seq(0,1,length.out = 101), knots = seq(0,1,length.out=n.knots)[-c(1,n.knots)], degree = 2, derivs = 1, intercept = F)
  quinn.env <<- new_environment(list(data.X=data.X,data.z=data.z,n.hidden=n.hidden,n.knots=n.knots,n.var=n.var,hyper.a=hyper.a,isp=isp,msp=msp))
  
  B <- vector("list",2)
  B[[1]] <- rep(0,(n.var+1)*n.hidden)
  B[[2]] <- rep(0,(n.hidden+1)*n.knots)
  logs <- rep(0,n.var+2)
  init <- c(unlist(B),logs)
  fn <- .loglik
  gr <- .loglik.grad
  
  if(!is.null(seed)) set.seed(seed)
  control <- .update_control(control)
  eps <- control$stepsize
  npar <- length(init)
  nz <- length(data.z)
  M <- control$metric
  if(is.null(M)) M <- rep(1, len=npar)
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
  theta.out <- matrix(NA, nrow=iter, ncol=npar)
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
    r.cur <- r.plus <- r.minus <-  rnorm(npar,0,1)
    H0 <- .calculate.H(theta=theta.cur, r=r.cur, fn=fn2)
    
    ## Draw a slice variable u in log space
    logu <-
      log(runif(1)) + .calculate.H(theta=theta.cur,r=r.cur, fn=fn2)
    j <- 0; n <- 1; s <- 1; divergent <- 0
    ## Track steps and divergences; updated inside .buildtree
    info <- as.environment(list(n.calls=0, divergent=0))
    while(s==1) {
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
      ## The Welford running variance calculation, see
      ## https://www.johndcook.com/blog/standard_deviation/
      if(m== w1){
        ## Initialize algorithm from end of first fast window
        m1 <- theta.out[m,]; s1 <- rep(0, len=npar); k <- 1
      } else if(m==anw){
        ## If at end of adaptation window, update the mass matrix to the estimated
        ## variances
        M <- as.numeric(s1/(k-1)) # estimated variance
        ## Update density and gradient functions for new mass matrix
        if(any(!is.finite(M))){
          warning("Non-finite estimates in mass matrix adaptation -- reverting to unit")
          M <- rep(1, length(M))
        }
        rotation <- .rotate_space(fn=fn, gr=gr, M=M,  y.cur=theta.out[m,])
        fn2 <- rotation$fn2; gr2 <- rotation$gr2; chd <- rotation$chd;
        theta.cur <- rotation$x.cur
        ## Reset the running variance calculation
        k <- 1; m1 <- theta.out[m,]; s1 <- rep(0, len=npar)
        ## Calculate the next end window. If this overlaps into the final fast
        ## period, it will be stretched to that point (warmup-w3)
        aws <- 2*aws
        anw <- .compute_next_window(m, anw, warmup, w1, aws, w3)
        ## Find new reasonable eps since it can change dramatically when M
        ## updates
        eps <- .find.epsilon(theta=theta.cur, fn=fn2, gr=gr2, eps=.01, verbose=FALSE)
        if(!is.null(control$verbose))
          print(paste(m, ": new range(M) is:",
                      round(min(M),5), round(max(M),5), ", pars",
                      which.min(M), which.max(M), ", eps=", eps))
      } else {
        k <- k+1; m0 <- m1; s0 <- s1
        ## Update M and S
        m1 <- m0+(theta.out[m,]-m0)/k
        s1 <- s0+(theta.out[m,]-m0)*(theta.out[m,]-m1)
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
  theta.out <- theta.out[seq(1, nrow(theta.out), by=thin),]
  lp__matrix <- lp__matrix[seq(1, nrow(lp__matrix), by=thin),]
  lp__ <- rowMeans(lp__matrix)
  suppressWarnings(waic <- waic(lp__matrix)$estimates[3])
  
  warm <- warmup/thin
  sampler_params <- sampler_params[seq(1, nrow(sampler_params), by=thin),]
  ndiv <- sum(sampler_params[-(1:warm),5])
  if(ndiv>0)
    message(paste0("There were ", ndiv, " divergent transitions after warmup"))
  msg <- paste0("Final acceptance ratio=", sprintf("%.2f", mean(sampler_params[-(1:warm),1])))
  if(useDA) msg <- paste0(msg,", and target=", adapt_delta)
  message(msg)
  if(useDA) message(paste0("Final step size=", round(eps, 3),
                           "; after ", warmup, " warmup iterations"))
  time.total <- difftime(Sys.time(), time.start, units='secs')
  .print.mcmc.timing(time.warmup=time.warmup, time.total=time.total)
  return(list(par=theta.out, lp__=lp__, waic=waic))
}


quinn_pred <- function(pred.model,newX,tau)
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
  isp <- iSpline(z.grid, knots = seq(0,1,length.out=n.knots)[-c(1,n.knots)], degree = 2, intercept = F)
  isp <- isp[rep(1:nrow(isp),X.shape[1]),]
  newX <- newX[rep(1:X.shape[1],each = n.z),]
  if(X.shape[2]==1) dim(newX) <- c(length(newX),1)
  cdf_est <- numeric(nrow(newX))
  for(i in 1:nsamp)
  {
    theta <- post.samp[i,]
    W <- B <- vector("list",2)
    B[[1]] <- matrix(theta[1:((n.var+1)*n.hidden)],nrow=n.var+1,ncol=n.hidden)
    B[[2]] <- matrix(theta[((n.var+1)*n.hidden+1):((n.var+1)*n.hidden+(n.hidden+1)*n.knots)],nrow=n.hidden+1,ncol=n.knots)
    logs <- theta[((n.var+1)*n.hidden+(n.hidden+1)*n.knots+1):((n.var+1)*n.hidden+(n.hidden+1)*n.knots+n.var+2)]
    s <- exp(logs)
    W[[1]] <- s[1:(n.var+1)]*B[[1]]
    W[[2]] <- s[n.var+2]*B[[2]]
    enn <- exp(.nn(newX,W,.tanh))
    cdf_est <- cdf_est + 1/nsamp*rowSums(isp*enn)/rowSums(enn)
  }
  dim(cdf_est) <- c(n.z,X.shape[1])
  q_est <- matrix(nrow = X.shape[1], ncol = length(tau))
  for(i in 1:X.shape[1])
  {
    q_est[i,] <- approx(cdf_est[,i],z.grid,xout = tau)$y 
  }
  return(q_est)
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
.nn <- function(X, W, f)
{
  fout <- .fwd(f(.fwd(X,W[[1]])),W[[2]])
  return(fout)
}

##Full log-Likelihood
.loglik <- function(theta)
{
  invisible(list2env(as.list(quinn.env),environment()))
  W <- B <- vector("list",2)
  B[[1]] <- matrix(theta[1:((n.var+1)*n.hidden)],nrow=n.var+1,ncol=n.hidden)
  B[[2]] <- matrix(theta[((n.var+1)*n.hidden+1):((n.var+1)*n.hidden+(n.hidden+1)*n.knots)],nrow=n.hidden+1,ncol=n.knots)
  logs <- theta[((n.var+1)*n.hidden+(n.hidden+1)*n.knots+1):((n.var+1)*n.hidden+(n.hidden+1)*n.knots+n.var+2)]
  s <- exp(logs)
  W[[1]] <- s[1:(n.var+1)]*B[[1]]
  W[[2]] <- s[n.var+2]*B[[2]]
  enn <- exp(.nn(data.X,W,.tanh))
  sp <- predict(msp,newx=data.z)
  loglik <- sum(log(rowSums(enn*sp)))-sum(log(rowSums(enn)))+sum(dnorm(B[[1]],log=T))+
    sum(dnorm(B[[2]],log=T))+sum(dhalfnorm(s,hyper.a,log=T))+sum(logs)
  return(loglik)
}

##Data log-Likelihood
.loglik.data <- function(theta)
{
  invisible(list2env(as.list(quinn.env),environment()))
  W <- B <- vector("list",2)
  B[[1]] <- matrix(theta[1:((n.var+1)*n.hidden)],nrow=n.var+1,ncol=n.hidden)
  B[[2]] <- matrix(theta[((n.var+1)*n.hidden+1):((n.var+1)*n.hidden+(n.hidden+1)*n.knots)],nrow=n.hidden+1,ncol=n.knots)
  logs <- theta[((n.var+1)*n.hidden+(n.hidden+1)*n.knots+1):((n.var+1)*n.hidden+(n.hidden+1)*n.knots+n.var+2)]
  s <- exp(logs)
  W[[1]] <- s[1:(n.var+1)]*B[[1]]
  W[[2]] <- s[n.var+2]*B[[2]]
  enn <- exp(.nn(data.X,W,.tanh))
  sp <- predict(msp,newx=data.z)
  loglik_array <- log(rowSums(enn*sp))-log(rowSums(enn))
  return(loglik_array)
}

##Gradient of log-Likelihood
.loglik.grad = function(theta)
{
  invisible(list2env(as.list(quinn.env),environment()))
  Z <- W <- B <- hidden <- grad_B <- denom <- numer <- vector("list",2)
  TT <- vector("list",4)
  grad_logs <- numeric(n.var+2)
  B[[1]] <- matrix(theta[1:((n.var+1)*n.hidden)],nrow=n.var+1,ncol=n.hidden)
  B[[2]] <- matrix(theta[((n.var+1)*n.hidden+1):((n.var+1)*n.hidden+(n.hidden+1)*n.knots)],nrow=n.hidden+1,ncol=n.knots)
  logs <- theta[((n.var+1)*n.hidden+(n.hidden+1)*n.knots+1):((n.var+1)*n.hidden+(n.hidden+1)*n.knots+n.var+2)]
  s <- exp(logs)
  W[[1]] <- s[1:(n.var+1)]*B[[1]]
  W[[2]] <- s[n.var+2]*B[[2]]
  hidden[[1]] <- .fwd(data.X,W[[1]])
  hidden[[2]] <- .fwd(.tanh(hidden[[1]]),W[[2]])
  Z[[1]] <- cbind(1,data.X)
  Z[[2]] <- cbind(1,.tanh(hidden[[1]]))
  sp <- predict(msp,newx = data.z)
  enn <- exp(hidden[[2]])
  snn <- sp*enn
  TT[[1]] <- (tcrossprod(snn,W[[2]][-1,])*.tanh(hidden[[1]],T))/rowSums(snn)
  TT[[2]] <- (tcrossprod(enn,W[[2]][-1,])*.tanh(hidden[[1]],T))/rowSums(enn)
  TT[[3]] <- snn/rowSums(snn)
  TT[[4]] <- enn/rowSums(enn)
  grad_B[[1]] <- crossprod(t(s[1:(n.var+1)]*t(Z[[1]])),TT[[1]])-crossprod(t(s[1:(n.var+1)]*t(Z[[1]])),TT[[2]])-B[[1]]
  grad_B[[2]] <- crossprod(s[n.var+2]*Z[[2]],TT[[3]])-crossprod(s[n.var+2]*Z[[2]],TT[[4]])-B[[2]]
  grad_logs[1:(n.var+1)] <- diag(tcrossprod(crossprod(Z[[1]],TT[[1]]),W[[1]]))-diag(tcrossprod(crossprod(Z[[1]],TT[[2]]),W[[1]])) - 2*hyper.a^2/pi*s[1:(n.var+1)]^2 + 1
  grad_logs[n.var+2] <- sum(TT[[3]]*hidden[[2]])- sum(TT[[4]]*hidden[[2]]) - 2*hyper.a^2/pi*s[n.var+2]^2 + 1
  grad_theta <- c(unlist(grad_B),grad_logs)
  return(grad_theta)
}
