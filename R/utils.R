# Helper functions for HMC/NUTS sampler
# 
# source: https://github.com/Cole-Monnahan-NOAA/adnuts/
#
# Modifed function for finding initial step size


## Test whether a "U-turn" has occured in a branch of the binary tree. Returns TRUE if no U-turn,
## FALSE if one occurred

.test.nuts <- function(theta.plus, theta.minus, r.plus, r.minus){
  theta.temp <- theta.plus-theta.minus
  res <- (crossprod(theta.temp,r.minus) >= 0) *
    (crossprod(theta.temp, r.plus) >= 0)
  return(res)
}

## A recursive function that builds a leapfrog trajectory using a balanced
## binary tree.

.buildtree <- function(theta, r, logu, v, j, eps, H0, fn, gr,
                       delta.max=1000, info = environment() ){
  if(j==0){
    ## ## Useful code for debugging. Returns entire path to global env.
    ## if(!exists('theta.trajectory'))
    ##   theta.trajectory <<- data.frame(step=0, t(theta))
    ## base case, take one step in direction v
    r <- r+(v*eps/2)*gr(theta)
    theta <- theta+(v*eps)*r
    r <- r+(v*eps/2)*gr(theta)
    ## verify valid trajectory. Divergences occur if H is NaN, or drifts
    ## too from from true H.
    H <- .calculate.H(theta=theta, r=r, fn=fn)
    n <- logu <= H
    s <- logu < delta.max + H
    if(!is.finite(H) | s == 0){
      info$divergent <- 1; s <- 0
    }
    ## Acceptance ratio in log space: (Hnew-Hold)
    logalpha <- H-H0
    alpha <- min(exp(logalpha),1)
    info$n.calls <- info$n.calls + 1
    ## theta.trajectory <<-
    ##   rbind(theta.trajectory, data.frame(step=tail(theta.trajectory$step,1),t(theta)))
    return(list(theta.minus=theta, theta.plus=theta, theta.prime=theta, r.minus=r,
                r.plus=r, s=s, n=n, alpha=alpha, nalpha=1))
  } else {
    ## recursion - build left and right subtrees
    xx <- .buildtree(theta=theta, r=r, logu=logu, v=v, j=j-1, eps=eps,
                     H0=H0, fn=fn, gr=gr, info=info)
    theta.minus <- xx$theta.minus
    theta.plus <- xx$theta.plus
    theta.prime <- xx$theta.prime
    r.minus <- xx$r.minus
    r.plus <- xx$r.plus
    alpha <- xx$alpha
    nalpha <- xx$nalpha
    s <- xx$s
    if(!is.finite(s)) s <- 0
    nprime <- xx$n
    ## If it didn't fail, update the above quantities
    if(s==1){
      if(v== -1){
        yy <- .buildtree(theta=theta.minus, r=r.minus, logu=logu, v=v,
                         j=j-1, eps=eps, H0=H0,
                         fn=fn, gr=gr, info=info)
        theta.minus <- yy$theta.minus
        r.minus <- yy$r.minus
      } else {
        yy <- .buildtree(theta=theta.plus, r=r.plus, logu=logu, v=v,
                         j=j-1, eps=eps, H0=H0,
                         fn=fn, gr=gr, info=info)
        theta.plus <- yy$theta.plus
        r.plus <- yy$r.plus
      }
      ### Update elements:
      ## If both slice variables fail you get 0/0.
      nprime <- yy$n+ xx$n
      if(!is.finite(nprime)) {nprime <- 0}
      ## choose whether to keep this theta
      if(nprime>0)
        if(runif(1) <= yy$n/nprime)
          theta.prime <- yy$theta.prime
      alpha <- xx$alpha+yy$alpha
      nalpha <- xx$nalpha+yy$nalpha
      ## check for valid proposal
      b <- .test.nuts(theta.plus=theta.plus,
                      theta.minus=theta.minus, r.plus=r.plus,
                      
                      r.minus=r.minus)
      s <- yy$s*b
    }
    return(list(theta.minus=theta.minus, theta.plus=theta.plus,
                theta.prime=theta.prime,
                r.minus=r.minus, r.plus=r.plus, s=s, n=nprime,
                alpha=alpha, nalpha=nalpha))
  }
}

## Calculate the log joint density (Hamiltonian) value for given position and
## momentum variables.

.calculate.H <- function(theta, r, fn){
  fn(theta)-(1/2)*sum(r^2)
} 

## Estimate a reasonable starting value for epsilon (step size) for a given
## model, for use with Hamiltonian MCMC algorithms.

.find.epsilon = function(theta, fn, gr, eps = 0.01, verbose = T){
  q = theta
  epsilon = epsilon_ = eps
  a_min = 0.25
  a_cross = 0.5
  a_max = 0.75
  d = 2
  p = rnorm(length(q), 0, 1)
  current_E = .calculate.H(theta=q, r=p, fn=fn)
  p = p + epsilon * gr(q) / 2
  q = q + epsilon * p
  p = p + epsilon * gr(q) / 2
  proposed_E = .calculate.H(theta=q, r=p, fn=fn)
  diff_E = proposed_E - current_E
  direction = 2*(diff_E > log(a_cross)) - 1
  if(!is.finite(direction)) direction <- -1
  k <- 1
  while(!is.finite(diff_E) | direction*diff_E > direction*log(a_cross)){
    epsilon = epsilon_
    epsilon_ = d^direction*epsilon
    current_E = .calculate.H(theta=q, r=p, fn=fn)
    p = p + epsilon_ * gr(q) / 2
    q = q + epsilon_ * p
    p = p + epsilon_ * gr(q) / 2
    proposed_E = .calculate.H(theta=q, r=p, fn=fn)
    diff_E = proposed_E - current_E
    k = k + 1
    if(k>50) {
      stop("More than 50 iterations to find reasonable eps. Possibly bad starting value.")
    }
  }
  ep = sort(c(epsilon,epsilon_))
  epsilon = ep[1]
  epsilon_ = ep[2]
  counter = 1
  while((diff_E > log(a_max)) || (diff_E < log(a_min))){
    epsilon_m = (epsilon+epsilon_)/2
    current_E = .calculate.H(theta=q, r=p, fn=fn)
    p = p + epsilon_m * gr(q) / 2
    q = q + epsilon_m * p
    p = p + epsilon_m * gr(q) / 2
    proposed_E = .calculate.H(theta=q, r=p, fn=fn)
    diff_E = proposed_E - current_E
    if(!is.finite(diff_E)){
      epsilon = ep[1]
      break
    }
    if(diff_E > log(a_max)){
      epsilon = epsilon_m
    }else if(diff_E < log(a_min)){
      epsilon_ = epsilon_m
    }else{
      epsilon = epsilon_m
      break
    }
    counter = counter + 1    
    if(counter>1000) {
      break
    }
  }
  if(verbose) message(paste("Reasonable epsilon=", epsilon, "found after", k, "steps"))
  return(invisible(epsilon))
}

## Compute the next window size in mass matrix adaptation
.compute_next_window <- function(i, anw, warmup, w1, aws, w3){
  anw <- i+aws
  if(anw == (warmup-w3) ) return(anw)
  ## Check that the next anw is not too long. This will be the anw for the
  ## next time this is computed. If the next one is too long, extend this
  ## one to the very end.
  nwb <- anw+2*aws
  if(nwb >= warmup-w3){
    ## if(i != warmup-w3)
    ##   message(paste("Extending last slow window from", anw, "to", warmup-w3))
    anw <- warmup-w3
  }
  return(anw)
}

## Check whether adaptation is in the slow phase
.slow_phase <- function(i, warmup, w1, w3){
  ## After w1, before start of w3
  x1 <- i>= w1 # after initial fast window
  x2 <- i<= (warmup-w3) # but before last fast window
  x3 <- i < warmup # definitely not during sampling
  return(x1 & x2 & x3)
}

## Update algorithm for mass matrix.
.rotate_space <- function(fn, gr, M,  y.cur){
  ## Rotation done using choleski decomposition
  ## First case is a dense mass matrix
  if(is.matrix(M)){
    chd <- t(chol(M))               # lower triangular Cholesky decomp.
    chd.inv <- solve(chd)               # inverse
    ## Define rotated fn and gr functions
    fn2 <- function(x) fn(chd %*% x)
    gr2 <- function(x) {as.vector( gr(chd %*% x) %*% chd )}
    ## Now rotate back to "x" space using the new mass matrix M
    x.cur <- as.numeric(chd.inv %*% y.cur)
  } else if(is.vector(M)){
    chd <- sqrt(M)
    fn2 <- function(x) fn(chd * x)
    gr2 <- function(x) as.vector(gr(chd * x) ) * chd
    ## Now rotate back to "x" space using the new mass matrix M. M is a
    ## vector here. Note the big difference in efficiency without the
    ## matrix operations.
    x.cur <- (1/chd) * y.cur
  } else {
    stop("Mass matrix must be vector or matrix")
  }
  ## Redefine these functions
  ## Need to adjust the current parameters so the chain is
  ## continuous. First rotate to be in Y space.
  return(list(gr2=gr2, fn2=fn2, x.cur=x.cur, chd=chd))
}

## Update the control list.
.update_control <- function(control){
  default <- list(adapt_delta=0.65, metric=NULL, stepsize=NULL,
                  adapt_mass=TRUE, max_treedepth=6, w1=75, w2=50, w3=25)
  if(is.matrix(control$metric) & !is.null(control$adapt_mass)){
    if(control$adapt_mass==TRUE){
      warning("Mass matrix adaptation disabled if metric is a matrix")
    }
    control$adapt_mass <- FALSE
  }
  new <- default
  if(!is.null(control))
    for(i in names(control))  new[[i]] <- control[[i]]
  if(is.matrix(new$metric)) new$adapt_mass <- FALSE
  return(new)
}

## Print MCMC progress to console.
.print.mcmc.progress <- function(iteration, iter, warmup, chain){
  i <- iteration
  refresh <- max(10, floor(iter/10))
  if(i==1 | i==iter | i %% refresh ==0){
    i.width <- formatC(i, width=nchar(iter))
    out <- paste0('Chain ',chain,', Iteration: ', i.width , "/", iter, " [",
                  formatC(floor(100*(i/iter)), width=3), "%]",
                  ifelse(i <= warmup, " (Warmup)", " (Sampling)"))
    message(out)
  }
}

## Print MCMC timing to console
.print.mcmc.timing <- function(time.warmup, time.total){
  x <- ' Elapsed Time: '
  message(paste0(x, sprintf("%.1f", time.warmup), ' seconds (Warmup)'))
  message(paste0(x, sprintf("%.1f", time.total-time.warmup), ' seconds (Sampling)'))
  message(paste0(x, sprintf("%.1f", time.total), ' seconds (Total)'))
}