library(pcalg)
library(pryr, quietly = TRUE)
library(graph, quietly = TRUE)
library(mvtnorm)
library(gtools)
source("DAG_given_ordering.R")
source("score.R")


order_sampler <- function(X, m, c2, log = FALSE, init.order = NULL, mtd = 'ztest', alpha = 0.05, threshold = 1e-1){
  
  # We transform X into matrix.
  if(is.matrix(X) != TRUE) {
    X <- data.matrix(X, rownames.force = NA)
  }
  # X should not contain NA, NaN, Inf, character, because it should be numeric.
  sanity_check_X(X)
  
  markov_sample <- list()
  score_trend <- numeric(m+1)
  G <- list()
  n <- nrow(X)
  p <- ncol(X)
  XtX = crossprod(X)
  # initialize Graph 
  if (is.null(init.order)){
    init.order <- sample(p, p)  # first try 0 matrix.
  }
  #####################
  if (mtd == 'ztest'){
    G_init = DAG_from_Ordering(X, init.order, mtd='ztest', alpha = alpha)
  }else{
    G_init = DAG_from_Ordering(X, init.order, mtd='chol', threshold = threshold)
  }
  
  
  for(t in 1:m){
    cat(t, "\rSim ")
    # sample next ordering.
    
    # sample transposition.
    A = sample(p, 2) 
    
    # adjacent transposition
    # A = sample(p-1, 2)
    # A[2] = A[1] + 1
    
    next.order = init.order
    next.order[A[1]] = init.order[A[2]]
    next.order[A[2]] = init.order[A[1]]
    if (mtd == 'ztest'){
      G_next = DAG_from_Ordering(X, next.order, mtd='ztest', alpha = alpha)
    }else{
      G_next = DAG_from_Ordering(X, next.order, mtd='chol', threshold = threshold)
    }
    
    posterior <- posterior_ratio(XtX, n, G_next, G_init, c2, c1 = 1/sqrt(2), alpha = 0.999,
                                 gamma = 1, kappa = 0, log = TRUE)
    
    if(runif(1) < exp(posterior)){
      markov_sample[[t]] = next.order
      G[[t]] = G_next
      G_init <- G_next
      init.order <- next.order
      score_trend[t+1] <- posterior + score_trend[t]
    }else{
      markov_sample[[t]] = init.order
      score_trend[t+1] <-  score_trend[t]
      G[[t]] <- G_init
    }
    
    
  }
  return(list(mc=markov_sample, score=score_trend, G=G))
}
sanity_check_X <- function(X){
  # X should not contain NA, NaN, Inf, character, because it should be numeric.
  if(any(is.character(X)) == TRUE){
    stop("X should not contain a character.")
  }
  if(any(is.nan(X)) == TRUE){
    stop("X should not contain a NaN.")
  }
  if(any(is.na(X)) == TRUE){
    stop("X should not contain a NA.")
  }
  if(any(is.infinite(X)) == TRUE){
    stop("X should not contain a infinite value.")
  }
}



'
order_sampler <- function(X, m, c2, log = FALSE, init.order = NULL){
  
  # We transform X into matrix.
  if(is.matrix(X) != TRUE) {
    X <- data.matrix(X, rownames.force = NA)
  }
  # X should not contain NA, NaN, Inf, character, because it should be numeric.
  sanity_check(X)
  
  
  markov_sample <- list()
  score_trend <- numeric(m+1)
  p <- ncol(X)
  # initialize Graph 
  if (is.null(init.order)){
    init.order <- sample(p, p)  # first try 0 matrix.
  }
  #####################
  
  G_init = DAG_from_Ordering(X, init.order, mtd = "chol")
  
  for(t in 1:m){
    cat(t, "\rSim ")
    # sample next ordering.
    
    # sample transposition.
    A = sample(p, 2) 
    
    # adjacent transposition
    # A = sample(p-1, 2)
    # A[2] = A[1] + 1
    
    next.order = init.order
    next.order[A[1]] = init.order[A[2]]
    next.order[A[2]] = init.order[A[1]]
    
    G_next <- DAG_from_Ordering(X, next.order, mtd = "chol")
    posterior <- posterior_ratio(X, G_next, G_init, c2, c1 = 1/sqrt(2), alpha = 0.999,
                                 gamma = 1, kappa = 0, log = TRUE)
    
    if(runif(1) < exp(posterior)){
      markov_sample[[t]] = next.order
      G_init <- G_next
      init.order <- next.order
      score_trend[t+1] <- posterior + score_trend[t]
    }else{
      markov_sample[[t]] = init.order
      score_trend[t+1] <-  score_trend[t]
    }
    
    
  }
  return(list(mc=markov_sample, score=score_trend))
}

'