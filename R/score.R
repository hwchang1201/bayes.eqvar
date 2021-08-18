posterior_ratio <- function(XtX, n, G, H, c2, c1 = 1/sqrt(2), alpha = 0.999,
                            gamma = 1, kappa = 0, log = FALSE){  # Input
  # XtX : crossprod(X) - computation purpose.
  # n : nrows(X)
  # G : graph or adjacent matrix.(only valid DAG)
  # H : graph or adjacent matrix.(only valid DAG)
  
  # hyperparameters not yet defined.
  
  # We transform X into matrix and G into adjacency matrix.
  if(is.matrix(G) != TRUE) {
    G = as(G, "matrix")
  }
  if(is.matrix(H) != TRUE) {
    H = as(H, "matrix")
  }
  # X should not contain NA, NaN, Inf, character.
  # G should be valid DAG -> ith column shows parent of ith node.
  # H should be valid DAG -> ith column shows parent of ith node.
  # sanity_check(G, H)
  
  p = ncol(XtX)
  # by posterior modularity, we only need to check only small number of nodes.
  # If we use add/delete/swap, we have only one node to check with difference 1(add/delete), difference 2(swap)
  # In this algorithm we only use add/delete, we have only one node to check with difference 1(add/delete).
  
  if (sum(H != G) == 0){
    val = 0
  } else{
    eta = c1 * p^(c2) * sqrt(1 + alpha / gamma)
    val = (sum(H) - sum(G)) * log(eta)
    val = val - (alpha * p * n + kappa)/2 * (log(sum_SSE(XtX,G)) - log(sum_SSE(XtX,H)))
  }
  
  # log = FALSE : then exponentiate the value.
  if(log == FALSE){
    val = exp(val)
  }
  # OUTPUT : 
  # ratio = posterior probability ratio pi_n(G)/pi_n(H) .
  # sum_SSE_G : sum of SSE by graph G.
  return(val)
}



sum_SSE <-function(XtX, G){
  p = ncol(XtX)
  sumSSE = 0
  for(i in 1:p){
    PA = which(G[, i] == 1)
    if(sum(PA) == 0){
      sumSSE = sumSSE + XtX[i,i]
    }else{
      sumSSE = sumSSE + XtX[i,i] - crossprod(XtX[i, PA], solve(XtX[PA, PA], XtX[PA, i]))
    }
  }
  as.numeric(sumSSE)
}


sanity_check <- function(G, H){
  
  # G should be valid DAG.
  if(isValidGraph(G, type = "dag") == FALSE){
    stop("the first graph is not a valid DAG.")
  }
  # H should be valid DAG.
  if(isValidGraph(H, type = "dag") == FALSE){
    stop("the second graph is not a valid DAG.")
  }
}