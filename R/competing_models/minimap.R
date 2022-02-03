# decomposable score P(x_i | U, {B_i, w_i})

score<- function(V, j, S, c3, c4) {
  # decomposable posterior score (logarithmic value)
  # V : crossprod(X)
  # j : a target index
  # S : a set of parent indices
  # c3, c4 : hyper-parameter
  if (length(S) == 0){
    - c4 * log(as.numeric(V[j,j]))
  }else{
    -c3 * length(S) - c4 * log(as.numeric(V[j,j] - crossprod(V[j, S],  solve(V[S, S], V[S, j]))))
  }
  
}


CI_test <- function(V, TO, n, lv=0.01){
  # fisher Z transformation.
  # n : number of samples
  # alpha : level
  p = nrow(V)
  adj = matrix(0, p, p)
  for (i in 1:(p-1)){
    for(j in (i+1):p){
      S = TO[setdiff(1:(j-1), i)]
      if(length(S) == 0){
        rho = V[TO[i],TO[j]] / sqrt(V[TO[j],TO[j]] * V[TO[i],TO[i]]) # sample corr
      }else{
        rho = as.numeric(V[TO[i],TO[j]] - crossprod(V[TO[j], S],  solve(V[S, S], V[S, TO[i]]))) /
          sqrt(as.numeric(V[TO[j],TO[j]] - crossprod(V[TO[j], S],  solve(V[S, S], V[S, TO[j]])))  *
                 as.numeric(V[TO[i],TO[i]] -crossprod(V[TO[i], S],  solve(V[S, S], V[S, TO[i]]))) ) # sample corr
      }
      z =  log((1+rho)/ (1-rho)) / 2 # fisher z-score
      if( pnorm(sqrt(n-length(S)+3) * abs(z)) > 1 - lv/2){
        adj[TO[i], TO[j]] = 1
      }
    }
  }
  return(adj)
}




# minimal I map MCMC.

minimap <- function(X, c3, c4, m, burnin, TO_init=NULL, lv=0.01){
  # V : crossprod(X)
  # c3, c4 : hyper-parameter
  # m : MC iterations
  V = crossprod(X)
  n = nrow(X)
  p = ncol(X)
  if (is.null(TO_init)){
    TO_init = sample(1:p, p)
  }else if(length(TO_init) != p){
    stop("input right dimension of init_TO")
  }
  
  p_init = 0
  DAG_sum = matrix(0, p, p)
  
  DAG_init = CI_test(V, TO_init, n, lv) # init_G
  score_list = numeric(m+1)
  score_list[1] = p_init
  #DAG_list = list()
  #order_list = list()
  #DAG_list[[1]] = DAG_init
  #order_list[[1]] = TO_init
  for (t in 1:m){
    k = sample(1:(p-1), 1)
    TO = TO_init
    TO[k] = TO_init[k+1]
    TO[k+1] = TO_init[k]
    p_next = p_init - score(V, TO_init[k], which(DAG_init[, TO_init[k]] == 1), c3, c4) -
      score(V, TO_init[k+1], which(DAG_init[, TO_init[k+1]] == 1), c3, c4)
    DAG_next = DAG_init
    if(DAG_init[TO_init[k],TO_init[k+1]]){
      DAG_next[TO_init[k],TO_init[k+1]] = 0
      DAG_next[TO_init[k+1],TO_init[k]] = 1
    }
    if(k != 1){
      # need checking.
      for (i in 1:(k-1)){
        for(j in c(k, k+1)){
          
          S = TO[setdiff(1:(j-1), i)]
          if(length(S) == 0){
            rho = V[TO[i],TO[j]] / sqrt(V[TO[j],TO[j]] * V[TO[i],TO[i]]) # sample corr
          }else{
            rho = as.numeric(V[TO[i],TO[j]] - crossprod(V[TO[j], S],  solve(V[S, S], V[S, TO[i]]))) /
              sqrt(as.numeric(V[TO[j],TO[j]] - crossprod(V[TO[j], S],  solve(V[S, S], V[S, TO[j]])))  *
                     as.numeric(V[TO[i],TO[i]] -crossprod(V[TO[i], S],  solve(V[S, S], V[S, TO[i]]))) ) # sample corr
          }
          z =  log((1+rho)/ (1-rho)) / 2 # fisher z-score
          if( pnorm(sqrt(n-length(S)+3) * abs(z)) > 1 - lv/2){
            DAG_next[TO[i], TO[j]] = 1
          }else{
            DAG_next[TO[i], TO[j]] = 0
          }
        }
      }
    }
    p_next = p_next + score(V, TO[k], which(DAG_next[, TO[k]] == 1), c3, c4) +
      score(V, TO[k+1], which(DAG_next[, TO[k+1]] == 1), c3, c4)
    accept = rbinom(1, 1, min(1, exp(p_next - p_init)))
    
    if(accept){
      p_init = p_next
      DAG_init = DAG_next
      TO_init = TO
      #DAG_list[[t+1]] = DAG_next
      #order_list[[t+1]] = TO
      score_list[t+1] = p_next
      if(t > burnin){
        DAG_sum = DAG_sum + DAG_next
      }
    }else{
      if(t > burnin){
        DAG_sum = DAG_sum + DAG_init
        #DAG_list[[t+1]] = DAG_init
        #order_list[[t+1]] = TO_init
        score_list[t+1] = p_init
      }
    }
  }
  #return(list("Ave"=DAG_sum/(m-burnin), "lst" = DAG_list, "trajactory" = score_list))
  return(list("Ave"=DAG_sum/(m-burnin),  "traj" = score_list))
}







'
#### Checking..

# model parameter.
n_rep = 30
n = 1000
p = 20
alpha = 0.99
gamma = 0.01
kappa = 0
c1 = 1
c2 = 3
c3 = (alpha*n + kappa)/2
c4 = (c2 * log(p) + 0.5 * log(1 + alpha / gamma) + log(c1))

source("utils/data_generator.R")
source("../utils/data_generator.R")

# b = 0 case mixing behavior
set.seed(20211220)
dat = data_generator_unif(n, p, 1, 1, 0.1, 0.5)
true = dat$G[[1]]
X = dat$X[[1]]
V = crossprod(X)

m = 5000

true_est = CI_test(V, 1:p, n)
est_max_score = 0
for(i in 1:p){
  est_max_score = est_max_score + score(V, i, which(true_est[, i] == 1), c3, c4)
}

true_score = 0
for(i in 1:p){
  true_score = true_score + score(V, i, which(true[, i] == 1), c3, c4)
}

init_perms = matrix(0, nrow = n_rep, ncol = p)
init_scores = numeric(n_rep)

for(j in 1:n_rep){
  TO_init <- sample(1:p, p)
  init_perms[j, ] <- TO_init
  DAG_init = CI_test(V, TO_init, n)
  for(i in 1:p){
    init_scores[j] = init_scores[j] + score(V, i, which(DAG_init[, i] == 1), c3, c4)
  }
}


for(j in 1:n_rep){
  
  minimap_estimator <- minimap(X, c3, c4, m, 0, init_perms[j, ])  
  
  if(j == 1){
    plot(1:m, minimap_estimator$traj[2:(m+1)] + init_scores[j], type="s", ylim=c(min(init_scores)*1.1, true_score*0.9) , ylab = c("log-posterior"),
         xlab=c("number of iterations"), col="gray")
    lines(1:m, y=rep(true_score,m), col="red")
  }else{
    lines(1:m, minimap_estimator$traj[2:(m+1)]+ init_scores[j], type="s", col="gray")
  }
  
}


# b = 0.3 case mixing behavior
set.seed(20211220)
dat = data_generator_misspecified(n, p, 1, 1, 0.1, 0.5, 0.3)
true = dat$G[[1]]
X = dat$X[[1]]
V = crossprod(X)

m = 50000

true_est = CI_test(V, 1:p, n)
est_max_score = 0
for(i in 1:p){
  est_max_score = est_max_score + score(V, i, which(true_est[, i] == 1), c3, c4)
}

true_score = 0
for(i in 1:p){
  true_score = true_score + score(V, i, which(true[, i] == 1), c3, c4)
}

init_perms = matrix(0, nrow = n_rep, ncol = p)
init_scores = numeric(n_rep)

for(j in 1:n_rep){
  TO_init <- sample(1:p, p)
  init_perms[j, ] <- TO_init
  DAG_init = CI_test(V, TO_init, n)
  for(i in 1:p){
    init_scores[j] = init_scores[j] + score(V, i, which(DAG_init[, i] == 1), c3, c4)
  }
}


for(j in 1:n_rep){
  
  minimap_estimator <- minimap(X, c3, c4, m, burnin, init_perms[j, ])  
  
  if(j == 1){
    plot(1:m, minimap_estimator$traj[2:(m+1)] + init_scores[j], type="s", ylim=c(min(init_scores)*1.1, true_score*0.9) , ylab = c("log-posterior"),
         xlab=c("number of iterations"), col="gray")
    lines(1:m, y=rep(true_score,m), col="red")
  }else{
    lines(1:m, minimap_estimator$traj[2:(m+1)]+ init_scores[j], type="s", col="gray")
  }
  
}


'