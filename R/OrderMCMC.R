## Code refactorization.
## Frequent check, consistency.

## Major changes
## 1. c3, c4 inside the function as function arguments.
## 2. clear explanation of functions : input and output.

EPS = 1e-12

## variable selection X_q ~ X_a;
stepwise <- function(V, q, a, sum_sse_rest, c3, c4, max_num_par){
  # V : the gram matrix of crossprod(X), p x p matrix.
  # q : index of target node.
  # a : a set of independent node variable indices.
  # sum_sse_rest : sum of sse of other nodes.
  # c3, c4 : model parameters.
  # max_num_par : maximum number of parents.
  
  m = length(a) # the number of potential parent nodes
  par = NULL
  sel = integer(m) # 1 if selected, 0 otherwise.
  num_par = 0
  score = -c3 * log(sum_sse_rest + V[q, q])
  sse.q = V[q, q]
  while (1){
    sses = numeric(m)
    for (i in 1:m){
      if (sel[i] == 1){
        sses[i] = sse.q
      }else{
        pa = c(par, a[i])
        sses[i] = V[q,q] - crossprod(V[q, pa],  solve(V[pa, pa], V[pa, q]))
      }
    }
    best_score = -(length(par) + 1) * c4 - c3 * log(sum_sse_rest + min(sses))
    if (best_score > score + EPS){
      add = which.min(sses)[1]
      sel[add] = 1
      par = c(par, a[add])
      score = best_score
      sse.q = min(sses)
      num_par = num_par + 1
      if (num_par >= max_num_par){
        break
      }
      #	cat("add", add, score, sse.q, "\n")
    }else{
      break
    }
  }
  while (length(par) >= 2){
    s = length(par)
    sses = numeric(s)
    for (i in 1:s){
      pa = par[-i]
      sses[i] = V[q,q] - crossprod(V[q, pa],  solve(V[pa, pa], V[pa, q]))
    }
    best_score = -(length(par) - 1) * c4 - c3 * log(sum_sse_rest + min(sses))
    if (best_score > score + EPS){
      del = which.min(sses)[1]
      par = par[-del]
      score = best_score
      sse.q = min(sses)
      #	cat("del", del, score, sse.q, "\n")
    }else{
      break
    }
  }
  # Output
  # model: influential variables of X_q among X_a
  # score: the score of model given the rest of sse.
  # sse: sum of square of X_q from the model.
  return (list("model"=par, "score"=score, "sse"=sse.q))
}

## TopDown algorithm using empirical Bayesian scoring criterion
top_down <- function(V, sse_pvector, c3, c4, max_num_par){
  # V : the gram matrix of crossprod(X), p x p matrix.
  # sse_pvector : p vector of sses.
  # c3, c4 : model parameters
  # max_num_par : maximum number of parents.
  
  p = ncol(V)
  perm = which.min(diag(V))[1]
  mm = max(diag(V))
  sel = numeric(p)
  sel[perm] = 1
  G = matrix(0, nrow=p, ncol=p)
  cat("1\r")
  for (t in 2:p){
    cat(t,"\r")
    sses = numeric(p)
    G1 = matrix(0, nrow=p, ncol=p)
    for (i in 1:p){
      if (sel[i] == 1){
        sses[i] = mm + 1
        next
      }
      s = c(perm, i)
      out = stepwise(V, i, perm, sum(sse_pvector[-i]), c3, c4, max_num_par)
      sses[i] = out$sse
      G1[out$model, i] = 1
    }
    node = which.min(sses)[1]
    sel[node] = 1
    sse_pvector[node] = sses[node]
    perm = c(perm, node)
    G[,node] = G1[,node]
  }
  # Output
  # TO: topological ordering.
  # SSE: p vector of sses.
  # DAG: estimated DAG with respect to TO.
  return(list("TO"=perm, "SSE"=sse_pvector, "DAG"=G))
}

## since Bayesian scoring criterion is not decomposable, we apply TD recursively
top_down_iterative <- function(V, c3, c4, max_num_par){
  # V : the gram matrix of crossprod(X), p x p matrix.
  # c3, c4 : model parameters
  # max_num_par : maximum number of parents.
  
  p = ncol(V)
  cat("BayTD iteration 1\n")
  td = top_down(V, diag(V), c3, c4, max_num_par)
  perm = numeric(p)
  count = 1
  while ( sum(abs(perm - td$TO)) > 0.5 ){
    count = count + 1
    cat("BayTD iteration", count, "\n")
    perm = td$TO
    td = top_down(V, td$SSE, c3, c4, max_num_par)
  }
  # td : top_down output
  return(td)
}

## similar to stepwise; returns detailed candidate parent nodes information
stepwise_range <- function(V, q, a, min_sse_pvector, max_sse_pvector, c3, c4, max_num_par){
  # V : the gram matrix of crossprod(X), p x p matrix.
  # q : index of target node.
  # a : a set of independent node variable indices. (potential parents)
  # min_sse_pvector : p vector of lower bound sse of p nodes.
  # max_sse_pvector : p vector of upper bound sse of p nodes.
  # c3, c4 : model parameters
  # max_num_par : maximum number of parents.
  
  p = ncol(V)
  m = length(a)
  if (m == 0){
    return (list("parents"=NULL, "remove"=NULL, "dSSE"=NULL, "SSE"=V[q,q]))
  }
  par = NULL # parent set among a.
  sel = integer(m) # 0-1 vector for the selection. 
  num_par = 0
  sse = sum(min_sse_pvector[-q]) # lower bound for sse sums except the q-th variable. (For the forward selection, we want to pick parents as many as possible)
  score = -c3 * log(sse + V[q, q]) # starting from the null parent.
  sse.q = V[q, q]
  while (1){
    sses = numeric(m) # the i-th element indicates sse.q when the i-th variable is included as a parent.
    for (i in 1:m){
      if (sel[i] == 1){
        sses[i] = sse.q # if the element is already selected as a parent, then return the previously updated sse.q. (with the current par)
      }else{
        pa = c(par, a[i]) # include the i-th element in the parent set. 
        sses[i] = V[q,q] - crossprod(V[q, pa],  solve(V[pa, pa], V[pa, q])) 
      }
    }
    best_score = -(length(par) + 1) * c4 - c3 * log(sse + min(sses)) # the best score among m possible element to include.
    if (best_score > score + EPS){
      add = which.min(sses)[1] # the index of the variable to include.
      sel[add] = 1 # update the 0-1 vector.
      num_par = num_par + 1
      par = c(par, a[add]) # update the index set.
      score = best_score # total score
      sse.q = min(sses) # the sse of the q-th element (target).
      if(num_par >= max_num_par){
        break
      }
    }else{
      break # no forward inclusion lower the score. (Done with the forward phase.)
    }
  }
  # the backward phase.
  cpar = par # the forward phase output. (! This is what we save.)
  mp = length(cpar) # the number of the variables from the forward phase
  par_id = integer(p) # parent id. 
  for (i in 1:mp){ par_id[par[i]] = i }
  sse = sum(max_sse_pvector[-q])
  score = -length(par) * c4 - c3 * log(sse + sse.q)
  rml = NULL
  dsse = NULL
  csse = sse.q

  while (length(par) >= 1){
    s = length(par)
    sses = numeric(s)
    for (i in 1:s){
      pa = par[-i]
      if (length(pa) == 0){
        sses[i] = V[q,q]
      }else{
        sses[i] = V[q,q] - crossprod(V[q, pa],  solve(V[pa, pa], V[pa, q]))
      }
    }
    best_score = -(length(par) - 1) * c4 - c3 * log(sse + min(sses))
    if (best_score > score + EPS){
      del = which.min(sses)[1]
      rml = c(rml, par_id[par[del]])
      dsse = c(dsse, min(sses) - sse.q)
      par = par[-del]
      score = best_score
      sse.q = min(sses)
    }else{
      break
    }
  }



  return (list("parents"=cpar, "remove"=rml, "dSSE"=dsse, "SSE"=csse))
}

## returns the maximum DAG (parent sets) given the ordering
candidate_parents <- function(V, TO, min_sse_pvector, max_sse_pvector, c3, c4, max_num_par){
  # V : the gram matrix of crossprod(X), p x p matrix.
  # TO : topological ordering (p vector)
  # min_sse_pvector : p vector of lower bound sse of p nodes.
  # max_sse_pvector : p vector of upper bound sse of p nodes.
  # c3, c4 : model parameters
  # max_num_par : maximum number of parents.
  
  p = ncol(V)
  parents = vector(mode="list", length=p)
  for (j in 1:p){
    i = TO[j]
    a = NULL
    if (j >= 2){a = TO[1:(j-1)]}
    parents[[i]] = stepwise_range(V, i, a, min_sse_pvector, max_sse_pvector, c3, c4, max_num_par)
  }

  # parents : the candidate parent sets of each node (p element list)
  return(parents)
}

## find an approximation of the best DAG given the candidate parent sets.
order_to_DAG <- function(parents, c3, c4){
  # parents : the candidate parent sets of each node (p element list)
  # c3, c4 : model parameters

  p = length(parents)
  sum_sse = 0
  n_edge = 0
  sse_pvector_withsum = numeric(p+1)
  for (i in 1:p){
    sum_sse = sum_sse + parents[[i]]$SSE
    sse_pvector_withsum[i] = parents[[i]]$SSE
    n_edge = n_edge + length(parents[[i]]$parents)
  }
  score = -n_edge * c4 - c3 * log(sum_sse) # We start from the (nodewise) forward selection output.
  #cat("Initial score", score, "\n")
  # We perform the global-wise backward selection.
  rem = integer(p) # each element indicates the number of removed parents.
  while (1){
    min_d = Inf
    x = NULL
    for (i in 1:p){
      if (length(parents[[i]]$remove) <= rem[i]){next} # if there is no variable to remove.
      diff = parents[[i]]$dSSE[rem[i]+1] # The sse decrease of the next parent removal of i-th element.
      if (diff < min_d){
        min_d = diff
        x = i # the index whose parent is subject to be removed.
      }
    }
    if (is.null(x)){break}
    new_score = -(n_edge - 1)*c4 - c3 * log(sum_sse + min_d)
    if (new_score > score + EPS){
      #cat(x, parents[[x]]$remove[rem[x] + 1], score, new_score, "\n")
      rem[x] = rem[x] + 1
      sse_pvector_withsum[x] = sse_pvector_withsum[x] + min_d
      sum_sse = sum_sse + min_d
      n_edge = n_edge - 1
      score = new_score
    }else{
      break
    }
  }
  sse_pvector_withsum[p+1] = sum_sse
  G = matrix(0, p, p)
  for (i in 1:p){
    pa = parents[[i]]$parents
    if (rem[i] > 0){
      r = parents[[i]]$remove
      pa = pa[ -r[1:(rem[i])] ]
    }
    if (length(pa) > 0){
      G[pa,i] = 1
    }
  }

  # DAG : adjacency matrix of given topological ordering.
  # score : the score of DAG.
  # SSE : p+1 vector of c(sse of p nodes, sum of p sses).
  return(list("DAG"=G, "score"=score, "SSE"=sse_pvector_withsum))
}


rao_black_estimate <- function(V, TO, sse_pvector_withsum, DAG, c3, c4){
  # V : the gram matrix of crossprod(X), p x p matrix.
  # TO : topological ordering (p vector)
  # sse_pvector_withsum : p+1 vector of c(sse of p nodes, sum of p sses).
  # DAG : adjacency matrix of given DAG.
  # c3, c4 : model parameters
  p = ncol(V)
  G_hat = matrix(0, p, p)
  score = -sum(DAG) * c4 - c3 * log(sse_pvector_withsum[p+1])
  for (j in TO[2:p]){ # child
    pa = which(DAG[, j] == 1) # parent of j in G
    for(i in TO[1:(which(TO == j) - 1)]){ #i: potential parent of j
      if(DAG[i, j] == 1){
        pa_j = pa[!(pa %in% i)]
        if(length(pa_j) == 0){
          sse_next_j = V[j,j]
        }else{
          sse_next_j = V[j,j] - crossprod(V[j, pa_j],  solve(V[pa_j, pa_j], V[pa_j, j]))
        }
        new_score = -(sum(DAG) - 1) * c4 - c3 * log(sse_pvector_withsum[p+1] - sse_pvector_withsum[j] + sse_next_j)
        G_hat[i, j] = 1 / (exp(new_score-score) + 1)
      }else{
        pa_j = c(pa, i)
        sse_next_j = V[j,j] - crossprod(V[j, pa_j],  solve(V[pa_j, pa_j], V[pa_j, j]))
        new_score = -(sum(DAG) + 1) * c4 - c3 * log(sse_pvector_withsum[p+1] - sse_pvector_withsum[j] + sse_next_j)
        G_hat[i, j] = 1 / (1 + exp(score - new_score))
      }
    }

  }
  # G_hat : rao blackwellization of input DAG.
  return(G_hat)
}

create_min_sse_pvector <- function(V, max_num_par){
  p = ncol(V)
  min_sse_pvector = numeric(p)
  for(q in 1:p){
    # regress Xi on the rest. 
    # do the screening
    # cat(q, "\r")
    a = (1:p)[-q]
    par = NULL
    # we will just use screening
    sses = numeric(p-1)
    for (i in 1:(p-1)){
      pa = c(par, a[i])
      sses[i] = V[q,q] - crossprod(V[q, pa],  solve(V[pa, pa], V[pa, q]))
    }
    pa = a[order(sses)][1:max_num_par]
    min_sse_pvector[q] = V[q,q] - crossprod(V[q, pa],  solve(V[pa, pa], V[pa, q]))
  }
  return(min_sse_pvector)
}

order_mcmc <- function(V, B, N, c3, c4, max_num_par = ncol(V), highdim = F){
  # V: gram matrix of data, crossprod(X).
  # B: the number of burnin step.
  # N: the number of total MC iterations.
  # c3, c4 : model parameters
  # max_num_par : maximum number of parents.
  if(highdim == T){
    if(max_num_par >= ncol(V)){
      stop("input valid max_num_par less than n.")
    }
  }
  
  p = ncol(V)
  #print(paste("The total number of MCMC iterations:", N))
  #print(paste("The number of burn-in iterations:", B))
  
  init = top_down_iterative(V, c3, c4, max_num_par)
  perm = init$TO
  ITD = init$DAG # output 1. ITD
  # perm = sample(1:p, p)
  # ITD = matrix(0, p, p)
  
  if(highdim == T){
    min_sse_pvector = create_min_sse_pvector(V, max_num_par)
  }else{
    min_sse_pvector = 1/diag(solve(V))
  }
  
  max_sse_pvector = diag(V)
  sup_parents = candidate_parents(V, perm, min_sse_pvector, max_sse_pvector, c3, c4, max_num_par)
  dag = order_to_DAG(sup_parents, c3, c4)

  # simple average.
  G = matrix(0, p, p)
  # RB estimate.
  G_hat_sum = matrix(0, p, p)
  #sse = sse_DAG(V, dag$DAG)
  counts = integer(3)
  #cat(0, dag$score, sum(dag$DAG),  0, 0, "\n")
  for (t in 1:N){
    cat(t,"\r")
    if(t == (B+1)){
      sse_pvector_withsum = dag$SSE
      score = -sum(dag$DAG) * c4 - c3 * log(sse_pvector_withsum[p+1])
      G_hat = rao_black_estimate(V, perm, sse_pvector_withsum, dag$DAG, c3, c4)
    }

    a = sample(1:(p-1), 1)
    b = a + 1
    new_perm = perm
    new_perm[a] = perm[b]
    new_perm[b] = perm[a]

    pa_a = NULL
    if (a > 1){pa_a = new_perm[1:(a-1)]}
    pa_b = new_perm[1:(b-1)]
    parents_a = stepwise_range(V, new_perm[a], pa_a, min_sse_pvector, max_sse_pvector, c3, c4, max_num_par)
    parents_b = stepwise_range(V, new_perm[b], pa_b, min_sse_pvector, max_sse_pvector, c3, c4, max_num_par)
    dd = HD(parents_a$parents, sup_parents[[new_perm[a]]]$parents) + HD(parents_b$parents, sup_parents[[new_perm[b]]]$parents)
    if (dd < EPS){
      perm = new_perm

      if(t > B){
        G_hat[perm[a], perm[b]] = 0
        j = perm[a]
        G_pa_a = which(dag$DAG[,j] == 1)
        pa_j = c(G_pa_a, perm[b])
        sse_next_j = V[j,j] - crossprod(V[j, pa_j],  solve(V[pa_j, pa_j], V[pa_j, j]))
        new_score = -(sum(dag$DAG) + 1) * c4 - c3 * log(sse_pvector_withsum[p+1] - sse_pvector_withsum[j] + sse_next_j)
        G_hat[perm[b], perm[a]] = 1 / (1 + exp(score - new_score))
      }

      #cat(t, "S", dag$score, sum(dag$DAG), perm[a], perm[b], G_hat[perm[b], perm[a]], new_score, score,  sse_next_j, sse[j], sse[p+1], c3, c4, sum(dag$DAG), "\n", file='tmp.log', append=TRUE)
      #cat(t, "S", dag$score, sum(dag$DAG), perm[a], perm[b], G_hat[perm[b], perm[a]], new_score, score,  "\n", file='tmp.log', append=TRUE)
      counts[1] = counts[1] + 1
    }else{
      new_parents = sup_parents
      new_parents[[new_perm[a]]] = parents_a
      new_parents[[new_perm[b]]] = parents_b
      new_dag = order_to_DAG(new_parents, c3, c4)
      ll = new_dag$score - dag$score
      if (runif(1) < exp(ll)){
        #cat(parents_a$parents, " and ", sup_parents[[new_perm[a]]]$parents, "\n")
        #cat(parents_b$parents, " and ", sup_parents[[new_perm[b]]]$parents, "\n")

        perm = new_perm
        sup_parents = new_parents
        dag = new_dag

        if(t > B){
          #print(cbind(sse, sse_DAG(V, dag$DAG)))
          sse_pvector_withsum = new_dag$SSE
          score = -sum(dag$DAG) * c4 - c3 * log(sse_pvector_withsum[p+1])
          G_hat = rao_black_estimate(V, perm, sse_pvector_withsum, dag$DAG, c3, c4)
        }

        ############# this is slightly different from my previous code; I am only printing an iteration if the LL changes #####
        if (abs(ll) > EPS){
          counts[2] = counts[2] + 1
          #cat(t, new_dag$score, sum(new_dag$DAG),  perm[a], perm[b], "\n")
        }else{
          counts[1] = counts[1] + 1
        }
        #cat(t, "A", new_dag$score, sum(new_dag$DAG),  perm[a], perm[b], "\n", file='tmp.log', append=TRUE)
      }else{
        counts[3] = counts[3] + 1
        #cat(t, "R", new_dag$score, sum(new_dag$DAG),  perm[a], perm[b], "\n", file='tmp.log', append=TRUE)
      }
    }
    #print(G_hat[3,1])
    if(t > B){
      G_hat_sum = G_hat_sum + G_hat
    }
    G = G + dag$DAG
  }
  #cat("MCMC summary", counts, "\n")

  # RB : rao blackwellization version of DAG
  # Ave : MC sample average of DAG
  # ITD : DAG from ITD.
  return(list("RB"=G_hat_sum/(N-B), "Ave"=G/N, "ITD" = ITD))
}

order_mcmc_mixing <- function(V, N, method, TO, c3, c4, max_num_par = ncol(V), highdim = F){
  # V : the gram matrix of crossprod(X), p x p matrix.
  # N: the number of total MC iterations.
  # method : neighborhood relation in c('adj_trs', 'ran_trs','ran-to-ran')
  # TO : topological ordering (p vector)
  # c3, c4 : model parameters

  p = ncol(V)
  perm = TO
  min_sse_pvector = 1/diag(solve(V))
  max_sse_pvector = diag(V)
  sup_parents = candidate_parents(V, perm, min_sse_pvector, max_sse_pvector, c3, c4, max_num_par)
  dag = order_to_DAG(sup_parents, c3, c4)
  G = matrix(0, p, p)
  counts = integer(3)
  score_trajectory = numeric(N)
  if((method=='ran-to-ran')|(method=='ran_trs')){
    for(t in 1:N){
      if(method=='ran_trs'){
        ran_trs_idx = sample(1:p, 2, replace = FALSE)
        new_perm = perm
        new_perm[ran_trs_idx[1]] = perm[ran_trs_idx[2]]
        new_perm[ran_trs_idx[2]] = perm[ran_trs_idx[1]]
      }else{
        ran_to_ran_idx = sample(1:p, 2, replace = FALSE) # first index element to second index.
        a = perm[ran_to_ran_idx[1]]
        b = perm[-ran_to_ran_idx[1]]
        if(ran_to_ran_idx[2] == 1){
          new_perm = c(a, b)
        }else if(ran_to_ran_idx[2] == p){
          new_perm = c(b, a)
        }else{
          new_perm = c(b[1:(ran_to_ran_idx[2]-1)], a, b[(ran_to_ran_idx[2]):(p-1)])
        }
      }

      new_parents = candidate_parents(V, new_perm, min_sse_pvector, max_sse_pvector, c3, c4, max_num_par)
      new_dag = order_to_DAG(new_parents, c3, c4)
      if(sum(abs(new_dag$DAG - dag$DAG)) < EPS){
        perm = new_perm
        #	cat(t, "S", dag$score, sum(dag$DAG), perm[a], perm[b], "\n", file='tmp.log', append=TRUE)
        score_trajectory[t] = dag$score
        counts[1] = counts[1] + 1
      }else{
        ll = new_dag$score - dag$score
        if (runif(1) < exp(ll)){
          perm = new_perm
          dag = new_dag
          counts[2] = counts[2] + 1
          #cat(t, "A", new_dag$score, sum(new_dag$DAG),  perm[a], perm[b], "\n", file='tmp.log', append=TRUE)
          score_trajectory[t] = new_dag$score
        }else{
          counts[3] = counts[3] + 1
          #	cat(t, "R", new_dag$score, sum(new_dag$DAG),  perm[a], perm[b], "\n", file='tmp.log', append=TRUE)
          score_trajectory[t] = dag$score
        }
      }
      G = G + dag$DAG
    }
  }else{
    for (t in 1:N){
      a = sample(1:(p-1), 1)
      b = a + 1
      new_perm = perm
      new_perm[a] = perm[b]
      new_perm[b] = perm[a]

      pa_a = NULL
      if (a > 1){pa_a = new_perm[1:(a-1)]}
      pa_b = new_perm[1:(b-1)]
      parents_a = stepwise_range(V, new_perm[a], pa_a, min_sse_pvector, max_sse_pvector, c3, c4, max_num_par)
      parents_b = stepwise_range(V, new_perm[b], pa_b, min_sse_pvector, max_sse_pvector, c3, c4, max_num_par)
      dd = HD(parents_a$parents, sup_parents[[new_perm[a]]]$parents) + HD(parents_b$parents, sup_parents[[new_perm[b]]]$parents)
      if (dd < EPS){
        perm = new_perm
        #	cat(t, "S", dag$score, sum(dag$DAG), perm[a], perm[b], "\n", file='tmp.log', append=TRUE)
        score_trajectory[t] = dag$score
        counts[1] = counts[1] + 1
      }else{
        new_parents = sup_parents
        new_parents[[new_perm[a]]] = parents_a
        new_parents[[new_perm[b]]] = parents_b
        new_dag = order_to_DAG(new_parents, c3, c4)
        ll = new_dag$score - dag$score
        if (runif(1) < exp(ll)){
          perm = new_perm
          sup_parents = new_parents
          dag = new_dag
          counts[2] = counts[2] + 1
          #cat(t, "A", new_dag$score, sum(new_dag$DAG),  perm[a], perm[b], "\n", file="tmp.log", append=TRUE)
          score_trajectory[t] = new_dag$score
        }else{
          counts[3] = counts[3] + 1
          #	cat(t, "R", new_dag$score, sum(new_dag$DAG),  perm[a], perm[b], "\n", file='tmp.log', append=TRUE)
          score_trajectory[t] = new_dag$score
        }
      }
      G = G + dag$DAG
    }

  }


  # cat("MCMC summary", counts, "\n")


  # Ave : MC sample average of DAG
  # traj : score trajectory of MC samples.
  return(list('Ave'=G/N, 'traj'=score_trajectory))
}


# Hamming distance
HD <- function(x, y){
  # x, y : same sized vector.
  if (length(x) != length(y)){return(100)}
  if (is.null(x)){return (0)}
  x = x[order(x)]
  y = y[order(y)]
  #	cat("HD1", x, "\n")
  #	cat("HD2", y, "\n")
  return(sum(abs(x-y)))
}


score_DAG <- function(V, G, c3, c4){
  sse = 0
  for (i in 1:p){
    pa = which(G[,i]==1)
    if (length(pa) == 0){
      sse = sse + V[i,i]
    }else{
      sse = sse + V[i,i] - crossprod(V[i, pa],  solve(V[pa, pa], V[pa, i]))
    }
  }
  score = -sum(G) * c4 - c3 * log(sse)
  return(score)
}




'
## this function is only used for sanity check.


sse_DAG <- function(V, G){
  sse = 0
  sses = numeric(p+1)
  for (i in 1:p){
    pa = which(G[,i]==1)
    if (length(pa) == 0){
      sses[i] = V[i,i]
      sse = sse + sses[i]
    }else{
      sses[i] =  V[i,i] - crossprod(V[i, pa],  solve(V[pa, pa], V[pa, i]))
      sse = sse + sses[i]
    }
  }
  sses[p+1] = sse
  return(sses)
}


order_mcmc <- function(V, N){
  init = top_down_iterative(V)
  perm = init$TO
  min_sse = 1/diag(solve(V))
  max_sse = diag(V)
  sup_parents = candidate_parents(V, perm, min_sse, max_sse)
  dag = order_to_DAG(sup_parents)
  G = matrix(0, p, p)
  counts = integer(3)
  for (t in 1:N){
    a = sample(1:(p-1), 1)
    b = a + 1
    new_perm = perm
    new_perm[a] = perm[b]
    new_perm[b] = perm[a]

    pa_a = NULL
    if (a > 1){pa_a = new_perm[1:(a-1)]}
    pa_b = new_perm[1:(b-1)]
    parents_a = stepwise_range(V, new_perm[a], pa_a, min_sse, max_sse)
    parents_b = stepwise_range(V, new_perm[b], pa_b, min_sse, max_sse)
    dd = HD(parents_a$parents, sup_parents[[new_perm[a]]]$parents) + HD(parents_b$parents, sup_parents[[new_perm[b]]]$parents)
    if (dd < EPS){
      perm = new_perm
      #	cat(t, "S", dag$score, sum(dag$DAG), perm[a], perm[b], "\n", file="tmp.log", append=TRUE)
      counts[1] = counts[1] + 1
    }else{
      new_parents = sup_parents
      new_parents[[new_perm[a]]] = parents_a
      new_parents[[new_perm[b]]] = parents_b
      new_dag = order_to_DAG(new_parents)
      ll = new_dag$score - dag$score
      if (runif(1) < exp(ll)){
        perm = new_perm
        sup_parents = new_parents
        dag = new_dag
        counts[2] = counts[2] + 1
        cat(t, "A", new_dag$score, sum(new_dag$DAG),  perm[a], perm[b], "\n", file="tmp.log", append=TRUE)
      }else{
        counts[3] = counts[3] + 1
        #	cat(t, "R", new_dag$score, sum(new_dag$DAG),  perm[a], perm[b], "\n", file="tmp.log", append=TRUE)
      }
    }
    G = G + dag$DAG
  }
  cat("MCMC summary", counts, "\n")
  return(G/N)
}


order_mcmc_RB <- function(V, N){
  init = top_down_iterative(V)
  perm = init$TO
  min_sse = 1/diag(solve(V))
  max_sse = diag(V)
  sup_parents = candidate_parents(V, perm, min_sse, max_sse)
  dag = order_to_DAG(sup_parents)
  sse = dag$SSE

  G = matrix(0, p, p)
  # RB estimate
  G_hat_sum = matrix(0, p, p)
  #sse = sse_DAG(V, dag$DAG)
  score = -sum(dag$DAG) * c4 - c3 * log(sse[p+1])
  G_hat = rao_black_estimate(dag$DAG, V, perm, sse)

  counts = integer(3)
  cat(0, dag$score, sum(dag$DAG),  0, 0, "\n")
  for (t in 1:N){
    a = sample(1:(p-1), 1)
    b = a + 1
    new_perm = perm
    new_perm[a] = perm[b]
    new_perm[b] = perm[a]

    pa_a = NULL
    if (a > 1){pa_a = new_perm[1:(a-1)]}
    pa_b = new_perm[1:(b-1)]
    parents_a = stepwise_range(V, new_perm[a], pa_a, min_sse, max_sse)
    parents_b = stepwise_range(V, new_perm[b], pa_b, min_sse, max_sse)
    dd = HD(parents_a$parents, sup_parents[[new_perm[a]]]$parents) + HD(parents_b$parents, sup_parents[[new_perm[b]]]$parents)
    if (dd < EPS){
      perm = new_perm
      G_hat[perm[a], perm[b]] = 0
      j = perm[a]
      G_pa_a = which(dag$DAG[,j] == 1)
      pa_j = c(G_pa_a, perm[b])
      sse_next_j = V[j,j] - crossprod(V[j, pa_j],  solve(V[pa_j, pa_j], V[pa_j, j]))
      new_score = -(sum(dag$DAG) + 1) * c4 - c3 * log(sse[p+1] - sse[j] + sse_next_j)
      G_hat[perm[b], perm[a]] = 1 / (1 + exp(score - new_score))

      #cat(t, "S", dag$score, sum(dag$DAG), perm[a], perm[b], G_hat[perm[b], perm[a]], new_score, score,  sse_next_j, sse[j], sse[p+1], c3, c4, sum(dag$DAG), "\n", file="tmp.log", append=TRUE)
      #cat(t, "S", dag$score, sum(dag$DAG), perm[a], perm[b], G_hat[perm[b], perm[a]], new_score, score,  "\n", file="tmp.log", append=TRUE)
      counts[1] = counts[1] + 1
    }else{
      new_parents = sup_parents
      new_parents[[new_perm[a]]] = parents_a
      new_parents[[new_perm[b]]] = parents_b
      new_dag = order_to_DAG(new_parents)
      ll = new_dag$score - dag$score
      if (runif(1) < exp(ll)){
        #cat(parents_a$parents, " and ", sup_parents[[new_perm[a]]]$parents, "\n")
        #cat(parents_b$parents, " and ", sup_parents[[new_perm[b]]]$parents, "\n")

        perm = new_perm
        sup_parents = new_parents
        dag = new_dag
        sse = new_dag$SSE
        #print(cbind(sse, sse_DAG(V, dag$DAG)))
        score = -sum(dag$DAG) * c4 - c3 * log(sse[p+1])
        G_hat = rao_black_estimate(dag$DAG, V, perm, sse)

        ############# this is slightly different from my previous code; I am only printing an iteration if the LL changes #####
        if (abs(ll) > EPS){
          counts[2] = counts[2] + 1
          cat(t, new_dag$score, sum(new_dag$DAG),  perm[a], perm[b], "\n")
        }else{
          counts[1] = counts[2] + 1
        }
        #cat(t, "A", new_dag$score, sum(new_dag$DAG),  perm[a], perm[b], "\n", file="tmp.log", append=TRUE)
      }else{
        counts[3] = counts[3] + 1
        #cat(t, "R", new_dag$score, sum(new_dag$DAG),  perm[a], perm[b], "\n", file="tmp.log", append=TRUE)
      }
    }
    #print(G_hat[3,1])
    G_hat_sum = G_hat_sum + G_hat
    G = G + dag$DAG
  }
  cat("MCMC summary", counts, "\n")
  return(list("RB"=G_hat_sum/N, "Ave"=G/N))
}


'
