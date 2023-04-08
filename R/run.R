rm(list = ls())
library(foreach, quietly = TRUE)
library(doParallel, quietly = TRUE)
library(doRNG, quietly = TRUE)
library(iterators, quietly = TRUE)
library(pcalg, quietly = TRUE)
library(optparse, quietly = TRUE)
library(ggplot2)
library(mvtnorm)
library(gtools)
source("competing_models/GDS/startups/startupGDS.R", chdir = TRUE)
source("competing_models/EqVarDAG/startupEqVarDAG.R", chdir = TRUE)
source("utils/metric.R", chdir = TRUE)
source("utils/data_generator.R", chdir = TRUE)
source("OrderMCMC.R")
source("competing_models/minimap.R")


# We have 7 experiment types in this code. Each type of experiment would give us
# replication of the results (Figures and Tables) in the paper.

##############################################################################
### experiment_type == "among" (Table 1 and Table 2.)
# The experiment is to compare the performance of our method ("ORDER") with that
# of other competing methods ("TD", "LISTEN", "GDS").
# Output: log files of performance table.

# Example of running code.
# Rscript run.R -x "among" -m "ORDER" -t 'unif' -N 30 -n 1000 -p 40 -s 1 -e 0.1 -l 0.5

# Rscript run.R
# -x "among" : it indicates that the experiment is among competing models.
# -m c("ORDER","TD", "LISTEN","GDS") : choosing among 4 models
# -t c("unif", "mis", "gauss", "same") : choosing among 4 options for signals.
# -N 30  : Choosing the number of core for parallel computings. (int)
# -n 1000  : Choosing the number of samples (int)
# -p 40  : Choosing the number of nodes (int)
# -s 1  : Choosing the value of error variance (numeric)
# -e 0.1  : Choosing the value of edge probability (numeric in (0, 1))
# -l 0.3  : Choosing the value of lower bound of signal size (numeric)
# -i 0.3  : Choosing the value of degree of misspecification. (numeric in (0, 1); only when t == "mis")
##############################################################################


##############################################################################
### experiment_type == "mixing" (Figure 1.)
# The experiment is to draw a mixing diagram to show chains with three
# neighborhoods are reasonably mixing well in moderate number of nodes.
# Output: a mixing plot.


# Example of running code.
# Rscript run.R -x "mixing" -m "ORDER"


# Rscript run.R
# -x "mixing" : it indicates that the experiment is to draw a mixing diagram.
# -m c("ORDER","MINIMAP") : choosing between 2 models

##############################################################################

##############################################################################
# experiment_type == "misspecify" (Figure 2.)
# The experiment is to draw boxplots of hamming distance according to the degree
# of misspecification for four models; "ORDER","TD", "LISTEN", "MINIMAP"

# Example of running code.
# Rscript run.R -x "misspecify"
# Output: boxplots of hamming distance.
##############################################################################

##############################################################################
# experiment_type == "highdim" 
# The experiment is to show high dimensional experiments

# Example of running code.
# Rscript run.R -x "highdim" -k c(1,2,3,4,5,6,7)
# Output: log files of performance table.
##############################################################################

##############################################################################
# experiment_type == "bias" 
# The experiment is to measure bias of the assumption

# Example of running code.
# Rscript run.R -x "bias" -i c(0, 0.1, ..., 0.9)
# Output: values of performance table.
##############################################################################

##############################################################################
# experiment_type == "td" 
# The experiment is to compare performance between TD and ITD

# Example of running code.
# Rscript run.R -x "td"
# Output: boxplots of log posterior difference (ITD - TD).
##############################################################################

##############################################################################
# experiment_type == "increasing" 
# The experiment is to show weakly increasing error variance

# Example of running code.
# Rscript run.R -x "increasing" -i c(0, 0.1, ..., 0.9)
# Output: values in performance table.
##############################################################################



option_list = list(
  make_option(c("-x", "--experiment_type"), type="character", default=NULL,
              help="type of experiment: among, mixing, misspecify, highdim, bias, td, increasing",
              metavar="character"),
  # experiment_type == "among"
  make_option(c("-m", "--model"), type="character", default=NULL,
              help="model names",  metavar="character"),
  make_option(c("-t", "--type"), type="character", default=NULL,
              help="signal type",  metavar="character"),
  make_option(c("-N", "--n_rep"), type="integer", default=NULL,
              help="the number of repetition in simulation."),
  make_option(c("-n", "--n"), type="integer", default=NULL,
              help="the number of sample of data X in simulation."),
  make_option(c("-p", "--p"), type="integer", default=NULL,
              help="the number of variable of data X in simulation."),
  make_option(c("-s", "--sig"), type="numeric", default=1,
              help="the error variance in simulation."),
  make_option(c("-e", "--edge_prob"), type="numeric", default=NULL,
              help="probability of edge creating for each couple of nodes"),
  make_option(c("-l", "--lower_signal"), type="numeric", default=NULL,
              help="numeric in (0,1) to generate signals. +- Unif[l, 1]"),
  make_option(c("-i", "--eps"), type="numeric", default=0,
              help=" the value of degree of misspecification in (0,1). Unif[1-i, 1+i]"),
  make_option(c("-k", "--size"), type="integer", default=NULL,
              help="from 1 to 7 in increasing complexity.",
              metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.element(opt$experiment_type, c("among", "mixing", "misspecify", "highdim", "bias", "td", "increasing")) != TRUE){
  stop("Please select proper experiment type!", call.=FALSE)
}
# we fix the hyperparameters for "ORDER".
alpha = 0.99
gamma = 0.01
kappa = 0
c1 = 1
c2 = 3


set.seed(20220103)
if (opt$experiment_type == "among"){

  if (is.null(opt$type) | is.null(opt$n_rep) | is.null(opt$n) | is.null(opt$p)|
      is.null(opt$edge_prob) |is.null(opt$lower_signal) ){
    stop("Please fill out the options!", call.=FALSE)
  }
  if ((opt$edge_prob > 1) | (opt$edge_prob < 0) | (opt$lower_signal > 1) | (opt$lower_signal < 0)){
    stop("Please fill in right range of edge_prob or lower_signal!", call.=FALSE)
  }
  # data-generating parameters
  model = opt$model # c("ORDER", "TD", "LISTEN", "GDS_EEV", "ITD")
  n = opt$n # the number of samples
  p = opt$p # the number of covariates
  n_rep = opt$n_rep # the number of simulations (repetition)
  sig = opt$sig # the error variance in simulation.
  edge_prob = opt$edge_prob # edge probability in generating true DAGs.
  lower_signal = opt$lower_signal # minimum signal size for weighted adjacency matrix.
  type = opt$type # c("unif","mis","gauss"), data_generating process.
  eps = opt$eps # The degree of misspecification. Only when type == "mis"

  # model parameter for "ORDER". (See the paper.)
  c3 = (alpha*p*n + kappa)/2
  c4 = (c2 * log(p) + 0.5 * log(1 + alpha / gamma) + log(c1))
  
  # the number of MC samples.
  if (p >= 40){
    m = 3000
  }else{
    m = 1000
  }
  # burnin = m %/% 2
  burnin = 2500
  # dataset generating process
  if (type == 'gauss'){ # signal from gaussian.
    dat = data_generator_gaussian(n, p, n_rep, sig, edge_prob)
  }else if(type == 'mis'){ # error variances are not equal with degree == 'eps'.
    if ((eps <= 0) | (eps > 1) ){
      stop("Please select the value of degree of misspecification (-i) correctly.", call.=FALSE)
    }
    dat = data_generator_misspecified(n, p, n_rep, sig, edge_prob, lower_signal, eps)
  }else if(type == 'unif'){ # signal from uniform.
    dat = data_generator_unif(n, p, n_rep, sig, edge_prob, lower_signal)
  }else{
    stop("Please select proper signal type!(-t).", call.=FALSE)
  }

  # for parallel computing.
  cluster = makeCluster(n_rep)
  registerDoParallel(cluster)


  run = foreach(rep = 1 : n_rep,.packages = c("pcalg")) %dorng%
    {
      X = dat$X[[rep]]
      true_adj = dat$G[[rep]]

      # estimating process
      start = Sys.time()
      if(model == "ORDER"){
        V = crossprod(X)
        estimated_adj = order_mcmc(V, burnin, m, c3, c4)
      }else if(model == "LISTEN"){
        LISTEN = EqVarDAG_BU(X)
        estimated_adj = 1 * LISTEN$adj
      }else if(model == "TD"){
        TD = EqVarDAG_TD(X)
        estimated_adj = 1 * TD$adj
      }else if(model == "GDS_EEV"){
        GDS_EEV = GDS(X)
        estimated_adj = GDS_EEV$Adj
      }else if(model == "ITD"){
        V = crossprod(X)
        estimated_adj = top_down_iterative(V, c3, c4,p)$DAG
      }else{
        stop("insert proper model name!")
      }
      end = Sys.time()
      timeit <- difftime(end, start, units = "secs")[[1]]

      # performance measure
      if(model=='ORDER'){
        hd = hammingDistance(true_adj, estimated_adj$RB)
        recall = calculate_recall(true_adj, estimated_adj$RB)
        flip = calculate_flip(true_adj, estimated_adj$RB)
        fdr = calculate_fdr(true_adj, estimated_adj$RB)
        est_edge = sum(estimated_adj$RB)
        true_edge = sum(true_adj)
      }else{
        hd = hammingDistance(true_adj, estimated_adj)
        recall = calculate_recall(true_adj, estimated_adj)
        flip = calculate_flip(true_adj, estimated_adj)
        fdr = calculate_fdr(true_adj, estimated_adj)
        est_edge = sum(estimated_adj)
        true_edge = sum(true_adj)
      }

      return(c(hd, recall, flip, fdr,est_edge, true_edge, timeit))
      }


  mat <- matrix(0, n_rep, 7)
  for(i in 1:n_rep){
    mat[i, ] <-run[[i]]
  }

  hd = mat[, 1]
  recall = mat[, 2]
  flip = mat[, 3]
  fdr = mat[, 4]
  est_edge = mat[, 5]
  true_edge = mat[, 6]
  timeit = mat[, 7]

  model_type = paste(model,'_',type,sep = '')
  save_metric(model_type, n, p, n_rep, edge_prob,lower_signal)
  log_metric(model_type, n, p, n_rep, edge_prob,lower_signal) # Output: log_file.

}else if(opt$experiment_type == "mixing"){
  if (is.null(opt$model)){
    stop("Please fill out the model, either ORDER or MINIMAP", call.=FALSE)
  }
  model = opt$model
  # data_generating information.
  n = 1000
  p = 20
  n_rep = 30
  sig = 1
  #edge_prob =  3 / (p - 1) / 2
  edge_prob = 0.1
  lower_signal = 0.5
  eps = 0
  
  # model parameter for "ORDER". (See the paper.)
  c3 = (alpha*n*p + kappa)/2
  c4 = (c2 * log(p) + 0.5 * log(1 + alpha / gamma) + log(c1))
  c5 = (alpha*n + kappa)/2
  data = data_generator_misspecified(n, p, 1, sig, edge_prob, lower_signal, eps)
  true = data$G[[1]]
  X = data$X[[1]]
  V = crossprod(X)
  
  init_perms = matrix(0, nrow = n_rep, ncol = p)
  for(iter in 1:n_rep){
    init_perms[iter, ] = sample(1:p,p)
  }
  init_scores = numeric(n_rep)
  cluster = makeCluster(n_rep)
  registerDoParallel(cluster)
  
  if (model == "MINIMAP"){
    
    m = 10000
    true_est = CI_test(V, 1:p, n, lv= 0.0001)
    est_max_score = 0
    for(i in 1:p){
      est_max_score = est_max_score + score(V, i, which(true_est[, i] == 1), c5, c4)
    }
    
    true_score = 0
    for(i in 1:p){
      true_score = true_score + score(V, i, which(true[, i] == 1), c5, c4)
    }
    
    
    for(iter in 1:n_rep){
      TO_init = init_perms[iter, ]
      DAG_init = CI_test(V, TO_init, n, lv=0.0001)
      for(i in 1:p){
        init_scores[iter] = init_scores[iter] + score(V, i, which(DAG_init[, i] == 1), c5, c4)
      }
    }
    print(min(init_scores))
    dir.create(file.path("./Figure1"), showWarnings = FALSE)
    filename=paste("Figure1/",model,"_",p,"_",n,"_",m,"_",eps,".pdf", sep = '')
    pdf(file = filename, width = 5, height=6)
    
    run <- foreach(rep = 1 : n_rep,.packages = c("pcalg")) %dorng%
      {
        TO_init = init_perms[rep, ]
        estimate = minimap(X, c5, c4, m, 0, TO_init, lv=0.0001)
        return(estimate)
      }
    grid =1:m
    
    #par(mar=c(5,7,4,2)+.4)
    par(mar=c(4,4,1,2))

    plot(grid, (run[[1]]$traj[2:(m+1)] + init_scores[1])/10000, type='s', ylim=c(min(init_scores) * 1.01/10000, true_score*0.99/10000) , ylab = c(" "),
         xaxt = 'n', xlab=c("Number of iterations"), col='gray', cex.lab = 2, cex.axis =2, las = 1)
    
    #ylab <- as.list(expression(paste("Log-posterior \u00D7"~10^-4)))
    #mtext("(a)\n\n\n\n\n\n\n\n\n\n\n", side = 2, line = 5, cex = 2.5, las = 2)
    #ylab[[2]] = "  " 
    #mtext(do.call(expression, ylab),side=2, line = c(5,1), cex = 2)
    
    par(cex.axis=2)
    axis(side=1, at=c(0,5000,10000))
    
    
    
    for (ttt in 2:n_rep){
      lines(grid, (run[[ttt]]$traj[2:(m+1)] + init_scores[ttt])/10000, type='s', col='gray')
    }
    lines(grid, y=rep(true_score/10000, m), col='red')
    
    
    legend("bottomright",
           legend = c("True model"),
           lty = c(1),
           col = c(2),
           lwd = 1,
           cex = 2)
    dev.off() # Output: 3 mixing plots.
  }else if(model == "ORDER"){
    
    
    true_score = score_DAG(V, true, c3, c4)
    #n_rep number of the initial permutations
    min_sse = 1/diag(solve(V))
    max_sse = diag(V)
    
    for(iter in 1:n_rep){
      sup_parents = candidate_parents(V, init_perms[iter, ], min_sse, max_sse, c3, c4, ncol(V))
      dag = order_to_DAG(sup_parents, c3, c4)
      init_scores[iter] = score_DAG(V, dag$DAG, c3, c4)
    }
    
    dir.create(file.path("./Figure1"), showWarnings = FALSE)
    for(method in c('adj_trs', 'ran_trs', 'ran-to-ran')){
    #for(method in c('adj_trs')){
      cluster = makeCluster(n_rep)
      registerDoParallel(cluster)
      if(method == 'adj_trs'){
        N = 5000
        filename=paste("Figure1/",method,"_",p,"_",n,"_",N,"_",eps,".pdf", sep = '')
      }else{
        N = 1500
        filename=paste("Figure1/",method,"_",p,"_",n,"_",N,"_",eps,".pdf", sep = '')
      }
      pdf(file = filename, width = 5, height=6)
      
      
      run <- foreach(rep = 1 : n_rep,.packages = c("pcalg")) %dorng%
        {
          estimate = order_mcmc_mixing(V, N, method, init_perms[rep, ], c3, c4)
          return(estimate)
        }
      if(method == 'adj_trs'){
        grid = 2*(1:N)
      }else{
        grid = p*(1:N) / 3
      }
      
      #par(mar=c(5,7,4,2)+.4)
      par(mar=c(4,4,1,2))
      plot(grid, run[[1]]$traj/10000, type='s', ylim=c(min(init_scores) * 1.01/10000, true_score*0.99/10000) , ylab = c("  "),
           xaxt = 'n', xlab=c("Effective no. of iterations"), col='gray', cex.lab = 2, cex.axis =2, las=1)
      #ylab <- as.list(expression(paste("Log-posterior \u00D7"~10^-4)))
      #ylab <- as.list(expression(paste("Log-posterior \u00D7"~10^-3)))
      
      #ylab[[2]] = "  " 
      #mtext(do.call(expression, ylab),side=2, line = c(5,1), cex = 2 )
      #mtext("(d)\n\n\n\n\n\n\n\n\n\n\n", side = 2, line = 5, cex = 2.5, las = 2)
      
      par(cex.axis=2)
      axis(side=1, at=c(0,5000,10000))
      
      for (ttt in 2:n_rep){
        lines(grid, run[[ttt]]$traj/10000, type='s', col='gray')
      }
      lines(grid, y=rep(true_score/10000,N), col='red')
      
      legend("bottomright",
             legend = c("True model"),
             lty = c(1),
             col = c(2),
             lwd = 1,
             cex = 2)
      dev.off() # Output: 3 mixing plots.
    }
  }else{
    stop("Please fill out the model, either ORDER or MINIMAP", call.=FALSE)
  }
  
}else if(opt$experiment_type == "misspecify"){
  dir.create(file.path("./Figure2"), showWarnings = FALSE)
  dir.create(file.path("./Figure2/output_hd"), showWarnings = FALSE)

  # data_generating information.
  n = 500
  p = 40
  n_rep = 30
  sig = 1
  edge_prob =  3 / (p - 1) / 2
  lower_signal = 0.3

  # model parameter for "ORDER". (See the paper.)
  c3 = (alpha*p*n + kappa)/2
  c4 = (c2 * log(p) + 0.5 * log(1 + alpha / gamma) + log(c1))
  c5 = (alpha*n + kappa)/2
  m = 3000
  

  for (model in c("ORDER", "TD", "LISTEN", "MINIMAP")){
    for (eps in 0:9/10){
      dat = data_generator_misspecified(n, p, n_rep, sig, edge_prob, lower_signal, eps)
      cluster = makeCluster(n_rep)
      registerDoParallel(cluster)

      run = foreach(rep = 1 : n_rep,.packages = c("pcalg")) %dorng%
        {
          X = dat$X[[rep]]
          true_adj = dat$G[[rep]]
          if(model == "ORDER"){
            V = crossprod(X)
            estimated_adj = order_mcmc(V, m %/% 2, m, c3, c4)
          }else if(model == "LISTEN"){
            LISTEN = EqVarDAG_BU(X)
            estimated_adj = 1 * LISTEN$adj
          }else if(model == "TD"){
            TD = EqVarDAG_TD(X)
            estimated_adj = 1 * TD$adj
          }else if(model == "MINIMAP"){
            V = crossprod(X)
            TO_init = top_down_iterative(V, c3, c4, p)$TO
            estimated_adj = minimap(X, c5, c4, m, m %/% 2, TO_init, lv=0.0001)
          }

          # calculate hamming distance
          if(model=="ORDER"){
            hd = hammingDistance(true_adj, estimated_adj$RB)
          }else if(model == "MINIMAP"){
            hd = hammingDistance(true_adj, estimated_adj$Ave)
          }else{
            hd = hammingDistance(true_adj, estimated_adj)
          }
          return(hd)
        }

      stopCluster(cluster)
      
      hd <- numeric(n_rep)
      for(i in 1:n_rep){
        hd[i] <-run[[i]]
      }
      save_hd(model,lower_signal, eps)
      log_hd(model, lower_signal, eps)
    }
  }
  note = c(1)
  for (tt in c("0","0.1", "0.2", "0.3","0.4","0.5","0.6","0.7","0.8","0.9")){
    for (ss in c("ORDER", "TD", "LISTEN", "MINIMAP" )){
      load(paste("./Figure2/output_hd/",ss, "_0.3_",tt ,".Rdata", sep=""))
      #print(paste(ss, tt, mean(hd), sd(hd)))
      note = c(note, hd)
    }
  }
  HammingDistance = note[2:length(note)]
  # create a data frame
  Heterogeneity=rep(c("0","0.1", "0.2", "0.3","0.4","0.5","0.6","0.7","0.8","0.9"),
                       each=120)
  Method= rep(c("Proposed","TD","LISTEN", "MINIMAP"),each=30)
  data=data.frame(Heterogeneity, Method, HammingDistance)

  # grouped boxplot
  filename=paste("Figure2/misspecification.pdf", sep = '')
  pdf(file = filename, width = 10,height = 3)
  myplot = ggplot(data, aes(x=Heterogeneity, y=HammingDistance,
                            fill= factor(Method, levels = c("Proposed", "TD", "LISTEN", "MINIMAP")))) +
    geom_boxplot() + labs(fill = "Method")
  myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 legend.position = c(0.08, 0.8))
  dev.off() # Output: figure 2.
}else if(opt$experiment_type == "highdim"){

  model = "ORDER"
  k = opt$size
  m = 3000
  n = 30 * (k + 1)
  p = 7 * 2^(k-1)
  D = 0.2 * sqrt(n)
  edge_prob = D / (p-1)
  if(n < p){highdim = T}else{highdim = F}
  
  c3 = (alpha*p*n + kappa)/2
  c4 = (c2 * log(p) + 0.5 * log(1 + alpha / gamma) + log(c1))
  sig = 1
  lower_signal = 0.5
  n_rep = 30
  dat = data_generator_unif(n, p, n_rep, sig, edge_prob, lower_signal)
  
  cluster = makeCluster(n_rep)
  registerDoParallel(cluster)
  
  run = foreach(rep = 1 : n_rep,.packages = c("pcalg")) %dorng%
    {
      X = dat$X[[rep]]
      true_adj = dat$G[[rep]]
      
      # estimating process
      start = Sys.time()
      V = crossprod(X)
      if(highdim == T){
        estimated_adj = order_mcmc(V, m %/% 2, m, c3, c4, n%/%2, highdim = T)
      }else{
        estimated_adj = order_mcmc(V, m %/% 2, m, c3, c4)
      }
      end = Sys.time()
      timeit <- difftime(end, start, units = "secs")[[1]]
      
      hd = hammingDistance(true_adj, estimated_adj$RB)
      recall = calculate_recall(true_adj, estimated_adj$RB)
      flip = calculate_flip(true_adj, estimated_adj$RB)
      fdr = calculate_fdr(true_adj, estimated_adj$RB)
      est_edge = sum(estimated_adj$RB)
      true_edge = sum(true_adj)
      return(c(hd, recall, flip, fdr,est_edge, true_edge, timeit))
    }
  mat <- matrix(0, n_rep, 7)
  for(i in 1:n_rep){
    mat[i, ] <-run[[i]]
  }
  
  hd = mat[, 1]
  recall = mat[, 2]
  flip = mat[, 3]
  fdr = mat[, 4]
  est_edge = mat[, 5]
  true_edge = mat[, 6]
  timeit = mat[, 7]
  
  type = "hd_setting"
  model_type = paste(model,"_",type,sep = "")
  save_metric(model_type, n, p, n_rep, edge_prob,lower_signal)
  log_metric(model_type, n, p, n_rep, edge_prob,lower_signal) # Output: log_file.

}else if(opt$experiment_type == "bias"){
  p = 7
  mode = "unif" # either "unif" or "IG"
  eps = opt$eps
  deg = 3
  # data generator
  
  n = 100 * p
  n_rep = 30
  sig = 1
  edge_prob =  3 / (p - 1) / 2
  lower_signal = 0.3
  
  c3 = (alpha*n*p + kappa)/2
  c4 = (c2 * log(p) + 0.5 * log(1 + alpha / gamma) + log(c1))
  c5 = (alpha*p + kappa)/2
  if(mode == "unif"){
    dat = data_generator_misspecified(n, p, n_rep, sig, edge_prob, lower_signal, eps)
  }else if(mode == "IG"){
    dat = data_generator_invgamma(n, p, n_rep, sig, edge_prob, lower_signal, deg)
  }
  
  cluster = makeCluster(n_rep)
  registerDoParallel(cluster)
  perm = permutations(p,p)
  
  run = foreach(rep = 1 : n_rep,.packages = c("pcalg")) %dorng%
    {
      X = dat$X[[rep]]
      true_adj = dat$G[[rep]]
      V = crossprod(X)
      max_sse_pvector = diag(V)
      min_sse_pvector = 1/diag(solve(V))
      
      score_neq = numeric(nrow(perm))
      score_eq = numeric(nrow(perm))
      DAG = matrix(0, p, p)
      DAG_EQ = matrix(0, p, p)
      
      for(i in 1:nrow(perm)){
        if(i == 1){
          sup_parents = candidate_parents(V, perm[i, ], min_sse_pvector, max_sse_pvector, c3, c4, p)
          dag = order_to_DAG(sup_parents, c3, c4)
          init_score = dag$score
          
          dag_eq = CI_test(V, perm[i, ], n, lv = 0.0001) # init_G
          init_score_eq = 0
          for(j in 1:p){
            init_score_eq = init_score_eq + score(V, j, which(dag_eq[, j] == 1), c5, c4)
          }
        }
        
        
        
        sup_parents = candidate_parents(V, perm[i, ], min_sse_pvector, max_sse_pvector, c3, c4, p)
        dag = order_to_DAG(sup_parents, c3, c4)
        print(dag$score)
        score_neq[i] = exp(dag$score - init_score)
        DAG = DAG + score_neq[i] * dag$DAG
        
        dag_eq = CI_test(V, perm[i, ], n, lv = 0.05) # init_G
        for(j in 1:p){
          score_eq[i] = score_eq[i] + score(V, j, which(dag_eq[, j] == 1), c5, c4)
        }
        score_eq[i] = exp(score_eq[i] - init_score_eq)
        DAG_EQ = DAG_EQ + score_eq[i] * dag_eq
      }
      DAG_EQ = DAG_EQ / sum(score_eq)
      DAG = DAG / sum(score_neq)
      
      result = matrix(0, 2,7)
      
      hd = hammingDistance(true_adj, DAG)
      hd_undirected = hammingDistance(true_adj + t(true_adj), DAG + t(DAG))/2
      recall = calculate_recall(true_adj, DAG)
      flip = calculate_flip(true_adj, DAG)
      fdr = calculate_fdr(true_adj, DAG)
      est_edge = sum(DAG)
      true_edge = sum(true_adj)
      
      result[1, ] = c(hd, hd_undirected, recall, flip, fdr,est_edge, true_edge)
      
      hd = hammingDistance(true_adj, DAG_EQ)
      hd_undirected = hammingDistance(true_adj + t(true_adj), DAG_EQ + t(DAG_EQ))/2
      recall = calculate_recall(true_adj, DAG_EQ)
      flip = calculate_flip(true_adj, DAG_EQ)
      fdr = calculate_fdr(true_adj, DAG_EQ)
      est_edge = sum(DAG_EQ)
      
      result[2, ] = c(hd, hd_undirected, recall, flip, fdr,est_edge, true_edge)
      
      return(result)
    }
  
  
  mat <- matrix(0, n_rep, 7)
  mat_eq <- matrix(0, n_rep, 7)
  for(i in 1:n_rep){
    mat[i, ] <-run[[i]][1, ]
    mat_eq[i, ] <-run[[i]][2, ]
  }
  
  hd = mat[, 1]
  hd_undirected = mat[, 2]
  recall = mat[, 3]
  flip = mat[, 4]
  fdr = mat[, 5]
  est_edge = mat[, 6]
  true_edge = mat[, 7]
  
  cat(paste("The result from model: ours" , "\n", "n_sample:", n,
            "n_var:", p ,"n_rep:", n_rep,"eps:", eps, "\n"))
  cat(paste("hd: ",sprintf("%.1f", mean(hd)), "$pm$", sprintf("%.1f", sd(hd)/sqrt(n_rep)), "\n", sep = ''))
  cat(paste("hd_undirected: ",sprintf("%.1f", mean(hd_undirected)), "$pm$", sprintf("%.1f", sd(hd_undirected)/sqrt(n_rep)), "\n", sep = ''))
  #cat(paste("tau: ",sprintf("%.1f", mean(tau)), "$pm$", sprintf("%.1f", sd(tau), "\n","\n", sep = ''))
  cat(paste("recall: ",sprintf("%.1f", mean(recall)*100), "$pm$", sprintf("%.1f", sd(recall)/sqrt(n_rep)*100), "\n", sep = ''))
  cat(paste("flip: ",sprintf("%.1f", mean(flip)*100), "$pm$", sprintf("%.1f", sd(flip)/sqrt(n_rep)*100), "\n", sep = ''))
  cat(paste("fdr: ",sprintf("%.1f", mean(fdr)*100), "$pm$", sprintf("%.1f", sd(fdr)/sqrt(n_rep)*100), "\n", sep = ''))
  cat(paste("est_edge: ",sprintf("%.1f", mean(est_edge)), "$pm$", sprintf("%.1f", sd(est_edge)/sqrt(n_rep)), "\n", sep = ''))
  cat(paste("true_edge: ",sprintf("%.1f", mean(true_edge)), "$pm$", sprintf("%.1f", sd(true_edge)/sqrt(n_rep)), "\n", sep = ''))
  
  
  
  
  hd = mat_eq[, 1]
  hd_undirected = mat_eq[, 2]
  recall = mat_eq[, 3]
  flip = mat_eq[, 4]
  fdr = mat_eq[, 5]
  est_edge = mat_eq[, 6]
  true_edge = mat_eq[, 7]
  
  
  cat(paste("The result from model: minimap" , "\n", "n_sample:", n,
            "n_var:", p ,"n_rep:", n_rep, "eps:", eps, "\n"))
  cat(paste("hd: ",sprintf("%.1f", mean(hd)), "$pm$", sprintf("%.1f", sd(hd)/sqrt(n_rep)), "\n", sep = ''))
  cat(paste("hd_undirected: ",sprintf("%.1f", mean(hd_undirected)), "$pm$", sprintf("%.1f", sd(hd_undirected)/sqrt(n_rep)), "\n", sep = ''))
  #cat(paste("tau: ",sprintf("%.1f", mean(tau)), "$pm$", sprintf("%.1f", sd(tau), "\n","\n", sep = ''))
  cat(paste("recall: ",sprintf("%.1f", mean(recall)*100), "$pm$", sprintf("%.1f", sd(recall)/sqrt(n_rep)*100), "\n", sep = ''))
  cat(paste("flip: ",sprintf("%.1f", mean(flip)*100), "$pm$", sprintf("%.1f", sd(flip)/sqrt(n_rep)*100), "\n", sep = ''))
  cat(paste("fdr: ",sprintf("%.1f", mean(fdr)*100), "$pm$", sprintf("%.1f", sd(fdr)/sqrt(n_rep)*100), "\n", sep = ''))
  cat(paste("est_edge: ",sprintf("%.1f", mean(est_edge)), "$pm$", sprintf("%.1f", sd(est_edge)/sqrt(n_rep)), "\n", sep = ''))
  cat(paste("true_edge: ",sprintf("%.1f", mean(true_edge)), "$pm$", sprintf("%.1f", sd(true_edge)/sqrt(n_rep)), "\n", sep = ''))
  
  
  
}else if(opt$experiment_type == "td"){
 
  diff = matrix(0,5,30)
  for (k in 1:5){
    n = 30 * (k + 1) # the number of samples
    p = 7 * 2^(k-1) # the number of covariates
    n_rep = 30 # the number of simulations (repetition)
    sig = 1 # the error variance in simulation.
    # edge_prob = 3/2/(p-1) # edge probability in generating true DAGs.
    edge_prob = 0.2 * sqrt(n)/(p-1) # edge probability in generating true DAGs.
    lower_signal = 0.1 # minimum signal size for weighted adjacency matrix.
    
    c3 = (alpha*p*n + kappa)/2
    c4 = (c2 * log(p) + 0.5 * log(1 + alpha / gamma) + log(c1))
    m = 3000
    
    #dat = data_generator_unif(n, p, n_rep, sig, edge_prob, lower_signal)
    dat = data_generator_gaussian(n, p, n_rep, sig, edge_prob)
    
    cluster = makeCluster(n_rep)
    registerDoParallel(cluster)
    
    
    run = foreach(rep = 1 : n_rep,.packages = c("pcalg")) %dorng%
      {
        X = dat$X[[rep]]
        true_adj = dat$G[[rep]]
        V = crossprod(X)
        min_sse_pvector = 1/diag(solve(V))
        max_sse_pvector = diag(V)
        max_num_par = p
        
        TD = EqVarDAG_TD(X)
        ITD = top_down_iterative(V, c3, c4, p)
        
        td_ordering = TD$TO
        itd_ordering = ITD$TO
        
        true_parent = candidate_parents(V, 1:p , min_sse_pvector, max_sse_pvector, c3, c4, max_num_par)
        TD_parent = candidate_parents(V, td_ordering, min_sse_pvector, max_sse_pvector, c3, c4, max_num_par)
        ITD_parent = candidate_parents(V, itd_ordering, min_sse_pvector, max_sse_pvector, c3, c4, max_num_par)
        
        TD_score = order_to_DAG(TD_parent, c3, c4)$score
        ITD_score = order_to_DAG(ITD_parent, c3, c4)$score
        true_score = order_to_DAG(true_parent, c3, c4)$score
        
        
        return(c(TD_score, ITD_score, true_score))
      }
    
    mat <- matrix(0, n_rep, 3)
    for(i in 1:n_rep){
      mat[i, ] <-run[[i]]
    }
    
    td_score = mat[, 1]
    itd_score = mat[, 2]
    true_score =mat[, 3]
    diff[k, ] = itd_score - td_score
  }
  
  dir.create(file.path("./TDvsITD"), showWarnings = FALSE)
  save(diff, file = "./TDvsITD/diff.Rdata")
  
  load("./TDvsITD/diff.Rdata")
  mg = max(diff, -diff)
  d = data.frame(
    logdiff = c(diff),
    p = as.factor(c(rep(7,30), rep(14,30),rep(28,30),rep(56,30),rep(112,30)))
  )
  means <- aggregate(logdiff ~  p, d, mean)
  means$logdiff <-  round(means$logdiff, 3)
  
  filename = "./TDvsITD/diff.pdf" 
  pdf(file = filename, width = 10, height=3)
  ggplot(d, aes(x=p, y=logdiff, fill = p)) +
    geom_boxplot()+ labs(x="p", y = "log posterior difference") +
    coord_cartesian(ylim = c(-50, mg)) +
    stat_summary(fun=mean, geom="point", shape=20, size=3, color="red", fill="red") +
    geom_text(data = means, aes(label = logdiff, y = logdiff, x = c(1.25,2.25,3.25,4.25,5.25)), size =5, color='red') +
    theme(axis.text=element_text(size=10)) +
    theme_classic()
  dev.off()
}else if(opt$experiment_type == "increasing"){
  
  
  # data_generating information.
  n = 1000
  p = 40
  n_rep = 30
  sig = 1
  edge_prob =  3 / (p - 1) / 2
  lower_signal = 0.3
  eps = opt$eps
  type = "unif" # or "unif"
  model = "ORDER"
  
  # model parameter for "ORDER". (See the paper.)
  alpha = 0.99
  gamma = 0.01
  kappa = 0
  c1 = 1
  c2 = 3
  c3 = (alpha*p*n + kappa)/2
  c4 = (c2 * log(p) + 0.5 * log(1 + alpha / gamma) + log(c1))
  
  if (p >= 40){
    m = 3000
  }else{
    m = 1000
  }
  
  
  dat = data_generator_increasing(n, p, n_rep, sig, edge_prob, lower_signal, eps, type)
  cluster = makeCluster(n_rep)
  registerDoParallel(cluster)
  
  
  run = foreach(rep = 1 : n_rep,.packages = c("pcalg")) %dorng%
    {
      X = dat$X[[rep]]
      true_adj = dat$G[[rep]]
      
      start = Sys.time()
      V = crossprod(X)
      estimated_adj = order_mcmc(V, m %/% 2, m, c3, c4)
      timeit = Sys.time() - start
      # calculate hamming distance
      
      hd = hammingDistance(true_adj, estimated_adj$RB)
      recall = calculate_recall(true_adj, estimated_adj$RB)
      flip = calculate_flip(true_adj, estimated_adj$RB)
      fdr = calculate_fdr(true_adj, estimated_adj$RB)
      est_edge = sum(estimated_adj$RB)
      true_edge = sum(true_adj)
      
      
      
      return(c(hd, recall, flip, fdr,est_edge, true_edge, timeit))
    }
  
  
  mat <- matrix(0, n_rep, 7)
  for(i in 1:n_rep){
    mat[i, ] <-run[[i]]
  }
  
  hd = mat[, 1]
  recall = mat[, 2]
  flip = mat[, 3]
  fdr = mat[, 4]
  est_edge = mat[, 5]
  true_edge = mat[, 6]
  timeit = mat[, 7]
  
  model_type = paste(model,"_",type,"_", eps,sep ="")
  save_metric(model_type, n, p, n_rep, edge_prob,lower_signal)
  log_metric(model_type, n, p, n_rep, edge_prob,lower_signal) # Output: log_file.
  
}

