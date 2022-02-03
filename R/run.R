rm(list = ls())
library(foreach, quietly = TRUE)
library(doParallel, quietly = TRUE)
library(doRNG, quietly = TRUE)
library(iterators, quietly = TRUE)
library(pcalg, quietly = TRUE)
library(optparse, quietly = TRUE)
library(ggplot2)
library(mvtnorm)
source("competing_models/GDS/startups/startupGDS.R", chdir = TRUE)
source("competing_models/EqVarDAG/startupEqVarDAG.R", chdir = TRUE)
source("utils/metric.R", chdir = TRUE)
source("utils/data_generator.R", chdir = TRUE)
source("OrderMCMC.R")
source("competing_models/minimap.R")


# We have 5 experiment types in this code. Each type of experiment would give us
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
# Rscript run.R -x "mixing"

# Rscript run.R
# -x "mixing" : it indicates that the experiment is to draw a mixing diagram.
##############################################################################

##############################################################################
# experiment_type == "misspecify" (Figure 2.)
# The experiment is to draw boxplots of hamming distance according to the degree
# of misspecification for three models; "ORDER","TD", "LISTEN"

# Example of running code.
# Rscript run.R -x "misspecify"
# Output: boxplots of hamming distance.
##############################################################################

##############################################################################
# experiment_type == "within" (Figure 4.)
# The experiment is to compare the performance within our method ("ORDER").
# ("RB", "unweighted","ITD")
# Output: log files of performance table.

# Example of running code.
# Rscript run.R -x "within"
# Output: boxplots of hamming distance.
##############################################################################

option_list = list(
  make_option(c("-x", "--experiment_type"), type="character", default=NULL,
              help="type of experiment: among, mixing, misspecify, real, within",
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
              help=" the value of degree of misspecification in (0,1). Unif[1-i, 1+i]")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.element(opt$experiment_type, c("among", "mixing", "misspecify", "within", "minimap")) != TRUE){
  stop("Please select proper experiment type!", call.=FALSE)
}

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
  alpha = 0.99
  gamma = 0.01
  kappa = 0
  c1 = 1
  c2 = 3
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
        estimated_adj = top_down_iterative(V, c3, c4)$DAG
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

}else if(opt$experiment_type == "minimap"){
  

  n = 1000
  p = 40
  n_rep = 30
  sig = 1
  edge_prob =  0.1
  lower_signal = 0.3
  m =  3000
  
  # model parameter for "ORDER". (See the paper.)
  alpha = 0.99
  gamma = 0.01
  kappa = 0
  
  c1 = 1
  c2 = 3
  c3 = (alpha*n*p + kappa)/2
  c4 = (c2 * log(p) + 0.5 * log(1 + alpha / gamma) + log(c1))
  c5 = (alpha*p + kappa)/2
  for (model in c("ORDER", "MINIMAP")){
    for (eps in 0:9/10){
      dat = data_generator_misspecified(n, p, n_rep, sig, edge_prob, lower_signal, eps)
      cluster = makeCluster(n_rep)
      registerDoParallel(cluster)
      
      run = foreach(rep = 1 : n_rep,.packages = c("pcalg")) %dorng%
        {
          X = dat$X[[rep]]
          true_adj = dat$G[[rep]]
          V = crossprod(X)
          if(model == "ORDER"){
            estimated_adj = order_mcmc(V, m %/% 2, m, c3, c4)
          }else{
            TO_init = top_down_iterative(V, c3, c4)$TO
            estimated_adj = minimap(X, c5, c4, m, m %/% 2, TO_init, lv=0.0001)
          }# calculate hamming distances
          if(model=='ORDER'){
            hd_directed = calculate_L1_err(true_adj, estimated_adj$RB, 1)
            hd_undirected = calculate_L1_err(true_adj, estimated_adj$RB, 0)
            fliped = calculate_flip(true_adj, estimated_adj$RB)
          }else{
            hd_directed = calculate_L1_err(true_adj, estimated_adj$Ave, 1)
            hd_undirected = calculate_L1_err(true_adj, estimated_adj$Ave, 0)
            fliped = calculate_flip(true_adj, estimated_adj$Ave)
          }
          
          return(c(hd_directed, hd_undirected, fliped))
        }
      
      
      mat <- matrix(0, n_rep, 3)
      for(i in 1:n_rep){
        mat[i, ] <-run[[i]]
      }
      
      hd_directed = mat[, 1]
      hd_undirected = mat[, 2]
      fliped = mat[, 3]
      save_minimap(model, eps)
      log_minimap(model, eps)
    }
  }

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
  alpha = 0.99
  gamma = 0.01
  kappa = 0
  c1 = 1
  c2 = 3
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
    pdf(file = filename)
    
    run <- foreach(rep = 1 : n_rep,.packages = c("pcalg")) %dorng%
      {
        TO_init = init_perms[rep, ]
        estimate = minimap(X, c5, c4, m, 0, TO_init, lv=0.0001)
        return(estimate)
      }
    grid =1:m
    #plot(grid, run[[1]]$traj[2:(m+1)] + init_scores[1], type='s', ylim=c(min(init_scores) * 1.03, est_max_score*0.97) , ylab = c("log-posterior"),
    #     xlab=c("iterations"), col='gray', cex.lab = 2, cex.axis =2)
    par(mar=c(5,7,4,2)+.4)
    
    plot(grid, (run[[1]]$traj[2:(m+1)] + init_scores[1])/10000, type='s', ylim=c(min(init_scores) * 1.01/10000, true_score*0.99/10000) , ylab = c(" "),
         xlab=c("The number of iterations"), col='gray', cex.lab = 2, cex.axis =2, las = 1)
    ylab <- as.list(expression(paste("Log-posterior \u00D7"~10^-4)))
    #mtext("(a)\n\n\n\n\n\n\n\n\n\n\n", side = 2, line = 5, cex = 2.5, las = 2)
    ylab[[2]] = "  " 
    mtext(do.call(expression, ylab),side=2, line = c(5,1), cex = 2)
    
    
    
    
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
      sup_parents = candidate_parents(V, init_perms[iter, ], min_sse, max_sse, c3, c4)
      dag = order_to_DAG(sup_parents, c3, c4)
      init_scores[iter] = score_DAG(V, dag$DAG, c3, c4)
    }
    
    dir.create(file.path("./Figure1"), showWarnings = FALSE)
    #for(method in c('adj_trs', 'ran_trs', 'ran-to-ran')){
    for(method in c('adj_trs')){
      cluster = makeCluster(n_rep)
      registerDoParallel(cluster)
      if(method == 'adj_trs'){
        N = 5000
        filename=paste("Figure1/",method,"_",p,"_",n,"_",N,"_",eps,".pdf", sep = '')
        pdf(file = filename)
      }else{
        N = 1500
        filename=paste("Figure1/",method,"_",p,"_",n,"_",N,"_",eps,".pdf", sep = '')
        pdf(file = filename)
      }
      
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
      
      par(mar=c(5,7,4,2)+.4)
      plot(grid, run[[1]]$traj/10000, type='s', ylim=c(min(init_scores) * 1.01/10000, true_score*0.99/10000) , ylab = c("  "),
           xlab=c("The effective number of iterations"), col='gray', cex.lab = 2, cex.axis =2, las=1)
      ylab <- as.list(expression(paste("Log-posterior \u00D7"~10^-4)))
      #ylab <- as.list(expression(paste("Log-posterior \u00D7"~10^-3)))
      
      ylab[[2]] = "  " 
      mtext(do.call(expression, ylab),side=2, line = c(5,1), cex = 2 )
      #mtext("(d)\n\n\n\n\n\n\n\n\n\n\n", side = 2, line = 5, cex = 2.5, las = 2)
      
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

  for (model in c("ORDER", "TD", "LISTEN")){
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
          }

          # calculate hamming distance
          if(model=='ORDER'){
            hd = hammingDistance(true_adj, estimated_adj$RB)
          }else{
            hd = hammingDistance(true_adj, estimated_adj)
          }
          return(hd)
        }


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
    for (ss in c("ORDER", "TD", "LISTEN" )){
      load(paste("./Figure2/output_hd/",ss, "_0.3_",tt ,".Rdata", sep=""))
      #print(paste(ss, tt, mean(hd), sd(hd)))
      note = c(note, hd)
    }
  }
  HammingDistance = note[2:length(note)]
  # create a data frame
  Heterogeneity=rep(c("0","0.1", "0.2", "0.3","0.4","0.5","0.6","0.7","0.8","0.9"),
                       each=90)
  Method= rep(c("Proposed","TD","LISTEN"),each=30)
  data=data.frame(Heterogeneity, Method, HammingDistance)

  # grouped boxplot
  filename=paste("Figure2/misspecification.pdf", sep = '')
  pdf(file = filename, width = 10,height = 3)
  myplot = ggplot(data, aes(x=Heterogeneity, y=HammingDistance,
                            fill= factor(Method, levels = c("Proposed", "TD", "LISTEN")))) +
    geom_boxplot() + labs(fill = "Method")
  myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 legend.position = c(0.08, 0.8))
  dev.off() # Output: figure 2.
}else if(opt$experiment_type == "within"){
  dir.create(file.path("./Figure4"), showWarnings = FALSE)

  # data_generating information.
  n = 500
  p = 40
  n_rep = 30
  sig = 1
  edge_prob = 0.1
  signal = 0.3

  # model parameter for "ORDER". (See the paper.)
  alpha = 0.99
  gamma = 0.01
  kappa = 0
  c1 = 1
  c2 = 3
  c3 = (alpha*p*n + kappa)/2
  c4 = (c2 * log(p) + 0.5 * log(1 + alpha / gamma) + log(c1))

  m = 10000

  dat = data_generator_unif(n, p, n_rep, sig, edge_prob, signal)

  cluster = makeCluster(n_rep)
  registerDoParallel(cluster)
  run = foreach(rep = 1 : n_rep,.packages = c("pcalg")) %dorng%
    {
      X = dat$X[[rep]]
      true_adj = dat$G[[rep]]
      # estimating process

      V = crossprod(X)
      out =  order_mcmc(V, m %/% 2, m, c3, c4)
      estimated_adj = out$RB
      unw_adj = out$Ave
      itd_adj = out$ITD

      L1_error_rb_0 = calculate_L1_err(estimated_adj, true_adj, 0)
      L1_error_rb_1 = calculate_L1_err(estimated_adj, true_adj, 1)
      L2_error_rb_0 = calculate_L2_err(estimated_adj, true_adj, 0)
      L2_error_rb_1 = calculate_L2_err(estimated_adj, true_adj, 1)
      L1_error_unw_0 = calculate_L1_err(unw_adj, true_adj, 0)
      L1_error_unw_1 = calculate_L1_err(unw_adj, true_adj, 1)
      L2_error_unw_0 = calculate_L2_err(unw_adj, true_adj, 0)
      L2_error_unw_1 = calculate_L2_err(unw_adj, true_adj, 1)
      L1_error_itd_0 = calculate_L1_err(itd_adj, true_adj, 0)
      L1_error_itd_1 = calculate_L1_err(itd_adj, true_adj, 1)
      L2_error_itd_0 = calculate_L2_err(itd_adj, true_adj, 0)
      L2_error_itd_1 = calculate_L2_err(itd_adj, true_adj, 1)
      return(c(L1_error_rb_0, L1_error_rb_1, L2_error_rb_0, L2_error_rb_1,
               L1_error_unw_0, L1_error_unw_1, L2_error_unw_0, L2_error_unw_1,
               L1_error_itd_0, L1_error_itd_1, L2_error_itd_0, L2_error_itd_1))
    }

  mat = matrix(0, 12, n_rep)
  for (i in 1:n_rep){
    mat[ , i] = run[[i]]
  }
  L1_error_rb_0 = mat[1, ]
  L1_error_rb_1 =  mat[2, ]
  L2_error_rb_0 =  mat[3, ]
  L2_error_rb_1 =  mat[4, ]
  L1_error_unw_0 =  mat[5, ]
  L1_error_unw_1 =  mat[6, ]
  L2_error_unw_0 =  mat[7, ]
  L2_error_unw_1 =  mat[8, ]
  L1_error_itd_0 =  mat[9, ]
  L1_error_itd_1 =  mat[10, ]
  L2_error_itd_0 =  mat[11, ]
  L2_error_itd_1 =  mat[12, ]

  line=paste("L1_error_rb_0: ", mean(L1_error_rb_0),"+-",sd(L1_error_rb_0)/sqrt(n_rep), "\n",
             "L1_error_unw_0: ", mean(L1_error_unw_0),"+-", sd(L1_error_unw_0)/sqrt(n_rep), "\n",
             "L1_error_itd_0: ", mean(L1_error_itd_0), "+-",sd(L1_error_itd_0)/sqrt(n_rep), "\n",
             "L1_error_rb_1: ", mean(L1_error_rb_1), "+-",sd(L1_error_rb_1)/sqrt(n_rep), "\n",
             "L1_error_unw_1: ", mean(L1_error_unw_1), "+-",sd(L1_error_unw_1)/sqrt(n_rep), "\n",
             "L1_error_itd_1: ", mean(L1_error_itd_1), "+-",sd(L1_error_itd_1)/sqrt(n_rep), "\n",
             "L2_error_rb_0: ", mean(L2_error_rb_0), "+-",sd(L2_error_rb_0)/sqrt(n_rep), "\n",
             "L2_error_unw_0: ", mean(L2_error_unw_0), "+-",sd(L2_error_unw_0)/sqrt(n_rep), "\n",
             "L2_error_itd_0: ", mean(L2_error_itd_0), "+-",sd(L2_error_itd_0)/sqrt(n_rep), "\n",
             "L2_error_rb_1: ", mean(L2_error_rb_1), "+-",sd(L2_error_rb_1)/sqrt(n_rep), "\n",
             "L2_error_unw_1: ", mean(L2_error_unw_1), "+-",sd(L2_error_unw_1)/sqrt(n_rep), "\n",
             "L2_error_itd_1: ", mean(L2_error_itd_1), "+-",sd(L2_error_itd_1)/sqrt(n_rep), "\n", sep = '')
  write(line,file=paste("./Figure4/log.txt",sep=''),append=TRUE)


  dt = list()
  dt$L1 <- data.frame(
    types = c(rep("RB" , n_rep) , rep("unweighted", n_rep), rep("ITD", n_rep)),
    Directed = c(L1_error_rb_1, L1_error_unw_1, L1_error_itd_1),
    Undirected = c(L1_error_rb_0, L1_error_unw_0, L1_error_itd_0)
  )
  dt$L2 <- data.frame(
    types = c(rep("RB" , n_rep) , rep("unweighted", n_rep), rep("ITD", n_rep)),
    Directed = c(L2_error_rb_1, L2_error_unw_1, L2_error_itd_1),
    Undirected = c(L2_error_rb_0, L2_error_unw_0, L2_error_itd_0)
  )
  
  i = 1
  for (d in dt){
    
    bilan <- aggregate(cbind(Directed, Undirected) ~ types , data=d , mean)
    rownames(bilan) <- bilan[,1]
    bilan <- as.matrix(bilan[,-1])
    bilan <- bilan[c(2,3,1),]
    #Plot boundaries
    lim <- 1.2*max(bilan)
    
    #A function to add arrows on the chart
    error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
      arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
    }
    
    #Then I calculate the standard deviation for each specie and condition :
    stdev <- aggregate(cbind(Directed, Undirected) ~ types , data=d , sd)
    rownames(stdev) <- stdev[,1]
    stdev <- as.matrix(stdev[,-1]) / sqrt(n_rep)
    stdev <- stdev[c(2,3,1),]
    #I am ready to add the error bar on the plot using my "error bar" function !
    
    filename=paste("Figure4/comparison_L",i,".pdf", sep = '')
    pdf(file = filename)
    ze_barplot <- barplot(bilan , beside=T , col=c("blue" , "skyblue", "red") ,
                          ylim=c(0,lim) , ylab=paste("L", i, " ", "Loss", sep=''), cex.axis = 2, cex.names = 2, cex.lab =1.5)
    error.bar(ze_barplot,bilan, stdev)
    legend("topright",
           legend = c("ITD", "unweighted", "RB"),
           fill = c("blue","skyblue","red"),
           cex = 2)
    dev.off()
    i = i + 1
  }
}





'
dt = list()
  dt$L1 <- data.frame(
    types = c(rep("RB" , n_rep) , rep("ITD", n_rep)),
    Directed = c(L1_error_rb_1,  L1_error_itd_1),
    Undirected = c(L1_error_rb_0, L1_error_itd_0)
  )
  dt$L2 <- data.frame(
    types = c(rep("RB" , n_rep) , rep("ITD", n_rep)),
    Directed = c(L2_error_rb_1, L2_error_itd_1),
    Undirected = c(L2_error_rb_0, L2_error_itd_0)
  )
  
  i = 1
  for (d in dt){
    
    bilan <- aggregate(cbind(Directed, Undirected) ~ types , data=d , mean)
    rownames(bilan) <- bilan[,1]
    bilan <- as.matrix(bilan[,-1])
    #bilan <- bilan[c(2,3,1),]
    #Plot boundaries
    lim <- 1.2*max(bilan)
    
    #A function to add arrows on the chart
    error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
      arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
    }
    
    #Then I calculate the standard deviation for each specie and condition :
    stdev <- aggregate(cbind(Directed, Undirected) ~ types , data=d , sd)
    rownames(stdev) <- stdev[,1]
    stdev <- as.matrix(stdev[,-1]) / sqrt(n_rep)
    #stdev <- stdev[c(2,3,1),]
    #I am ready to add the error bar on the plot using my "error bar" function !
    
    filename=paste("Figure4/comparison_L",i,".pdf", sep = '')
    pdf(file = filename)
    ze_barplot <- barplot(bilan , beside=T , col=c("blue" , "red") ,
                          ylim=c(0,lim) , ylab=paste("L", i, " ", "Loss", sep="""), cex.axis = 2, cex.names = 2, cex.lab =1.5)
    error.bar(ze_barplot,bilan, stdev)
    legend("topright",
           legend = c("ITD", "RB"),
           fill = c("blue","red"),
           cex = 2)
    dev.off()
    i = i + 1
  }
'
