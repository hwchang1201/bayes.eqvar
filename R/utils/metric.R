

hammingDistance <- function(G1,G2)
  # hammingDistance(G1,G2)
  #
  # Computes Hamming Distance between DAGs G1 and G2 with SHD(->,<-) = 1!!!!
  #
  # INPUT:  G1, G2     adjacency graph containing only zeros and ones: (i,j)=1 means edge from X_i to X_j.
  #
  # OUTPUT: hammingDis Hamming Distance between G1 and G2
  #
  # Copyright (c) 2012-2013  Jonas Peters [peters@stat.math.ethz.ch]
  # All rights reserved.  See the file COPYING for license terms.
{
  allMistakesOne <- FALSE
  if(allMistakesOne)
  {
    Gtmp <- (G1+G2)%%2
    Gtmp <- Gtmp + t(Gtmp)
    nrReversals <- sum(Gtmp == 2)/2
    nrInclDel <- sum(Gtmp == 1)/2
    hammingDis <- nrReversals + nrInclDel
  } else
  {
    hammingDis <- sum(abs(G1 - G2))
    # correction: dist(-,.) = 1, not 2
    hammingDis <- hammingDis - 0.5*sum(G1 * t(G1) * (1-G2) * t(1-G2) + G2 * t(G2) * (1-G1) * t(1-G1))
  }
  return(hammingDis)
}


calculate_recall <- function(true_adj, estimated_adj){
  # INPUT:  true_adj, estimated_adj
  # adjacency graph containing only zeros and ones: (i,j)=1 means edge from X_i to X_j.
  #
  # OUTPUT: recall(= TP/(TN + FN)) of estimated_adj
  # If all true edges are included in estimated edges, recall is 1.
  # (No matter how dense the estimated edges are.)
  ifelse(sum(true_adj),sum(true_adj*estimated_adj)/sum(true_adj),0)
}

calculate_flip <- function(true_adj, estimated_adj){
  # INPUT:  true_adj, estimated_adj
  # adjacency graph containing only zeros and ones: (i,j)=1 means edge from X_i to X_j.
  #
  # OUTPUT: the number of cases if true_adj[i,j] == 1 and estimated_adj[j,i] == 1
  #
  #
  ifelse(sum(true_adj),sum(true_adj*t(estimated_adj))/sum(true_adj),0)
}

calculate_fdr <- function(true_adj, estimated_adj){
  # INPUT:  true_adj, estimated_adj
  # adjacency graph containing only zeros and ones: (i,j)=1 means edge from X_i to X_j.
  #
  # OUTPUT: false discovery rate(= FP/(FP + TP)) of estimated_adj
  # If the estimated_adj is 0, then fdr = 0. Count ratio of false alarm.
  #
  ifelse(sum(estimated_adj),1-sum(true_adj*estimated_adj)/sum(estimated_adj),0)
}


calculate_L1_err <- function(G1, G2, directed=1){
  if (directed == 1){
    as.numeric(sum(abs(G1 - G2)))
  }else{
    as.numeric(sum(abs(G1 + t(G1) - G2 - t(G2)))/2)
  }
}

calculate_L2_err <- function(G1, G2, directed=1){
  if (directed == 1){
    as.numeric(sqrt(sum((G1 - G2)^2)))
  }else{
    as.numeric(sqrt(sum( (G1 + t(G1) - G2 - t(G2))^2 )/2 ))
  }
}





save_metric <- function(model, n_sample, n_var, n_rep, edge_prob, lower_signal){
  dir.create(file.path("./output"), showWarnings = FALSE)
  dir.create(file.path(paste("./output/",model,sep = '')), showWarnings = FALSE)
  save(list = c("hd", "recall", "flip", "fdr", "est_edge", "true_edge", "timeit"),
       file = paste("./output/",model,"/","metric_",n_sample,"_", n_var, "_",n_rep,"_",edge_prob, "_",lower_signal,
                    ".Rdata", sep = ''))
}

log_metric <- function(model, n_sample, n_var, n_rep, edge_prob, lower_signal){
  # load metrics
  dir.create(file.path("./output"), showWarnings = FALSE)
  dir.create(file.path("./output/log"), showWarnings = FALSE)
  load(paste("./output/",model,"/","metric_",n_sample,"_",n_var, "_",n_rep,
             "_",edge_prob, "_",lower_signal,".Rdata", sep = ''))
  line=paste("The result from model: ", model," // " , " n_sample: ", n_sample,
               " n_var: ", n_var ," n_rep: ", n_rep, "edge_prob: ",edge_prob, "lower_signal: ",lower_signal,"\n",
         "hd: ",sprintf("%.1f", mean(hd)), "$pm$", sprintf("%.1f", sd(hd)/sqrt(n_rep)), "\n",
         # "tau: ",sprintf("%.1f", mean(tau)), "$pm$", sprintf("%.1f", sd(tau), "\n",
         "recall: ",sprintf("%.1f", mean(recall)*100), "$pm$", sprintf("%.1f", sd(recall)/sqrt(n_rep)*100), "\n",
         "flip: ",sprintf("%.1f", mean(flip)*100), "$pm$", sprintf("%.1f", sd(flip)/sqrt(n_rep)*100), "\n",
         "fdr: ",sprintf("%.1f", mean(fdr)*100), "$pm$", sprintf("%.1f", sd(fdr)/sqrt(n_rep)*100), "\n",
         "est_edge: ",sprintf("%.1f", mean(est_edge)), "$pm$", sprintf("%.1f", sd(est_edge)/sqrt(n_rep)), "\n",
         "true_edge: ",sprintf("%.1f", mean(true_edge)), "$pm$", sprintf("%.1f", sd(true_edge)/sqrt(n_rep)), "\n",
         "timeit: ",sprintf("%.1f", mean(timeit)), "$pm$", sprintf("%.1f", sd(timeit)/sqrt(n_rep)), "\n","\n",sep='')
  write(line,file=paste("./output/log/", model,"_",lower_signal,".txt",sep=''),append=TRUE)
  cat(paste("The result from model:", model,"\n", "n_sample:", n_sample,
              "n_var:", n_var ,"n_rep:", n_rep, "\n"))
  cat(paste("hd: ",sprintf("%.1f", mean(hd)), "$pm$", sprintf("%.1f", sd(hd)/sqrt(n_rep)), "\n", sep = ''))
  #cat(paste("tau: ",sprintf("%.1f", mean(tau)), "$pm$", sprintf("%.1f", sd(tau), "\n","\n", sep = ''))
  cat(paste("recall: ",sprintf("%.1f", mean(recall)*100), "$pm$", sprintf("%.1f", sd(recall)/sqrt(n_rep)*100), "\n", sep = ''))
  cat(paste("flip: ",sprintf("%.1f", mean(flip)*100), "$pm$", sprintf("%.1f", sd(flip)/sqrt(n_rep)*100), "\n", sep = ''))
  cat(paste("fdr: ",sprintf("%.1f", mean(fdr)*100), "$pm$", sprintf("%.1f", sd(fdr)/sqrt(n_rep)*100), "\n", sep = ''))
  cat(paste("est_edge: ",sprintf("%.1f", mean(est_edge)), "$pm$", sprintf("%.1f", sd(est_edge)/sqrt(n_rep)), "\n", sep = ''))
  cat(paste("true_edge: ",sprintf("%.1f", mean(true_edge)), "$pm$", sprintf("%.1f", sd(true_edge)/sqrt(n_rep)), "\n", sep = ''))
  cat(paste("timeit: ",sprintf("%.1f", mean(timeit)), "$pm$", sprintf("%.1f", sd(timeit)/sqrt(n_rep)), "\n","\n", sep = ''))
}


save_hd <- function(model,lower_signal, eps){
  dir.create(file.path("./Figure2/output_hd"), showWarnings = FALSE)
  save(list = c("hd"),
       file = paste("./Figure2/output_hd/", model, "_" ,lower_signal, "_", eps,".Rdata", sep = ''))
}


log_hd <- function(model,lower_signal, eps){
  # load metrics
  dir.create(file.path("./Figure2/output_hd"), showWarnings = FALSE)
  load(paste("./Figure2/output_hd/", model, "_" ,lower_signal, "_", eps,".Rdata", sep = ''))
  line=paste("The result from model: ", model," // " , "lower_signal: ",lower_signal, "eps: ", eps,"\n",
             "hd: ",sprintf("%.1f", mean(hd)), "$pm$", sprintf("%.1f", sd(hd)/sqrt(n_rep)), "\n",sep='')
  write(line,file=paste("./Figure2/output_hd/log.txt",sep=''),append=TRUE)
}

save_minimap <- function(model, eps){
  dir.create(file.path("./Figure5"), showWarnings = FALSE)
  dir.create(file.path("./Figure5/output"), showWarnings = FALSE)
  save(list = c("hd_directed", "hd_undirected", "fliped"),
       file = paste("./Figure5/output/", model, "_" , eps,".Rdata", sep = ''))
}

log_minimap <- function(model, eps){
  # load metrics
  dir.create(file.path("./Figure5/output"), showWarnings = FALSE)
  load(paste("./Figure5/output/", model, "_" , eps,".Rdata", sep = ''))
  line=paste("The result from model: ", model," // " , "eps: ", eps,"\n",
             "hd_direct: ",sprintf("%.1f", mean(hd_directed)), "$pm$", sprintf("%.1f", sd(hd_directed)/sqrt(n_rep)), "\n",
             "hd_undirect: ",sprintf("%.1f", mean(hd_undirected)), "$pm$", sprintf("%.1f", sd(hd_undirected)/sqrt(n_rep)), "\n",
             "fliped: ",sprintf("%.1f", mean(fliped)), "$pm$", sprintf("%.1f", sd(fliped)/sqrt(n_rep)), "\n","\n",sep='')
  write(line,file=paste("./Figure5/output/log.txt",sep=''),append=TRUE)
}

