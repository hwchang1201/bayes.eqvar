GDS <- function(X, scoreName = "SEMSEV", pars = list(), check = "checkUntilFirstMinK", penFactor = 1, output = FALSE, startAt = "randomGraph",initG=NULL,smallK=F)
# Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
#
# INPUT:
#   X:         n*p matrix of data
#   scoreName: in principle, one can plug in different score functions. standart: "SEMSEV"
#   pars:      contains parameters for score functions. standart: list()
#   check:     standart: "checkUntilFirstMinK"
#   penFactor: for experimentation one can increase the penFactor of the BIC penalization. standart: 1
#   output:    if TRUE, some output is shown. standart: FALSE
#   startAt:   if "emptyGraph", the search algorithm starts at the empty graph. if "randomGraph", ist starts with a random sparse graph.
#
# OUTPUT:
#   list with
#     Adj:     adjacency matrix of final estimate
#     Score:   score of final estimate
#     B:       coefficients of final estimate
#     ... and some more 
{   
    #computeStatisticsFromData <- get(paste("computeStatisticsFromData",scoreName,sep=""))
    #pars <- computeStatisticsFromData(X,pars)
    pars <- computeStatisticsFromDataSEMSEV(X,pars)
    p <- dim(X)[2]
    kvec <- c(1*p,2*p,3*p,5*p,300) 
    if (smallK==T){
      kvec <- max(300,5*p)
    }
    numRandomRestarts <- length(kvec)
    StatesRestarts <- array(list(NULL), numRandomRestarts)
    
    for(starts in 1:numRandomRestarts)
    {
        k <- kvec[starts]
        if(startAt == "emptyGraph")
        {
            # initialize with empty matrix
            initialAdj <- diag(rep(0,p))
        }
        if(startAt == "randomGraph")
        {
          # initialize with random matrix
          initialAdj <- randomDAG(p,2/(p-1))
          initialAdj <- as(initialAdj, "matrix")
          initialAdj <- 1* (abs(initialAdj) > 0)
        }
        if(startAt == "Warm")
        {
          # initialize with a given graph
          initialAdj <- initG
        }
        #initializeState <- get(paste("initialize",scoreName,sep=""))
        #State <- initializeState(X,initialAdj,pars)
        State <- initializeSEMSEV(X,initialAdj,pars)    
        # preparation
        madeStep <- TRUE
        stepNr <- 0
        checkDAGs <- 0
        
        # optimize    
        while(madeStep == TRUE)
        {
            stepNr <- stepNr + 1
            if(output == TRUE)
            {
                cat("\r Step Nr:",stepNr, " ")
                if(p<8)
                {
                    cat("Current DAG:\n ")
                    show(State$Adj)
                }
            cat("Current Score:",State$Score,"\t")
            #cat("Current HD:", hammingDistance(trueG,State$Adj))
            }            

            ####
            # optimization by checking ALL neighbors (slow)
            ####
            if(check == "checkAll")
            {
                afterStep <- oneStepCheckingAll(State,scoreName,pars,X,penFactor,checkDAGs)
            }

            ####
            # some random neighbours
            ####
            if(check == "checkUntilFirstMinK")
            {
                afterStep <- oneStepCheckingUntilFirstMinK(State,scoreName,pars,X,k, penFactor,checkDAGs)
                if(output == TRUE)
                {
                    cat("step is done at ", afterStep$whichStep[1], afterStep$whichStep[2])
                }   
            }

            ####
            # making step after first neighbour with better score has been found
            ####
            if(check == "checkUntilFirst")
            {
                afterStep <- oneStepCheckingUntilFirst(State,scoreName,pars,X,penFactor,checkDAGs)
            }
            State <- afterStep$State
            madeStep <- afterStep$madeStep 
            checkDAGs <- afterStep$checkDAGs
            if(scoreName == "SEMIND")
            {
                if(State$Score < -0.05)
                {
                    madeStep <- FALSE
                }
            }
        }
        if(output == TRUE)
        {
            cat("Number of DAGs tested:",checkDAGs,"\n")
            cat("\n")
        }
        StatesRestarts[[starts]] <- State
    }
    allScores <- rep(-1,numRandomRestarts)
    for(starts in 1:numRandomRestarts)
    {  
        allScores[starts] <- StatesRestarts[[starts]]$Score
    }
    if(output == TRUE)
    {
        cat("All Scores: ")
        show(allScores)
    }
    best <- which.min(allScores)
    finalState <- StatesRestarts[[best]]
#    computeScoreSEMSEV(finalState,SigmaHat,output=TRUE)
    return(finalState)
}

