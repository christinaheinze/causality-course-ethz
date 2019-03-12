plotGraphfromAdj <- function(Adj, labels = 1:dim(Adj)[1], is_weighted = NULL)
    # Copyright (c) 2013 - 2013  Jan Ernest [ernest@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
    
    # Input:  Adj         - Adjacency matrix of a graph
    #         is_weighted - TRUE if the plot should contain edge weights
    # Output: Plot of the corresponding (un-)weighted directed graph
    
{
    library(igraph)
    
    if(!is.null(is_weighted)) 
    { 
        gr <- graph.adjacency(Adj, mode = "directed", weighted = is_weighted, diag = FALSE)	
        E(gr)$label <- E(gr)$weight 
    } else 
    {
        Adj[Adj != 0] <- 1 
        gr <- graph.adjacency(Adj, mode = "directed", weighted = is_weighted, diag = FALSE)
#        gr <- graph.adjacency(Adj, mode = "undirected", weighted = is_weighted, diag = FALSE)
    }
    
    V(gr)$label <- labels 
    
    plot(gr)
}




plotCausalOrderedDAGfromAdj <- function(Adj, labels = 1:dim(Adj)[1], main=NULL)
    # Copyright (c) 2013 - 2013  Jan Ernest [ernest@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
    
    # Input:  Adj         - Adjacency matrix of a graph
    #         labels      - an optional vector of labels of the nodes in the graph
    #         
    # Output: Plot of the corresponding (un-)weighted directed graph respecting the causal order. 
    
{
    library(Rgraphviz)
    G <- as(Adj, "graphNEL")
    
    z <- labels
    names(z) = nodes(G)
    nAttrs <- list()
    nAttrs$label <- z
    attrs <- list(node = list(shape = "ellipse", fixedsize = FALSE))
    plot(G, nodeAttrs = nAttrs, attrs = attrs, main=main)    
}





####
#gam Regression from mgcv package
####
library(mgcv)
train_gam <- function(X,y,pars = list(numBasisFcts = 10))
{
  if(is.null(X)){
    result <- list()
    result$Yfit <- as.matrix(rep(mean(y), length(y)))
    result$residuals <- as.matrix(y - result$Yfit)
    result$model <- NA
    result$df <- NA     
    result$edf <- NA     
    result$edf1 <- NA     
    result$p.values <- NA
    
    # for degree of freedom see mod_gam$df.residual
    # for aic see mod_gam$aic
    return(result)
  } else {
    if(dim(as.matrix(X))[2] == 0){
      result <- list()
      result$Yfit <- as.matrix(rep(mean(y), length(y)))
      result$residuals <- as.matrix(y - result$Yfit)
      result$model <- NA
      result$df <- NA     
      result$edf <- NA     
      result$edf1 <- NA     
      result$p.values <- NA
      
      # for degree of freedom see mod_gam$df.residual
      # for aic see mod_gam$aic
      return(result)
    } else {
      if(!("numBasisFcts" %in% names(pars) ))
      { 
        pars$numBasisFcts = 10
      }
      p <- dim(as.matrix(X))
      if(p[1]/p[2] < 3*pars$numBasisFcts)
      {
        pars$numBasisFcts <- ceiling(p[1]/(3*p[2]))
        cat("changed number of basis functions to    ", pars$numBasisFcts, "    in order to have enough samples per basis function\n")
      }
      dat <- data.frame(as.matrix(y),as.matrix(X))
      coln <- rep("null",p[2]+1)
      for(i in 1:(p[2]+1))
      {
        coln[i] <- paste("var",i,sep="")
      }
      colnames(dat) <- coln
      labs<-"var1 ~ "
      if(p[2] > 1)
      {
        for(i in 2:p[2])
        {
          labs<-paste(labs,"s(var",i,",k = ",pars$numBasisFcts,") + ",sep="")
          #      labs<-paste(labs,"s(var",i,") + ",sep="")
          # labs<-paste(labs,"lo(var",i,") + ",sep="")
        }
      }
      labs<-paste(labs,"s(var",p[2]+1,",k = ",pars$numBasisFcts,")",sep="")
      # labs<-paste(labs,"s(var",p[2]+1,")",sep="")
      # labs<-paste(labs,"s(var",p[2]+1,", bs = "cc")",sep="") #factor 2 faster
      # labs<-paste(labs,"s(var",p[2]+1,", bs = "cr")",sep="") # factor 2 + eps faster
      # labs<-paste(labs,"lo(var",p[2]+1,")",sep="")
      mod_gam <- FALSE
      try(mod_gam <- gam(formula=formula(labs), data=dat),silent = TRUE)
      if(typeof(mod_gam) == "logical")
      {
        cat("There was some error with gam. The smoothing parameter is set to zero.\n")
        labs<-"var1 ~ "
        if(p[2] > 1)
        {
          for(i in 2:p[2])
          {
            labs<-paste(labs,"s(var",i,",k = ",pars$numBasisFcts,",sp=0) + ",sep="")
          }
        }
        labs<-paste(labs,"s(var",p[2]+1,",k = ",pars$numBasisFcts,",sp=0)",sep="")
        mod_gam <- gam(formula=formula(labs), data=dat)
      }
      result <- list()
      result$Yfit <- as.matrix(mod_gam$fitted.values)
      result$residuals <- as.matrix(mod_gam$residuals)
      result$model <- mod_gam 
      result$df <- mod_gam$df.residual     
      result$edf <- mod_gam$edf     
      result$edf1 <- mod_gam$edf1     
      result$p.values <- summary.gam(mod_gam)$s.pv
      
      # for degree of freedom see mod_gam$df.residual
      # for aic see mod_gam$aic
      return(result)
    }
  }
}
