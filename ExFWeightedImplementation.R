## The heart of the algorithm is to iterate over all possible
## infection clusters of size two, storing their FI (force of
## infection, or, more exactly, the sum of the weighted degree
## of the cluster). The computed FI value is then scaled by the
## probability that it can form, by normalizing the probability
## that the path taken was actually chosen. The return value is
## the entropy of the normalized FI values.
##
## This file contains two implementation of the ExF:
##   ExFW and
##   ExFWD
## The first is for weighted graphs, the second is for weighted and directed
## graphs.
## To use:
## source("ExpectedForceWeighted.R") ## loads the functions into your current workspace
## ExFvals <- sapply(1:vcount(gr),ExFW, gr) ## computes the ExF of each node in the graph
##   or you can call it on one specific node, i.e.
## val <- ExF(qnode, gr) ## where qnode is the query node id


ExFW <- function(qnode,graph){
  if(! "weight" %in% list.edge.attributes(graph)){
    stop("Error: graph edges must have a weight attribute.") }
  .FI <- function(graph,clust){
    fi <- sum(E(graph)[adj(clust)]$weight) -
      sum(E(graph)[clust %--% clust]$weight)
    ##cat("clust",clust,"\t degree:",fi,"\n")
    return(fi) }
  graph <- simplify(graph) ## remove self edges
  ## Get all neighbors of the querry node at distance one and two:
  neigh <- graph.bfs(graph,qnode,order=FALSE,
                     dist=TRUE,unreachable=FALSE,
                     callback=function(graph,data,extra){
                       data["dist"]==3})$dist
  ## vector of nodes at distance one
  d.one.nodes <- which(neigh==1)
  n.d.one <- length(d.one.nodes)
  w.d.one <- sum(E(graph)[adj(qnode)]$weight)
  ##cat("w.d.one:",w.d.one,"\n")
  ## vector of nodes at distance two
  d.two.nodes <- which(neigh==2)

  ## SAFETY CHECKS
  if(length(E(graph)[adj(qnode)]) == 0){
    ## cat("no outbound edges\n")
    return(-1)
  }
  if(length(d.two.nodes)==0){
    ## cat("no nodes at distance two\n")
    return(0)
  }

  ## pre-allocate the vector of FI values
  guestimated.numFI <- 2*sum(n.d.one*length(d.two.nodes))
  allFI <- numeric(guestimated.numFI+5)
  numFI <- 0;  totalFI <- 0

  ## The iteration is over all nodes at distance one from the source,
  ##     within this loop we consider both all remaining d.one.nodes
  ##     and all d.two.nodes reachable from the current d.one.node.
  for(i in 1:n.d.one){
    ## we will need the following edge weight (sums) to scale the FI:
    w.q.i <- E(graph, c(qnode,d.one.nodes[i]))$weight
    w.d.i <- sum(E(graph)[adj(d.one.nodes[i])]$weight)
    ##cat(" *w.q.i:",w.q.i,"\n")
    ##cat("  w.d.i:",w.d.i,"\n")
    if(i<n.d.one){   ## all remaining d.one.nodes
      for(j in (i+1):n.d.one){
        ## Increase storage, if necessary: (code for optimization only)
        if(numFI>guestimated.numFI){
          guestimated.numFI <-  round(1.5 * guestimated.numFI)
          foo <- allFI
          allFI <- vector(mode="numeric",length=guestimated.numFI+5)
          allFI[1:numFI] <- foo[1:numFI]
        } ## END increase storage
        ## compute cluster FI
        clustFI <- .FI(graph, c(qnode, d.one.nodes[c(i,j)]))
        ## store once for each way it can form, scaling by the probability
        ## of this path; we need the following edge weight (sums)
        w.q.j <- E(graph, c(qnode,d.one.nodes[j]))$weight
        w.d.j <- sum(E(graph)[adj(d.one.nodes[j])]$weight)
        ##cat("    w.q.j:",w.q.j,"\n")
        ##cat("    w.d.j:",w.d.j,"\n")
        ## way 1: qnode -- i, qnode -- j
        wmult <- w.q.i/w.d.one * w.q.j/(w.d.one + w.d.i - 2*w.q.i)
        numFI <- numFI+1
        allFI[numFI] <- clustFI*wmult
        totalFI <- totalFI + clustFI*wmult
        ## way 2: qnode -- j, qnode -- i
        wmult <- w.q.j/w.d.one * w.q.i/(w.d.one + w.d.j - 2*w.q.j)
        numFI <- numFI+1
        allFI[numFI] <- clustFI*wmult
        totalFI <- totalFI + clustFI*wmult
        ## if an i -- j edge, two more ways:
        if(d.one.nodes[j] %in% neighborhood(graph,d.one.nodes[i],
                                            order=1)[[1]]){
          ## one more edge weight needed
          w.i.j <- E(graph, c(d.one.nodes[i],d.one.nodes[j]))$weight
          ##cat("      w.i.j:",w.i.j,"\n")
          ## way 3: qnode -- i, i -- j
          wmult <- w.q.i/w.d.one * w.i.j/(w.d.one + w.d.i - 2*w.q.i)
          numFI <- numFI+1
          allFI[numFI] <- clustFI*wmult
          totalFI <- totalFI + clustFI*wmult
          ## way 4: qnode -- j, j -- i
          wmult <- w.q.j/w.d.one * w.i.j/(w.d.one + w.d.j - 2*w.q.j)
          numFI <- numFI+1
          allFI[numFI] <- clustFI*wmult
          totalFI <- totalFI + clustFI*wmult
        }}} ## end all remaining d.one.nodes
    for(dtn in d.two.nodes){   ## all d.two.nodes
      if(dtn  %in% neighborhood(graph,d.one.nodes[i],order=1)[[1]]){
      ## If an edge to the current d.one.node
      ## increase storage, if necessary: (code for optimization only)
      if(numFI>guestimated.numFI){
        guestimated.numFI <-  round(1.5 * guestimated.numFI)
        foo <- allFI
        allFI <- vector(mode="numeric",length=guestimated.numFI+5)
        allFI[1:numFI] <- foo[1:numFI]
      } ## END increase storage
      ## compute cluster FI
      clustFI <- .FI(graph, c(qnode, d.one.nodes[i], dtn))
      ## one more edge weight needed
      w.i.j <- E(graph, c(d.one.nodes[i],dtn))$weight
      ##cat("  w.i.j:",w.i.j,"\n")
      wmult <- w.q.i/w.d.one * w.i.j/(w.d.one + w.d.i - 2*w.q.i)
      numFI <- numFI+1
      allFI[numFI] <- clustFI * wmult
      totalFI <- totalFI + clustFI * wmult
    }}
  } ## end looping over all nodes at distance one
  ## calculate the entropy, note that this clips allFI at numFI
  norm <- allFI[1:numFI]/totalFI
  ##cat("clust degs:",round(allFI[1:numFI],digits=5),"\n")
  ##cat("num clusts:",numFI,"\n")
  return(-sum(norm*log(norm)))
}

ExFWD <- function(qnode,graph){
  if(! is.directed(graph)){
    stop("Caution: This function expects a directed graph and may produce inconsistent results for undirected graphs.")
  }
  if(! "weight" %in% list.edge.attributes(graph)){
    stop("Error: graph edges must have a weight attribute.") }
  graph <- simplify(graph) ## remove self edges
  .FI <- function(graph,clust){
    fi <- sum(E(graph)[from(clust)]$weight) -
      sum(E(graph)[clust %->% clust]$weight)
    ##cat("clust",clust,"\t degree:",fi,"\n")
    return(fi) }
  ## Get all neighbors of the querry node at distance one and two:
  neigh <- graph.bfs(graph,qnode,order=FALSE,dist=TRUE,
                     neimode="out",unreachable=FALSE,
                     callback=function(graph,data,extra){
                       data["dist"]==3})$dist
  ## vector of nodes at distance one
  d.one.nodes <- which(neigh==1)
  n.d.one <- length(d.one.nodes)
  w.d.one <- sum(E(graph)[from(qnode)]$weight)
  ##cat("w.d.one:",w.d.one,"\n")
  ## vector of nodes at distance two
  d.two.nodes <- which(neigh==2)
  ## cat("node:",qnode,'\n')
  ## cat('\t',d.one.nodes,'\n')
  ## cat('\t',d.two.nodes,'\n')
  if(length(E(graph)[from(qnode)]) == 0){
    ## cat("no outbound edges\n")
    return(-1)
  }
  if(length(d.two.nodes)==0){
    ## cat("no nodes at distance two\n")
    return(0)
  }

  ## pre-allocate the vector of FI values
  guestimated.numFI <- 2*sum(n.d.one*length(d.two.nodes))
  allFI <- numeric(guestimated.numFI+5)
  numFI <- 0;  totalFI <- 0

  ## The iteration is over all nodes at distance one from the source,
  ##     within this loop we consider both all remaining d.one.nodes
  ##     and all d.two.nodes reachable from the current d.one.node.
  for(i in 1:n.d.one){
    ## we will need the following edge weight (sums) to scale the FI:
    w.q.i <- E(graph, c(qnode,d.one.nodes[i]))$weight
    if(are.connected(graph,d.one.nodes[i],qnode)){
      w.i.q <- E(graph, c(d.one.nodes[i],qnode))$weight
    } else {w.i.q <- 0 }
    w.d.i <- sum(E(graph)[from(d.one.nodes[i])]$weight)
    ##cat(" *w.q.i:",w.q.i,"\n")
    ##cat("  w.d.i:",w.d.i,"\n")
    if(i<n.d.one){   ## all remaining d.one.nodes
      for(j in (i+1):n.d.one){
        ## Increase storage, if necessary: (code for optimization only)
        if(numFI>guestimated.numFI){
          guestimated.numFI <-  round(1.5 * guestimated.numFI)
          foo <- allFI
          allFI <- vector(mode="numeric",length=guestimated.numFI+5)
          allFI[1:numFI] <- foo[1:numFI]
        } ## END increase storage
        ## compute cluster FI
        clustFI <- .FI(graph, c(qnode, d.one.nodes[c(i,j)]))
        ## store once for each way it can form, scaling by the probability
        ## of this path; we need the following edge weight (sums)
        w.q.j <- E(graph, c(qnode,d.one.nodes[j]))$weight
        if(are.connected(graph,d.one.nodes[j],qnode)){
          w.j.q <- E(graph, c(d.one.nodes[j],qnode))$weight
        } else { w.j.q <- 0 }
        w.d.j <- sum(E(graph)[from(d.one.nodes[j])]$weight)
        ##cat("    w.q.j:",w.q.j,"\n")
        ##cat("    w.d.j:",w.d.j,"\n")
        ## way 1: qnode -- i, qnode -- j
        wmult <- w.q.i/w.d.one * w.q.j/(w.d.one - w.q.i + w.d.i - w.i.q)
        numFI <- numFI+1
        allFI[numFI] <- clustFI*wmult
        totalFI <- totalFI + clustFI*wmult
        ## way 2: qnode -- j, qnode -- i
        wmult <- w.q.j/w.d.one * w.q.i/(w.d.one - w.q.j + w.d.j - w.j.q)
        numFI <- numFI+1
        allFI[numFI] <- clustFI*wmult
        totalFI <- totalFI + clustFI*wmult
        ## way 3: qnode -- i, i -- j
        if(are.connected(graph,d.one.nodes[i],d.one.nodes[j])){
          ## one more edge weight needed
          w.i.j <- E(graph, c(d.one.nodes[i],d.one.nodes[j]))$weight
          wmult <- w.q.i/w.d.one * w.i.j/(w.d.one - w.q.i + w.d.i - w.i.q)
          numFI <- numFI+1
          allFI[numFI] <- clustFI*wmult
          totalFI <- totalFI + clustFI*wmult
        }
        ## way 4: qnode -- j, j -- i
        if(are.connected(graph,d.one.nodes[j],d.one.nodes[i])){
          w.j.i <- E(graph, c(d.one.nodes[j],d.one.nodes[i]))$weight
          wmult <- w.q.j/w.d.one * w.j.i/(w.d.one - w.q.j + w.d.j - w.j.i)
          numFI <- numFI+1
          allFI[numFI] <- clustFI*wmult
          totalFI <- totalFI + clustFI*wmult
        }
      }} ## end all remaining d.one.nodes
    for(dtn in d.two.nodes){   ## all d.two.nodes
      ##if(dtn  %in% neighborhood(graph,d.one.nodes[i],order=1)[[1]]){
      ## If an edge to the current d.one.node
      if(are.connected(graph,d.one.nodes[i],dtn)){
        ## increase storage, if necessary: (code for optimization only)
        if(numFI>guestimated.numFI){
          guestimated.numFI <-  round(1.5 * guestimated.numFI)
          foo <- allFI
          allFI <- vector(mode="numeric",length=guestimated.numFI+5)
          allFI[1:numFI] <- foo[1:numFI]
        } ## END increase storage
        ## compute cluster FI
        clustFI <- .FI(graph, c(qnode, d.one.nodes[i], dtn))
        ## one more edge weight needed
        w.i.j <- E(graph, c(d.one.nodes[i],dtn))$weight
        ##cat("  w.i.j:",w.i.j,"\n")
        wmult <- w.q.i/w.d.one * w.i.j/(w.d.one - w.q.i + w.d.i - w.i.q)
        numFI <- numFI+1
        allFI[numFI] <- clustFI * wmult
        totalFI <- totalFI + clustFI * wmult
      }}
  } ## end looping over all nodes at distance one
  ## calculate the entropy, note that this clips allFI at numFI
  norm <- allFI[1:numFI]/totalFI
  ##cat("clust degs:",round(allFI[1:numFI],digits=5),"\n")
  ##cat("num clusts:",numFI,"\n")
  ##cat(round(allFI[1:numFI],digits=2),'\n')
  return(-sum(norm*log(norm)))
}
