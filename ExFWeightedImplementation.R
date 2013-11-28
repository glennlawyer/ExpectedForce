## ExF implementation.R
##
## requires: exffunction.cpp  
##
## The heart of the algorithm is to iterate over all possible
## infection clusters of size two, storing their FI (force of
## infection, or, more exactly, the sum of the weighted degree
## of the cluster). The computed FI value is then scaled by the
## probability that it can form, by normalizing the probability
## that the path taken was actually chosen. The return value is
## the entropy of the normalized FI values.
##
## This file contains two implementation of the ExF:
## ExFW.cannonical
## ExFW.fast
## along with:
## benchmarking code/useage example
## profiling code
## a short, development-oriented test case
##
## Both implementation have the same call structure and return value
## INPUTS:
##   qnode -- the querry node
##   graph -- an Igraph graph. Edges MUST have a weight attribute
## RETURN:
##   the expected force of the querry node
## This calling structure allows for natural vectorization
## i.e.  exfvals <- sapply(nodes, ExF.fast, graph=mygraph)
##
## note that the fast version requires that the graph nodes
## have an "index" attribute, which is best created by:
##    V(graph)$index <- 1:vcount(graph)
## this is because of how Igraph returns the local subgraph.
##
## The cannonical version is written in R, with an emphasis on
## readability. As a result it is very slow. This is also why
## it uses for loops rather than vectorization.
##
## The fast version uses C++ code for the required loops, resulting
## in a several order of magnitude speedup. The actual computation of
## the ExF should be just as easy to read, but the code requires
## several helper functions to set up the C++ version of the graph.
## The code is again written for clarity rather than speed. Several
## possible optimizations are noted in the code comments. An astute
## reader could doubtles find more. Certainly the code would be
## faster still if the graph was stored from the get-go in C++.
## Such an implementation would, however, bury the ExF algorithm
## in a rather larger class. 
##
## Note that the code has few error checks. Use at your own risk.
##
## copyright Aug 2013 by Glenn Lawyer
######################################
require(igraph)
require(Rcpp)
cat("Compiling/loading C++ code...")
sourceCpp(file="exfWfunction.cpp")  ## compile/loads the c++ implementation.
cat("finished compiling.\n")


ExFW.cannonical <- function(qnode,graph){
  if(! "weight" %in% list.edge.attributes(gr)){
    stop("Error: graph edges must have a weight attribute.") }
  .FI <- function(graph,clust){
    fi <- sum(E(graph)[adj(clust)]$weight) -
      sum(E(graph)[clust %--% clust]$weight)
    ##cat("clust",clust,"\t degree:",fi,"\n")
    return(fi) }
  ## Get all neighbors of the querry node at distance one and two:
  neigh <- graph.bfs(graph,qnode,order=FALSE,dist=TRUE,
                     callback=function(graph,data,extra){
                       data["dist"]==3})$dist
  ## vector of nodes at distance one
  d.one.nodes <- which(neigh==1)
  n.d.one <- length(d.one.nodes)
  w.d.one <- sum(E(graph)[adj(qnode)]$weight)
  ##cat("w.d.one:",w.d.one,"\n")
  ## vector of nodes at distance two
  d.two.nodes <- which(neigh==2)
  
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
    w.d.i <- sum(E(gr)[adj(d.one.nodes[i])]$weight)
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
        w.d.j <- sum(E(gr)[adj(d.one.nodes[j])]$weight) 
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


ExFW.fast <- function(qnode,graph){
  c(qnode,"\n")
  if(! "index" %in% list.vertex.attributes(graph)){
    msg <- paste("Graph verticies need an index attribute",
                 "to allow id of seed node in subgraphs.\n",
                 "try:\n   V(graph)$index <- 1:vcount(graph)\n")
    stop(msg)
  }
  if(! "weight" %in% list.edge.attributes(graph)){
    msg <- paste("Graph edges need a (positive) weight attribute",
                 "interpreted as transmission probabilities.\n")
    stop(msg)
  }
  ## get the local neighborhood of the node
  ## extract the edgelist and the id of the seed node in this edgelist
  local <- graph.neighborhood(graph,order=3,nodes=qnode)[[1]]
  seed <- which(V(local)$index==qnode)
  elist <- cbind(get.edgelist(local),E(local)$weight)
  elist <- rbind(elist,elist[,c(2,1,3)])
  sorder <- order(elist[,1],elist[,2])
  ## turn to cpp for the heavy lifting
  return(exfWcpp(elist[sorder,1],elist[sorder,2],elist[sorder,3],seed))
}


#########################################################
if(0){ ## development stuff
  source('ExFWImplementation.R')
  require(igraph)
  require(Rcpp)
  sourceCpp(file="exffunction.cpp")

  ## NOTE!! THIS FUNCTION RELIES ON GLOBAL VARIABLES
  ## as set in the rest of this code block, below.
  testplot <- function(seed){
    V(gr)$color <- "lightblue"
    V(gr)$color[seed] <- "red"
    E(gr)$label <- E(gr)$weight
    plot(gr,layout=lout,vertex.size=7)
    cat("exf:",seed,exfWcpp(elist[,1],elist[,2],elist[,3],seed),"\n")
  }
  
  gr <- ba.game(m=2,n=100,directed=FALSE)
  V(gr)$index <- 1:vcount(gr)
  E(gr)$weight <- round(runif(ecount(gr),min=1,max=10))
  lout <- layout.fruchterman.reingold(gr)

  testers <- which(degree(gr)==2) 
  
  elist <- cbind(get.edgelist(gr),E(gr)$weight)
  elist <- rbind(elist,elist[,c(2,1,3)])
  sorder <- order(elist[,1],elist[,2])
  elist <- elist[sorder,]

  sourceCpp(file="exfWfunction.cpp")  ## compile/loads the c++ implementation.
  exfWcpp(elist[,1],elist[,2],elist[,3],8)

  testers <- 1:100
  can <- sapply(testers,ExFW.cannonical,graph=gr)
  fas <- sapply(testers,ExFW.fast,graph=gr)
  all.equal(can,fas)
}
