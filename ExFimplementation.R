## ExF implementation.R
##
## requires: exffunction.cpp  
##
## The heart of the algorithm is to iterate over all possible
## infection clusters of size two, storing their FI (force of
## infection, or, more exactly, the degree of the cluster).
## Each FI value is stored once for each way the cluster can
## form. The return value is the entropy of the normalized
## FI values.
##
## This file contains two implementation of the ExF:
## ExF.cannonical
## ExF.fast
## benchmarking code/useage example
## profiling code
## a short, development-oriented test case
##
## Both implementation have the same call structure and return value
## INPUTS:
##   qnode -- the querry node
##   graph -- an Igraph graph.
## RETURN:
##   the expected force of the querry node
## The node is given first, for more natural vectorization
## i.e.  exfvals <- sapply(nodes, ExF.fast, graph=mygraph)
##
## note that the fast version requires that the graph nodes
## have an "indx" attribute, which is best created by:
##    V(graph)$indx <- 1:vcount(graph)
## this is because of how Igraph returns the local subgraph.
##
## The cannonical version is written in R. As a result it is very
## slow, even after careful profiling/tuning of the code.
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
sourceCpp(file="exffunction.cpp")  ## compile/loads the c++ implementation.
cat("finished\n")


ExF.cannonical <- function(qnode,graph){
  .FI <- function(graph,clust){
    return(length(unique(unlist(neighborhood(graph,clust,order=1)))) - 3) }
  ## Here, we know the cluster will have three elements, hence the "- 3"
  ## if the cluster size is unknown, better to use "- length(clust)"

  ## Get all neighbors of the querry node at distance one and two:
  neigh <- graph.bfs(graph,qnode,order=FALSE,dist=TRUE,
                     callback=function(graph,data,extra){
                       data["dist"]==3})$dist
  ## vector of nodes at distance one
  d.one.nodes <- which(neigh==1)
  n.d.one <- length(d.one.nodes)
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
        ## is there an edge between nodes i and j?
        mult <- 
          if(length(E(graph)[ d.one.nodes[i] %--% d.one.nodes[j] ])) 4 else 2
        ## store cluster FI the appropriate number of times
        allFI[seq(numFI+1,length.out=mult)] <- clustFI
        totalFI <- totalFI + mult*clustFI
        numFI <- numFI+mult
      }} ## end all remaining d.one.nodes
  for(dtn in d.two.nodes){   ## all d.two.nodes 
    if(length(E(graph)[ d.one.nodes[i] %--% dtn ])){
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
      numFI <- numFI+1
      allFI[numFI] <- clustFI
      totalFI <- totalFI + clustFI
    }}
  } ## end looping over all nodes at distance one
  ## calculate the entropy, note that this clips allFI at numFI
  norm <- allFI[1:numFI]/totalFI
  return(-sum(norm*log(norm)))
}


ExF.fast <- function(qnode,graph){
  c(qnode,"\n")
  if(! "indx" %in% list.vertex.attributes(graph)){
    msg <- paste("Graph verticies need an index attribute",
                 "to allow id of seed node in subgraphs.\n",
                 "try:\n   V(graph)$indx <- 1:vcount(graph)\n")
    stop(msg)
  }
  ## get the local neighborhood of the node
  ## extract the edgelist and the id of the seed node in this edgelist
  local <- graph.neighborhood(graph,order=3,nodes=qnode)[[1]]
  seed <- which(V(local)$indx==qnode)
  elist <- get.edgelist(local)
  elist <- rbind(elist,elist[,c(2,1)])
  sorder <- order(elist[,1],elist[,2])
  ## turn to cpp for the heavy lifting
  return(exfcpp(elist[sorder,1],elist[sorder,2],seed))
}


#########################################################
## test code
if(0){
  require(microbenchmark)
  
  gr <- ba.game(m=2,n=1000,directed=FALSE)
  V(gr)$indx <- 1:vcount(gr)
  ## if you have time to kill, use this set of testers for the
  ## benchmarking. The R code is slow.
  testers <- which(degree(gr)<max(degree(gr)*0.6)) ## all except hubs. 
  ## Otherwise, let me suggest:
  testers <- which(degree(gr)==5) 

  mb <- microbenchmark(foo <- sapply(testers,ExF.cannonical,graph=gr),
                       bar <- sapply(testers,ExF.fast,graph=gr),
                       times=10)
  print(mb)
  ## one run with 5 tests gives this:
  ##         median          uq         max
  ##1    57.08308    57.28477    57.50847
  ##2 12782.11993 12798.37445 13090.11240
  all.equal(foo,bar)

  testers <- which(degree(gr)<max(degree(gr)*0.6)) ## all except hubs. 
  mb <- microbenchmark(sapply(testers,expForce.opt,graph=gr,USE.NAMES=FALSE),
                       times=100)
  ##          lq   median       uq      max
  ##1 1.422428 1.429312 1.469794 1.535797
}

## profiling code
if(0){
  Rprof("exf.prof")
  sapply(testers,ExF.cannonical,graph=gr)
  Rprof(NULL)
  summaryRprof("exf.prof")

  ## or 
  require(profr)
  prof <- profr(lapply(testers,ExF.cannonical,graph=gr),0.01)
  plot(prof)
}


if(0){ ## development stuff
  source('ExFimplementation.R')
  require(igraph)
  require(Rcpp)
  sourceCpp(file="exffunction.cpp")

  ## NOTE!! THIS FUNCTION RELIES ON GLOBAL VARIABLES
  ## as set in the rest of this code block, below.
  testplot <- function(seed){
    V(gr)$color <- "lightblue"
    V(gr)$color[seed] <- "red"
    plot(gr,layout=lout)
    cat("exf:",seed,exf(elist[,1],elist[,2],seed),"\n")
  }
  
  gr <- ba.game(m=2,n=100,directed=FALSE)
  lout <- layout.fruchterman.reingold(gr)
  V(gr)$indx <- 1:vcount(gr)
  elist <- get.edgelist(gr)
  elist <- rbind(elist,elist[,c(2,1)])
  elist <- elist[order(elist[,1],elist[,2]),]
  
  vals <- rep(-1,vcount(gr))
  for(b in 1:vcount(gr)){
    vals[b] <- exfcpp(elist[,1],elist[,2],b)
  }
  
  testplot(1)
}
