##
##  Source the implementation file, which loads Igraph/Rcpp 
## and compiles the function as ExF.fast
##

source('ExFimplementation.R')

#' Plot the graph, highlighting the selected node
#' 
#' Useful for debugging. Shows the graph and the query node,
#' to give a sense of if the ExF of that node is reasonable.
#' 
#' The layout can be given as a function argument. The function
#' returns the layout so one can use the same layout for 
#' plotting mulitple nodes.
#' 
#' @param qnode The index of the query node.
#' @param graphObj The graph to plot.
#' @param lout The layout for the graph.
testplot <- function(qnode, graphObj, lout=NULL){
	if(is.null(lout)){ lout <- layout.fruchterman.reingold(gr) }
	V(graphObj)$color <- "lightblue"
	V(graphObj)$color[qnode] <- "red"
	plot(graphObj,layout=lout)
	cat("exf:",qnode,ExF.fast(qnode, graphObj),"\n")
	invisible(lout)
}

## Make a sample graph
gr <- ba.game(m=2,n=100,directed=FALSE)
lout <- layout.fruchterman.reingold(gr)

## Generate a test plot, saving the graph layout:
lout <- testplot(sample(1:vcount(gr),1), gr)

## get the ExF of a set of nodes:
##    (uncomment only one of these definitions of testNodes)
testNodes <- 1:vcount(gr) 												 ## all nodes
testNodes <- which(degree(gr)<max(degree(gr)*0.6)) ## all except hubs. 
testNodes <- which(degree(gr)==5)                  ## all with a certain degree

## generate the ExF of 
qnode <- 1           ## the first node
qnode <- testNodes[1]  ## the first test node
cat("The ExF of node",qnode,"is", ExF.fast(qnode, gr),"\n")

## generate the ExF of all testNodes
exfVals <- sapply(testNodes,ExF.fast,graph=gr)

## plot ExF vs Degree
pltorder <- order(exfVals)
plot(exfVals[pltorder], degree(gr)[testNodes][pltorder])

