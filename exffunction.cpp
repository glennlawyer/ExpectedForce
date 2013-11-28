/* exffunction.cpp
 * 
 * not a stand-alone! This file is sourced
 * from 
 * which compiles it into an R function
 *
 * Computing the ExF requires knowledge of the local
 network topology. We use R/Igraph to extract the
 node neighborhood. The code here parses this R data
 into a C++ representation, then uses this 
 representation to determine the ExF. Hence I divide
 the code into "helper functions" which build this
 representation, and the main function "ExF" which
 provides the interface to R.
 
 The graph is stored as an indexed edgelist. The
 index tracks the start/end of each node in the 
 "ego" portion of the edgelist.

 I get cleaner code by using a FULL edge list, 
 i.e.  containing both 1--2 and 2--1, then by
 a more compact list which implies bi-directionality;
 a compact list would probably be faster.
 
 Also, I am not sure if creating the edge index
 saves any time. Since the edge list is assumed
 to be sorted, it may be faster to binary search
 the edge list for the start positions, than to
 create/store the index. Again, however, this 
 creates more complex code.
 */

#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <queue>
#include <map>
#include <iostream>
#include <iterator> // for ostream_iterator, used for debugging output
 
typedef std::vector<int> svi;
typedef std::vector<int>::iterator svii;
typedef std::map<int,int> smi;
typedef std::map<int,int>::iterator smii;

//////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS 
// set_egostarts: Creates a node-based index to the edgelist
// get_neighbors: Finds neighbors from a seed at distance 1 and 2
// cluster_degree: Returns the degree of a given cluster.
//////////////////////////////////////////////////////////////////////


/* Given a reference to the index, and a 
   SORTED, FULL edgelist, create the index 
   and return 0 (success).

   @param[out] egostart, egoend. The index
   @param[in] egos, alters. The edgelist
   @return success code, in case we later
   wish to add a safety check
 */
int set_egostarts(svi & egostart, svi & egoend,
		  const svi & egos,  const svi & alters){
  egostart.clear(); egoend.clear();
  int i,nnum=0,cnt=0;
  egostart.push_back(cnt);
  while(nnum<egos[0]){ // no edge
    egostart.push_back(cnt);
    egoend.push_back(cnt);
    nnum++;
  } 
  for(i=1;i<egos.size();i++){
    cnt++;
    if(egos[i] != egos[i-1]){
      //std::cout<<egos[i]<<"  ";
      egostart.push_back(cnt); 
      egoend.push_back(cnt);
      nnum++;
      while(nnum<egos[i]){ // no edge
	//std::cout<<"*";
	nnum++;
      }}
  }
  egoend.push_back(++cnt);
  return 0;
}


/* Find the neighbors at distance one and two
   from a seed in a network.

   @param[out] dOne, dTwo. References to the vectors which 
   store the neighbor node ids
   @param[in] seed. The seed node
   @param[in] egostart,egoend, alters. The needed network topology
   @return success (usefull for debugging, or adding safety checks)
*/
int get_neighbors(svi & dOne, svi & dTwo,
		  const int seed,
		  const svi & egostart, const svi & egoend,
		  const svi & alters){
  dOne.clear(); dTwo.clear();  
  std::queue<int> q;
  smi visited; // maps, as I index by node #
  smi distance; 
  int i,curN,curD;
    
  // initialize the queue
  q.push(seed);
  visited[seed]=1;
  distance[seed]=0;
  //std::cout<<egostart.size();
  while(! q.empty()){ // run the BFS
    curN=q.front(); q.pop();  // grab the next node 
    //std::cout<<curN<<"  "<<distance[curN]<<std::endl;
    switch(distance[curN]){ // store the results, or exit if done
    case 1: dOne.push_back(curN); break;
    case 2: dTwo.push_back(curN); break;
    case 3: return 0; // STOP CONDITION (first node at dist 3
    default: break; // should hit this for the seed node, dist 0      
    }
    for(i=egostart[curN];i<egoend[curN];i++){ // add its children to the queue
      curD=alters[i];
      if(! visited[curD]){
	visited[curD]=1;
	distance[curD]=distance[curN]+1;
	q.push(curD);
      }}
  }
  // the BFS should only get here if it finds no nodes at distance 3
  //std::cout<<"found "<<dOne.size()+dTwo.size()<<" neighbors."<<std::endl;
  //std::cout<<"get_neighbors had a problem."<<std::endl;
  return 0; 
}


/* Find the degree of a cluster

   returned as a double so we don't have to
   cast it later.
*/ 
double cluster_degree(svi & clusterNodes,
		      const svi & egostart, const svi & egoend,
		      const svi & alters) {
  int degree=0;
  // for each element of clusterNodes,
  //   for each neighbor of that element,
  //      which is not in clusterNodes
  //         increase the degree
  for(svii clusteriter=clusterNodes.begin();
      clusteriter!=clusterNodes.end();clusteriter++){
    for(int i=egostart[*clusteriter];i<egoend[*clusteriter];i++){
      if(std::find(clusterNodes.begin(),clusterNodes.end(),alters[i]) == 
	 clusterNodes.end()){
	degree++;
      }}}
  return(degree);
}


////////////////////////////////////////////////////////
// MAIN FUNCTION
//

/* Given a FULL, SORTED edgelist and a seed node,
   return the ExF of the seed.
   It is easy to create a full, sorted edgelist from
   R/IGraph by i.e. 
   elist <- get.edgelist(graph)
   elist <- rbind(elist,elist[,c(2,1)]
   eorder <- order(elist[,1],elist[,2])
   
   Then call:
   exfcpp(elist[eorder,1],elist[eorder,2],seed)

   Note that the R function
   defined in 
   does this for you.
   
   Note also that only the local neighborhood of 
   the seed is needed (and that extracting this 
   using Igraph may change node indexing)
*/

// [[Rcpp::export]]
double exfcpp(svi _egos, svi _alters, int seed){
  // SAFETY: check if the seed is in the edgelist!
  // SAFETY: check that the edge list is complete and sorted
  /////////////////////////////////////////////////////////////
  // set up the graph structure 
  svi egostart,egoend; // indexes into egos, alters
  svi dOne, dTwo; // nodes at distance 1 (resp 2) from seed
  int all_ok;
  all_ok = set_egostarts(egostart, egoend, _egos, _alters);
  //std::cout<<"egostarts ok "<<all_ok<<std::endl;
  if(all_ok != 0){ return -1;}
  all_ok=get_neighbors(dOne, dTwo, seed, egostart, egoend, _alters);
  //std::cout<<"get neighbors "<<all_ok<<std::endl;
  if(all_ok != 0){ return -2;}
  /////////////////////////////////////////////////////////////
  // Initialize the vectors and etc to store the FI values 
  // for each cluster
  svi tmp(3); // stores the nodes in the cluster
  tmp[0]=seed;
  svii i,j; // iterators over the neighbors of the seed
  double clustFI, totalFI; // cluster FI and total FI
  std::vector<double>  FIvalues; // the vector of FI values
  FIvalues.reserve(1000); // faster to use a constant that to guestimate
  ///////////////////////////////////////////////////////////////
  // Iterate over all possible clusters of size 2 (plus the seed).
  // The iteration is over all nodes at distance one from the source,
  //     within this loop we consider both all remaining dOne nodes
  //     and all dTwo nodes reachable from the current dOne node.
  for(i=dOne.begin();i!=dOne.end();i++){ 
    tmp[1]=*i;
    for(j=i+1;j!=dOne.end();j++){ // the remaining dOne nodes
      tmp[2]=*j;
      clustFI=cluster_degree(tmp,egostart,egoend,_alters);
      // add it once for each way the cluster could form
      FIvalues.push_back(clustFI); totalFI+=clustFI;
      FIvalues.push_back(clustFI); totalFI+=clustFI;
      for(int edgeindx=egostart[*i];edgeindx<egoend[*i];edgeindx++){ 
	if(_alters[edgeindx]==*j){
	  FIvalues.push_back(clustFI); totalFI+=clustFI;
	  FIvalues.push_back(clustFI); totalFI+=clustFI;
	}}
    }
    // now search for all neighbors of i at distance two
    for(int neigh=egostart[*i];neigh<egoend[*i];neigh++){ 
      j=find(dTwo.begin(),dTwo.end(),_alters[neigh]);
      if(j != dTwo.end()){
	tmp[2]=*j;
	clustFI=cluster_degree(tmp,egostart,egoend,_alters);
	FIvalues.push_back(clustFI); totalFI+=clustFI;	
      }}
  } // end iteration over all clusters
  /////////////////////////////////////////////////////////////
  // compute "entropy" of the FI values
  double normalizedFI, ExF(0); 
  //std::cout<<"clust degs: ";
  for(std::vector<double>::iterator i=FIvalues.begin();i!=FIvalues.end();i++){
    //std::cout<<*i<<" ";
    normalizedFI =*i/totalFI;
    ExF -= (log(normalizedFI)*normalizedFI);
  }
  //std::cout<<std::endl;
  return(ExF);
}
  







