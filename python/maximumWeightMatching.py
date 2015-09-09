#!/usr/bin/env python
"""
This scirpt provides functions to subset a graph based on its maximum weight matching.
"""
epilog="""
15.04.15 Jonas Ibn-Salem <j.ibn-salem@uni-mainz.de>
"""

import networkx as nx
import matplotlib.pyplot as plt # only for drawing for debug



def plotGraph(G):
    
    nx.draw(G)
    #~ nx.draw_circular(G)
    #~ nx.draw_spectral(G)
    plt.show()

def getMaxWeightMatching(G):
    """Returns the maximum weighted matching of G as dict"""    
    
    # initialize matching as dict
    matching = {}
    
    # iterate over each connected component C of G (for better runtime performance)
    for C in nx.connected_component_subgraphs(G) :
        
        # TODO get maximum cardinality!!!
        
        # compute maximum weight matching for connected component C
        m = nx.max_weight_matching(C, maxcardinality=True) 
        
        # append it to the matching of the entire graph G
        matching.update(m)
    
    return matching
    

def getMaxWeightMatchingFromFile(inFile, outFile):
    """
Reads a graph G form input file inFile. Computes a subgraph GM of G as the maximum weighted matching of G.
Writes the GM to output file
"""    
    
    # parse graph from file    
    G = nx.read_weighted_edgelist(inFile)
    
    # get maximum weighted matching of G
    matching= getMaxWeightMatching(G)
    
    # get subgraph from the matching
    GM = nx.Graph()
    GM.add_edges_from(matching.items())
    
    # write subgraph to output file
    nx.write_edgelist(GM, outFile, delimiter='\t', data=False)
    

def getMaxWeightMatchingAsDict(x, y, weight):
    if not len(x) == len(y) == len(weight):
        error("Input arguments have to have same length")
    
    # convert input arguments to edges as tuple
    edges = [(x[i], y[i], weight[i]) for i in range(len(x))]
    
    G = nx.Graph()
    G.add_weighted_edges_from(edges)
    
    # get subgraph GM as the maximum weighted matching of G
    matching = getMaxWeightMatching(G)
    
    return matching
    


# test file
#~ TEST_FILE = "temp/genePairs.txt.noHeader"
#~ TEST_FILE_OUT = "temp/graphExample.csv.maxMatch.txt"

#~ TEST_FILE = "temp/paralogPairsFull.txt.noHeader"
#~ TEST_FILE_OUT = "temp/paralogPairsFull.txt.noHeader.maxMatch.txt"


#~ getMaxWeightMatchingFromFile(TEST_FILE, TEST_FILE_OUT)

#~ def getDict():
    #~ return {"a":"b", "c":"d"}

#~ TEST_FILE = "temp/graphExample.csv"
    #~ 
#~ G = nx.read_weighted_edgelist(TEST_FILE)
#~ edgeList = [ (u,v,edata['weight']) for u,v,edata in G.edges(data=True) if 'weight' in edata ]
#~ x = [e[0] for e in edgeList]
#~ y = [e[1] for e in edgeList]
#~ weight = [e[2] for e in edgeList]

