#!/usr/bin/env python
"""
This scirpt provides functions to subset a graph based on its maximum weight matching.

15.04.15 Jonas Ibn-Salem <j.ibn-salem@uni-mainz.de>
"""

import networkx as nx

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
    
