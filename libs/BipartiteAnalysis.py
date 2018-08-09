# -*- coding: utf-8 -*-
"""
Script to take a nx-style gml of a bipartite graph, calculate 
"""

import networkx as nx
import argparse
parser =  argparse.ArgumentParser()
parser.add_argument("graph")
parser.add_argument("out")
args = parser.parse_args()
path = args.graph
path_out = args.out

B = nx.read_gml(path)
print("READ")
nx.set_node_attributes(B,
                       values = dict(nx.degree(B)),
                       name = "grado" 
                       )
print("GRADO SETTED")
#F = B.subgraph((n for n in B if B.node[n]['grado']>1))
F = B
F = F.to_undirected()
print("UNDIRECTIZED")
print(F)
print(nx.is_bipartite(F))
#clustering coefficients

nx.set_node_attributes(F, 
                       values = nx.bipartite.clustering(F, mode = "dot"),
                       name = "ClusteringDot"
                       )

print("CLUSDOTED")
nx.set_node_attributes(F, 
                       values = nx.bipartite.robins_alexander_clustering(F),
                       name   = "ClusteringRobbinsAlexander"
                       )

print("ROBBINSED")
#if redundancy fails with  Cannot compute redundancy coefficient for a node that has fewer than two neighbors
redundancia= dict()
for nodo in F.nodes():
    if F.node[nodo]["grado"]<3:
        redundancia[nodo] = 0
    else:
        redundancia[nodo] = nx.bipartite.node_redundancy(F, nodes = [nodo])[nodo]

        
nx.set_node_attributes(F, 
                       values = redundancia,
                       name ="Redundancy" 
                       )
print("REDUNDANCED")
nx.write_gml(F, path_out)
#nx.write_graphml(F, path_out)
print("victoria!")
