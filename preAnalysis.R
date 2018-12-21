#read libraries
source("libs/functions_tripartite.R")

print("libraries loaded")

#load tripartite graph 
load(file = "results/directed_graphs_bi_tri.RData")
print("graph loaded")

#Drug Combination Analysis

combination_analysis = tripartite_full_function(g.tri, kores_o = 6, kores_i = 1)
save(combination_analysis, "results/combination_analysis.RData")
print("done")
