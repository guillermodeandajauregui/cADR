#read libraries
source("libs/functions_tripartite.R")

#load tripartite graph 
load(file = "results/directed_graphs_bi_tri.RData")

#Drug Combination Analysis

combination_analysis = tripartite_full_function(g.tri)
save(combination_analysis, "results/combination_analysis.RData")
