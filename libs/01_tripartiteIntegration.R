library(tidyverse)
library(igraph)
library(data.table)

#Network loading and prep
#tripartite integration

##########################
#Network loading and prep
##########################

#DRUG -> GO  network
dg.nw = read.graph(file = "results/anx_zubdrug_apw_g.gml", format = "gml")
#GO   -> ADR network
ga.nw = read.graph(file = "results/anx_GO_ADR_network.gml", format = "gml")
#drug dictionary 
load(file = "data/unique_dict.RData")

#add names to vertices
V(dg.nw)$name =V(dg.nw)$label
V(ga.nw)$name =V(ga.nw)$xref
##

#add node types and color
V(dg.nw)$type = ifelse(test = grepl(pattern = "GO:", x = V(dg.nw)$name),
                       yes = "GO",
                       no =  "DRUG"
)

V(ga.nw)$type = ifelse(test = grepl(pattern = "GO:", x = V(ga.nw)$name),
                       yes = "GO",
                       no =  "ADR"
)

V(dg.nw)$color = ifelse(test = grepl(pattern = "GO:", x = V(dg.nw)$name),
                        yes = "cornflowerblue",
                        no =  "blueviolet"
)

V(ga.nw)$color = ifelse(test = grepl(pattern = "GO:", x = V(ga.nw)$name),
                        yes = "cornflowerblue",
                        no =  "darksalmon"
)


#add drug common names 
m0 = V(dg.nw)[type=="DRUG"]$name
m1 = match(x = m0, table = my_unique_dict$pert_id)
mx = my_unique_dict$pert_iname[m1]
V(dg.nw)[type=="DRUG"]$drugName = mx
V(dg.nw)$drugName
rm(m0, m1, mx)

#
save(ga.nw, dg.nw, file = "results/drug_go_adr_bipartites.RData")

########################
#tripartite integration
########################

#make directed graphs 

is.directed(dg.nw)
dg.d = as.directed(dg.nw, "mutual")
length(E(dg.d))#517586
#keep only edges FROM DRUG TO GO
dg.d = subgraph.edges(dg.d, eids = E(dg.d)[V(dg.d)[type == "DRUG"]%->%V(dg.d)[type == "GO"]])
length(E(dg.d))# 258793
E(dg.d)
#keep only edges FROM GO TO ADR
ga.d = as.directed(ga.nw, "mutual")
length(E(ga.d))#840
#keep only edges FROM DRUG TO GO
ga.d = subgraph.edges(ga.d, eids = E(ga.d)[V(ga.d)[type == "GO"]%->%V(ga.d)[type == "ADR"]])
length(E(ga.d))# 419

#merge directed graphs
g.tri = igraph::union(dg.d, ga.d)

#in the merged graph, _1 indicates attributes from the DRUG->GO nw 
#                     _2 indicates attributes from the GO->ADR  nw 
#Edge attribute "value" is the weight in the DRUG->GO nw (PAEA score)
#Edge attribute "featureselection" is the weight in the GO->ADR nw (from Ma'yan original work)

#calculate simple tripartite parameters
V(g.tri)$degree.tri.in = degree(graph = g.tri, mode = "in")
V(g.tri)$degree.tri.out = degree(graph = g.tri, mode = "out")
V(g.tri)$betweenness.tri = betweenness(g.tri)

#make types and colors for tri

V(g.tri)$type = ifelse(is.na(V(g.tri)$type_1), yes = V(g.tri)$type_2, no = V(g.tri)$type_1)
V(g.tri)[type == "DRUG"]$color = "blueviolet"
V(g.tri)[type == "GO"]$color = "cornflowerblue"
V(g.tri)[type == "ADR"]$color = "darksalmon"

#remove "size"
g.tri = delete_vertex_attr(g.tri, "size")

save(dg.d, ga.d, g.tri, file = "results/directed_graphs_bi_tri.RData")
