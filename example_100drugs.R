#first figure gml 
my_nodes = sample(V(g.tri)[type=="DRUG"], 100)
my_nhood = unique(unlist(neighborhood(g.tri, order = 2, my_nodes, "out")))
my_xampl = induced_subgraph(g.tri, my_nhood)
my_xampl = remove_headless_go(my_xampl)
my_xampl = remove_tailess_go(my_xampl)
happyPlot(my_xampl)

#cytoscape is very picky with the gml, so we have to make some adaptations 
names2remove = grep(pattern = "_2", x = names(vertex.attributes(my_xampl)), value = TRUE)
ready2write = my_xampl
for(i in names2remove){
  ready2write = delete_vertex_attr(graph = ready2write, name = i)
}
names2remove = grep(pattern = "_1", x = names(vertex.attributes(ready2write)), value = TRUE)
for(i in names2remove){
  ready2write = delete_vertex_attr(graph = ready2write, name = i)
}
names(edge.attributes(ready2write))[1]<-"neglogp"
E(ready2write)$neglogp[is.na(E(ready2write)$neglogp)] <-777
E(ready2write)$neglogp[is.na(E(ready2write)$neglogp)] <-777
E(ready2write)$featureselection[is.na(E(ready2write)$featureselection)] <-777
write.graph(ready2write, "results/example_100.gml", "gml")
