#example 3: metformin, captopril, omeprazole

V(g.tri)[grep("metfor", V(g.tri)$drugName)]$drugName
grep("metfor", V(g.tri)$drugName)#1764

V(g.tri)[grep("captopril", V(g.tri)$drugName)]$drugName
grep("captopril", V(g.tri)$drugName)#2707

grep("omeprazole", V(g.tri)$drugName)#2298


my_nodes = unique(unlist(neighborhood(g.tri, 
                                      order = 2, 
                                      nodes = c(2707,2298,1764), 
                                      mode = "out")
                         )
                  )
my_example_3 = induced_subgraph(g.tri, my_nodes)

pair_composites_3 = function(g){
  #takes a drug subgraph for three drugs
  #identifies which ADRs are composite ADRs 
  m1 = neighborhood(graph = g, 
                    order = 2, 
                    nodes = V(g)[type=="DRUG"][1], 
                    mode = "out", 
                    mindist = 2)[[1]]
  
  m2 = neighborhood(graph = g, 
                    order = 2, 
                    nodes = V(g)[type=="DRUG"][2],  
                    mode = "out", 
                    mindist = 2)[[1]]
  
  m3 = neighborhood(graph = g, 
                    order = 2, 
                    nodes = V(g)[type=="DRUG"][3],  
                    mode = "out", 
                    mindist = 2)[[1]]
  
  i12  = intersect(m1,m2)
  i13  = intersect(m1,m3)
  i23  = intersect(m2,m3)
  
  ui123 = unique(c(i12, i13, i23))
  
  toDelete = setdiff(V(g)[type=="ADR"], ui123)
  #print(toDelete)
  return(toDelete)
}

pair_composites_3(my_example_3)
my_example_3 = remove_tailess_go(delete.vertices(my_example_3, pair_composites_3(my_example_3)))
happyPlot(my_example_3)

components(my_example_3)
length(V(my_example_3)[type=="ADR"])
length(V(my_example_3)[type=="GO"])
length(V(my_example_3)[type=="DRUG"])
components(find_CRM(my_example_3))


#cytoscape is very picky with the gml, so we have to make some adaptations 
names2remove = grep(pattern = "_2", x = names(vertex.attributes(my_example_3)), value = TRUE)
ready2write = my_example_3
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

degree(graph = my_example_3, v = V(my_example_3)[type=="GO"], mode = "in")
degree(graph = my_example_3, v = V(my_example_3)[type=="GO"], mode = "out")
degree(graph = my_example_3, v = V(my_example_3)[type=="GO"], mode = "out")

write.graph(ready2write, "results/example_3.gml", "gml")
