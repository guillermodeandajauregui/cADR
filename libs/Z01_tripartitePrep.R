library(tidyverse)
library(igraph)
library(data.table)

#ag: analyzed Ma'yan bipartite network 
#ABC drugpert: drugs affecting pathways in ag
Drugs_Apw   = lapply(X = ABC_drugpert, FUN = function(i){
  names(i[["Aind"]])
})

drug_apw_nw = data.frame(source = character(), 
                         target = character(), 
                         value = numeric()
)

for(i in seq_along(Drugs_Apw)){
  ii = Drugs_Apw[[i]]
  name_i = names(Drugs_Apw)[i]
  smallframe = data.frame(source = character(), 
                          target = character(), 
                          value = numeric()
  )
  for(j in ii){
    rowling = data.frame(source = name_i, 
                         target = j, 
                         value = z[name_i, j]
    )
    smallframe = rbind(smallframe, rowling)
    
  }
  drug_apw_nw = rbind(drug_apw_nw, smallframe)
  
}

save(drug_apw_nw, file = "results/drug_apw_nw.RData")

drug_apw_g = igraph::graph_from_data_frame(d = drug_apw_nw, directed = TRUE)
write.graph(drug_apw_g, file = "results/drug_apw_g.gml", format = "gml")

#run bipartite analysis
system(command = "libs/rgml2nx.sh results/drug_apw_g.gml results/nx_drug_apw_g.gml")
system(command = "python libs/BipartiteAnalysis.py results/nx_drug_apw_g.gml results/anx_drug_apw_g.gml")
