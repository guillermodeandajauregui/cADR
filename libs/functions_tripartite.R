library(igraph)
library(tidyverse)
library(data.table)

pair_subgraph = function(g, v1, v2){
  #takes a graph and two vertices (drugs)
  #returns subset of trigraph with origin in pair of nodes
  n1 = unlist(ego(graph = g, order = 2, nodes = v1, mode = "out"))
  n2 = unlist(ego(graph = g, order = 2, nodes = v2, mode = "out"))
  sg = induced_subgraph(g, vids = c(n1,n2))
  return(sg)
}

pair_composites = function(g){
  #takes a drug subgraph for a pair of drugs
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
  
  return(V(g)[intersect(m1,m2)])
}

remove_headless_go = function(g, v){
  #takes a (sub)graph 
  #and a set of nodes 
  #removes vertices without upstream 
  t1 = which(degree(graph = g, v = v, mode = "in")==0)
  t2 = v[t1]
  sg = delete.vertices(g, t2)
  return(sg)
}

remove_tailess_go = function(g, v = V(g)[type=="GO"]){
  #takes a (sub)graph 
  #and a set of nodes 
  #removes vertices without downstream
  t1 = which(degree(graph = g, v = v, mode = "out")==0)
  t2 = v[t1]
  sg = delete.vertices(g, t2)
  return(sg)
}

composite_subgraph = function(g){
  #takes a drug subgraph for a pair of drugs
  #removes non-composite ADRs 
  toRemove = setdiff(V(g)[type=="ADR"],pair_composites(g))
  #sg = delete.vertices(graph = g, v = !(pair_composites(g)))
  sg = delete.vertices(graph = g, v = toRemove)
  sg = remove_tailess_go(sg)
  return(sg)
}

find_CRM = function(g, precomposite = FALSE){
  #takes a drug subgraph for a pair of drugs
  #identifies the Composite Risk Modules 
  if(precomposite == TRUE){
    sg = composite_subgraph(g = g)
    sg = delete.vertices(sg, V(sg)[type=="DRUG"])
    return(sg)
  }else{
    sg = delete.vertices(g, V(g)[type=="DRUG"])
    return(sg)
  }
}

paths_to_ADR   = function(g, ADR){
  #takes the drug pair subgraph 
  #and an ADR 
  #returns a list of GOs connecting each drug to ADR 
  ll = lapply(V(g)[type=="DRUG"], function(x){
    q = intersect(neighbors(g, v = x, mode = "out"),
                  neighbors(g, v = ADR, mode = "in")
    )
    r = V(g)[q]
    return(r)
  })
  return(ll)
}

pair_ADR_Mode1 = function(g, ADR){
  #takes a drug subgraph for a pair of drugs
  #and an ADR
  #identifies the number of Mode 1 paths from DRUG to ADR
  timp = paths_to_ADR(g, ADR)
  tamp = intersect(timp[[1]], timp[[2]])
  tomp = V(g)[tamp]
  return(tomp)
}

pair_ADR_Mode2 = function(g, ADR){
  #takes a drug subgraph for a pair of drugs
  #and an ADR
  #identifies the number of Mode 2 paths from DRUG to ADR
  Mode_1 = pair_ADR_Mode1(g, ADR)
  Mode_2 = lapply(X = V(g)[type=="DRUG"], function(x){
  }
  )
}

happyPlot = function(g){
  #convenience function for my plotting needs
  plot(g, vertex.label = "", edge.arrow.size = .2)
}

heatmap.tidier = function(df){
  #takes a result data frame
  #returns a data.frame for heatmapping
  mx = as.matrix(df)
  #make matrix symmetric; hack to make clustering easier
  mx = Matrix::forceSymmetric(mx, uplo = "L")
  #clustering order
  ord = hclust(dist(mx))$order
  #make a reordered matrix 
  t6 = as.matrix(mx)[ord,ord]
  #remove lower triangle 
  t6[lower.tri(t6)] = NA
  #make tidy dataframe
  m7 = reshape2::melt(data = t6, na.rm = FALSE)
  return(m7)
}

rest.of.the.heatmap = function(things){
  #ggplot tile 
  p = ggplot(m7, aes(Var1,Var2, fill=value)) + geom_raster()
  p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p = p + labs(fill = filly.name, x = "", y = "")
  #p = p + scale_fill_gradient2(low = "black", mid = "white", high = "blue", limits=c(0,max.value))
  #p = p + lims(fill = c(0, max.value))
  p = p + scale_fill_gradient2(low = "white", mid = "steelblue", high = "violet", midpoint = 50, limits = c(0, max.value))
  p = p + labs(title = legenda, subtitle = sub_legend)
  if(d.names == FALSE){
    p = p + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
    p = p + theme(axis.ticks = element_blank())
  }
  plot(p)
  return(p)
}

#functions for whole network analysis

tripartite_pair_function = function(g, d1, d2){
  #takes a graph with DRUGS, GO. and ADRs
  #and two drugs 
  #returns a big list of results 
  
  result_list = list(ADR = list(), 
                     GO  = list(),
                     cADR   = list(),
                     cADRnw = list()
  )
  
  #make the subgraph for ease of managing 
  sg = pair_subgraph(g, d1, d2)
  
  result_list$sg = sg
  DRUGS = V(sg)[type=="DRUG"]
  namesDRUGS = V(sg)[type=="DRUG"]$name
  #print(namesDRUGS)
  #list of ADR 
  result_list$ADR = lapply(X = DRUGS, FUN = function(x){
    neighborhood(graph = sg, 
                 order = 2, 
                 nodes = x, 
                 mode = "out", 
                 mindist = 2)[[1]]
    
  }
  )
  
  
  #list of GO 
  
  result_list$GO = lapply(X = DRUGS, FUN = function(x){
    neighborhood(graph = sg, 
                 order = 1, 
                 nodes = x, 
                 mode = "out", 
                 mindist = 1)[[1]]
    
  }
  )
  
  #cADR
  
  cADR = V(sg)[intersect(result_list$ADR[[1]], result_list$ADR[[2]])]
  namesADR = names(cADR)
  
  
  result_list$cADR = lapply(X = cADR, FUN = function(adry){
    paths = lapply(X = namesDRUGS, function(droga){
      
      adr_neigh = neighbors(graph = sg, 
                            v = V(sg)[adry], 
                            mode = "in")
      
      dru_neigh = result_list$GO[[droga]]
      
      respar    = intersect(adr_neigh, dru_neigh)
      res       = V(sg)[respar]
      return(res)
      
    })
    mode_1 = V(sg)[intersect(paths[[1]], paths[[2]])] #if you are counting, count it twice 
    allpaths = V(sg)[union(paths[[1]], paths[[2]])]
    mode_2 = setdiff(allpaths, mode_1)
    mode_2 = V(sg)[mode_2]
    adry_list = list(paths = paths,
                     allpaths = allpaths,
                     mode_1 = mode_1, 
                     mode_2 = mode_2)
    return(adry_list)
  })
  
  ##name each result_list$cADR$"ADR"$paths with drug names 
  for(i in seq_along(result_list$cADR)){
    names(result_list$cADR[[i]]$paths) = namesDRUGS
  }
  
  #make subgraphs 
  
  result_list$cADRnw$nw = composite_subgraph(sg)
  
  #CRMs 
  result_list$cADRnw$components = find_CRM(result_list$cADRnw$nw, precomposite = TRUE)
  result_list$cADRnw$components = components(result_list$cADRnw$components)
  #return final list
  return(result_list)
}



tripartite_full_function = function(g){
  #takes a tripartite graph with DRUGS, GO, and ADRs
  #returns a big list
  
  #get drugs 
  all_drugs = V(g)[type=="DRUG"]
  drug_seq  = seq_along(all_drugs)
  names(drug_seq) = names(all_drugs) #so the lists are named already
  
  #iterate over drugs
  result_list = lapply(X = drug_seq, function(i){
    lapply(X = drug_seq, function(j){
      if(j<=i){return(NA)}else{#this is to avoid doing the calcs twice
        d1 = all_drugs[i]
        d2 = all_drugs[j]
        r  = tripartite_pair_function(g = g, d1 = d1, d2 = d2)
        return(r)
      }
    })
  })
  
  return(result_list)
}


