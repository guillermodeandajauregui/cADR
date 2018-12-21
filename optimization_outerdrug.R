tripartite_outer_list_function = function(g, drug_out, in_list){
  #takes a graph with DRUGS, GO. and ADRs
  #and two drugs 
  #returns a big list of results 
  #this will take a drug as the "outer" drug and another as the "inner" drug
  #to try and make the code more efficient, reducing the number of subsettings 
  
  
  result_list = list(ADR = list(), 
                     GO  = list(),
                     cADR   = list(),
                     cADRnw = list()
  )
  
  #extract neighborhoods of the outer drug 
  
  #make a subgraph 
  g.drug.out = make_ego_graph(g, 
                              order = 2, 
                              mode = "out", 
                              nodes = drug_out)[[1]]
  
  #get outer drug list of GOs (pathways)
  outer.pathways = names(V(g.drug.out)[type == "GO"])
  
  #get outer drug list of ADR
  outer.adr  = names(V(g.drug.out)[type == "ADR"])
  
  #######################################
  #Now, take the list of inner drugs 
  
  lapply(X = in_list, FUN = function(drug_in){
    #get drug_in subgraph
    g.drug.in = make_ego_graph(g, 
                                order = 2, 
                                mode = "out", 
                                nodes = drug_in)[[1]]
    
    #get in.drug list of GOs (pathways)
    inner.pathways = names(V(g.drug.in)[type == "GO"])
    
    #get in.drug list of ADR
    inner.adr  = names(V(g.drug.in)[type == "ADR"])
    
    #composite ADR 
    cADR = intersect(outer.adr, inner.adr)
    
    #paths to cADR 
    outer.paths = lapply(X = cADR, FUN = function(cadri){
      neighbors(graph = g.drug.out, v = cadri, mode = "in")
    })
    inner.paths = lapply(X = cADR, FUN = function(cadri){
      neighbors(graph = g.drug.in, v = cadri, mode = "in")
    })
    
    #mode one pathways 
    mode_1 = lapply(X = seq_along(cADR), FUN = function(i){
      intersect(x = outer.paths[[i]], inner.paths[[i]])
    })
    
    #mode two pathways
    mode_2 = lapply(X = seq_along(cADR), FUN = function(i){
      a = setdiff(x = outer.paths[[i]], inner.paths[[i]])
      b = setdiff(x = inner.paths[[i]], outer.paths[[i]])
      return(union(a,b))
    })
    
    
  
  
  #this is the parenthesis that closes the lapply
  })
  #######################################
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