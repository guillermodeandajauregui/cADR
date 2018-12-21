library(igraph)
library(tidyverse)
library(data.table)
library(plyr)

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

remove_headless_go = function(g, v = V(g)[type=="GO"]){
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



tripartite_full_function = function(g, kores_o = 3, kores_i = 3){
  #takes a tripartite graph with DRUGS, GO, and ADRs
  #returns a big list
  
  #get drugs 
  all_drugs = V(g)[type=="DRUG"]
  drug_seq  = seq_along(all_drugs)
  names(drug_seq) = names(all_drugs) #so the lists are named already
  
  #iterate over drugs
  result_list = parallel::mclapply(X = drug_seq, mc.cores = kores_o, function(i){
    parallel::mclapply(X = drug_seq, mc.cores = kores_i, function(j){
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

#result analysis functions

#extract values from result list
##example, with number of cADR

number.cadr = function(result_list){
  lapply(X = result_list, FUN = function(i){
    lapply(X = i, FUN = function(j){
      if(is.na(j)){return(NA)}else{
        return(length(j$cADR))
      }
    })
  })  
}

#take a simplified (Extracted) result list, return a matrix (for plotting)

resultList2Matrix = function(result_list){
  #list must have value L2 L1 order AND column names
  
  #melt to dataframe
  t3 = melt(result_list)
  #make intermediate mirror data frame
  tm = t3[c(1,3,2)]
  #and change colnames to make them rbindable
  colnames(tm) = c("value", "L2", "L1")
  #remove the rows where L1 == L2 in the intermidiate mirror
  tm = tm[-(which(tm$L1==tm$L2)),]
  #rbind
  t5 = rbind(t3, tm)
  #Remove the rows with NA values
  t5 = t5[-which(t5$L1 != t5$L2 & is.na(t5$value)),]
  t6 = dcast(t5, L1~L2)
  t7 = as.matrix(t6[,-1])
  rownames(t7) = t6$L1
  return(t7)
}

reorder.matrix2df = function(mx){
  ord = hclust(dist(mx))$order
  t6 = as.matrix(mx)[ord,ord]
  t7 = reshape2::melt(data = t6, na.rm = FALSE)
  return(t7)
}

#take a list of results and return a data frame
getNo.cADR.df = function(resultList){
  r1 = lapply(resultList, FUN = function(z){
    ldply(z, .id = "drug2", .fun = function(y){
      if(length(y)!=5){return(NA)} #full results always have 5 
      else{return(length(y[["cADR"]]
      )
      )
      }
    })  
  })
  
  r2 = dplyr::bind_rows(r1, .id = 'drug1')
  return(r2)
}

#mode 1
getNo.mode1.df = function(resultList){
  r1 = lapply(resultList, FUN = function(z){
    ldply(z, .id = "drug2", .fun = function(y){
      if(length(y)!=5){return(NA)} #full results always have 5 
      else{
        w = y[["cADR"]]
        r = sum(sapply(X = w, FUN = function(zz){
          length(zz$mode_1)
        }))
        return(r)
      }})  
  })
  
  r2 = dplyr::bind_rows(r1, .id = 'drug1')
  return(r2)
}

#mode 2
getNo.mode2.df = function(resultList){
  r1 = lapply(resultList, FUN = function(z){
    ldply(z, .id = "drug2", .fun = function(y){
      if(length(y)!=5){return(NA)} #full results always have 5 
      else{
        w = y[["cADR"]]
        r = sum(sapply(X = w, FUN = function(zz){
          length(zz$mode_2)
        }))
        return(r)
      }})  
  })
  
  r2 = dplyr::bind_rows(r1, .id = 'drug1')
  return(r2)
}

#no CRM

getNo.CRM.df = function(resultList){
  r1 = lapply(resultList, FUN = function(z){
    ldply(z, .id = "drug2", .fun = function(y){
      if(length(y)!=5){return(NA)} #full results always have 5 
      else{
        w = y[["cADRnw"]]$components$no
        return(w)
      }})  
  })
  
  r2 = dplyr::bind_rows(r1, .id = 'drug1')
  return(r2)
}

#symetrizer function

symmetrizer_function = function(df){
  #takes a ggplot ready 3 column data frame
  # Elem1 Elem2 Value 
  #that has A B X but not B A X
  #and makes it simmetrical in a very naive way
  
  #drop NAs that are NOT the same element
  # A B NA is dropped, but A A NA is kept
  idx2drop = which(df[,1]!=df[,2] & is.na(df[,3]))
  df2      = df[-idx2drop,]
  
  #make a second one, where element 1 and element 2 will be swapped
  df3      = df2[,c(2,1,3)]     
  ##drop the repeated ones in this one
  todrop2  = which(df3[,1]==df3[,2])
  df3      = df3[-todrop2,]
  
  #make colnames the same so we can paste
  colnames(df3) <- colnames(df2)
  
  #paste!
  df4     = rbind(df2, df3)
  
  #and make sure elements are characters, not factors
  
  df4[,1] = as.character(df4[,1])
  df4[,2] = as.character(df4[,2])
  
  return(df4)
}

#clustering analysis
clustering_analysis <- function(df, 
                                clustering.method = "complete", 
                                distance.method = "euclidean",
                                formula = drug1~drug2){
  #takes a ggplot ready 3 column data frame
  # Elem1 Elem2 Value
  #with all possible combination values
  # A B X and B A X must be in it
  #and will return a list 
  ##with an hclust object 
  ##and a copy of the original data frame, but reordered 
  
  #copy the names of the original dataframe
  originalNames = colnames(df)
  
  #cast from data frame
  mx = reshape2::acast(df, formula = formula)
  
  #clustering object 
  clustering.object = hclust(d = dist(x = mx, 
                                      method = distance.method),
                             method = clustering.method
  )
  
  #reorder matrix
  mx = mx[clustering.object$order, clustering.object$order]
  
  #melt matrix
  df2 = melt(mx)
  df2 = df2[,c(2,1,3)] #so they have the same order as before
  colnames(df2) = originalNames
  
  rl = list(clustering.object = clustering.object,
            ordered.df = df2)
  return(rl)
  
}

