#example 1: 
V(g.tri)[type=="GO"][262]
x = as.numeric(neighbors(g.tri, "GO:0005249", "all"))
y = as.numeric(V(g.tri)[type=="GO"][262])

red_test = induced_subgraph(g.tri, vids = V(g.tri)[c(x,y)])

#example 2 
test2= make_ego_graph(g.tri, 
                      order = 2, 
                      nodes = "C0037199", #sinusitis
                      mode = "in")[[1]]

tt   = pair_subgraph(g.tri, "BRD-K13078532", "BRD-K72222507") #hydrochlorothiazide y quinapril

unlist(neighbors(tt, "C0037199", "in"))

results.gkb$hydroxychloroquine$omeprazole$cADR$C0085593
results.gkb$hydroxychloroquine$omeprazole$cADR$C0043352

x = results.gkb$hydroxychloroquine$omeprazole$cADR
lapply(X = x, FUN = function(x){
  r = x[["allpaths"]]
})

results.gkb$hydroxychloroquine$omeprazole$cADR$C0035455$mode_1
results.gkb$hydroxychloroquine$omeprazole$cADR$C0035455$mode_2

tt = results.gkb$hydroxychloroquine$omeprazole$sg
happyPlot(tt)

ttt = neighborhood(graph = tt, order = 1, nodes = c("GO:0072350", "GO:1903293", "GO:0008635"))
ttt = unique(unlist(ttt))
ttt = induced.subgraph(graph = tt, vids = ttt)
happyPlot(ttt)

saveRDS(ttt, file = "ejemplo_redADR.RData")
