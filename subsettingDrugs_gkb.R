
#take the original g.tri

#remove terms with hyphens 
g.tri.sub = grep(pattern = "-", V(g.tri)[type=="DRUG"]$drugName, invert = TRUE)
g.tri.sub = neighborhood(graph = g.tri, 
                         order = 2, 
                         nodes = V(g.tri)[type=="DRUG"][g.tri.sub], 
                         mode = "out")
g.tri.sub = induced.subgraph(graph = g.tri, vids =  unlist(g.tri.sub))

#read TWOSIDES from pharmGKB
twosides = fread(input = "3003377s-twosides.tsv")
drugsGKB = unique(c(unique(twosides$drug1), unique(twosides$drug2)))
nomendrug = V(g.tri.sub)[type=="DRUG"]$drugName
#
nomen_gkb = nomendrug[nomendrug%in%drugsGKB]

#subset to keep only drugs in GKB-TWOSIDES
g.tri.gkb = neighborhood(g.tri, 
                         order = 2, 
                         nodes = V(g.tri)[V(g.tri)$drugName%in%nomen_gkb], 
                         mode = "out")

g.tri.gkb = induced.subgraph(graph = g.tri, vids =  unlist(g.tri.gkb))
table(V(g.tri.gkb)$type)

#rename drug nodes with generic drug name

V(g.tri.gkb)[type=="DRUG"]$drugSymbol = V(g.tri.gkb)[type=="DRUG"]$name
V(g.tri.gkb)[type=="DRUG"]$name = V(g.tri.gkb)[type=="DRUG"]$drugName
V(g.tri.gkb)[type=="DRUG"]$name[1:5]

#Analyze

proctor = proc.time()
results.gkb = tripartite_full_function(g = g.tri.gkb, kores_o = 4, kores_i = 1)
proctor = proc.time() - proctor
proctor

save(results.gkb, file = "results/results_gkb.RData")
