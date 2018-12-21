max(hm.cADR.df$V1, na.rm = TRUE)
hm.cADR.df[which(hm.cADR.df$V1>58),]

unique(grep("sartan", hm.cADR.df$drug1, value = TRUE))
unique(grep("pril", hm.cADR.df$drug1, value = TRUE))
unique(grep("carba", hm.cADR.df$drug1, value = TRUE))

fgkb.twoside = fread(input = "data/3003377s-twosides.tsv")
dict.adrs    = fread(input = "data/diccionario_srep36325-s2.csv")
filter(fgkb.twoside, drug1=="captopril", confidence==5)%>%head

CarbaCapto = results.gkb$captopril$carbamazepine
unique(c(names(CarbaCapto$ADR$captopril), names(CarbaCapto$ADR$carbamazepine)))
dict.adrs[dict.adrs$CUI%in%names(CarbaCapto$cADR),]
dict.adrs[dict.adrs$CUI%in%unique(c(names(CarbaCapto$ADR$captopril), names(CarbaCapto$ADR$carbamazepine))),]

FluoxPhen =  results.gkb$fluoxetine$phenelzine
sort(FluoxPhen$ADR$fluoxetine)
sort(FluoxPhen$ADR$phenelzine)
dict.adrs[dict.adrs$CUI%in%names(FluoxPhen$cADR),]

my_xampl = FluoxPhen$cADRnw$nw
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
write.graph(ready2write, "results/example_FluoxPhen_cadr.gml", "gml")

responseChem = c("C0043352",
"C0018524",
"C0423791",
"C0020443",
"C0004238",
"C0020458",
"C0037274")
dict.adrs[dict.adrs$CUI%in%responseChem,]


length(unique(FluoxPhen$cADR))

filter(hm.mode1.df, drug1 == "fluoxetine", drug2 == "phenelzine")
filter(hm.mode2.df, drug1 == "fluoxetine", drug2 == "phenelzine")
