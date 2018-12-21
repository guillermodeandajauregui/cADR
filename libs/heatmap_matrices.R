hm.cADR.df = getNo.cADR.df(results.gkb)
hm.mode1.df = getNo.mode1.df(results.gkb)
hm.mode2.df = getNo.mode2.df(results.gkb)
hm.CRM.df = getNo.CRM.df(results.gkb)

hm.cADR.mx = reshape2::acast(data = hm.cADR.df, formula = drug1~drug2, value.var = "V1")
hm.mode1.mx = reshape2::acast(data = hm.mode1.df, formula = drug1~drug2, value.var = "V1")
hm.mode2.mx = reshape2::acast(data = hm.mode2.df, formula = drug1~drug2, value.var = "V1")
hm.CRM.mx = reshape2::acast(data = hm.CRM.df, formula = drug1~drug2, value.var = "V1")

write.table(x = hm.cADR.mx, 
            file = "results/cADR_matrix.txt", 
            quote = FALSE, 
            sep = "\t", 
            row.names = TRUE, 
            col.names = NA)

write.table(x = hm.mode1.mx, 
            file = "results/mode1_matrix.txt", 
            quote = FALSE, 
            sep = "\t", 
            row.names = TRUE, 
            col.names = NA)

write.table(x = hm.mode2.mx, 
            file = "results/mode2_matrix.txt", 
            quote = FALSE, 
            sep = "\t", 
            row.names = TRUE, 
            col.names = NA)

write.table(x = hm.CRM.mx, 
            file = "results/CRM_matrix.txt", 
            quote = FALSE, 
            sep = "\t", 
            row.names = TRUE, 
            col.names = NA)
