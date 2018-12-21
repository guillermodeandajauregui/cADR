#heatmap no cADR
hm.cADR.df = getNo.cADR.df(results.gkb)
hm.cADR.df = symmetrizer_function(hm.cADR.df)
hm.cADR.df.clustering = clustering_analysis(hm.cADR.df)

p = ggplot(hm.cADR.df.clustering$ordered.df, aes(drug1, drug2, fill = V1)) +
  geom_raster() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) +
  theme(axis.text.y = element_text(size = 5)) +
  labs(fill = "cADR", x = "", y = "") + 
  scale_fill_gradient2(low = "white", 
                       mid = "turquoise", 
                       high = "tomato4", 
                       midpoint = 40, 
                       limits = c(0, max(hm.cADR.df$V1, 
                                         na.rm = TRUE)
                       )
  ) +
  labs(title = "cADR", subtitle = "")
p
ggsave(filename = "results/hm_cadr.pdf", dpi = 600, width = 240, height = 240, units = "mm")
#mode 1

hm.mode1.df = getNo.mode1.df(results.gkb)
hm.mode1.df = symmetrizer_function(hm.mode1.df)
hm.mode1.df.clustering = clustering_analysis(hm.mode1.df)

p = ggplot(hm.mode1.df.clustering$ordered.df, aes(drug1, drug2, fill = V1)) +
  geom_raster() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) +
  theme(axis.text.y = element_text(size = 5)) +
  labs(fill = "mode_1", x = "", y = "") + 
  scale_fill_gradient2(low = "white", 
                       mid = "skyblue", 
                       high = "sienna1", 
                       midpoint = 40, 
                       limits = c(0, max(hm.mode1.df$V1, 
                                         na.rm = TRUE)
                       )
  ) +
  labs(title = "mode_1", subtitle = "")
p
ggsave(filename = "results/hm_mode1.pdf", dpi = 600, width = 240, height = 240, units = "mm")
#mode 2

hm.mode2.df = getNo.mode2.df(results.gkb)
hm.mode2.df = symmetrizer_function(hm.mode2.df)
hm.mode2.df.clustering = clustering_analysis(hm.mode2.df)

p = ggplot(hm.mode2.df.clustering$ordered.df, aes(drug1, drug2, fill = V1)) +
  geom_raster() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) +
  theme(axis.text.y = element_text(size = 5)) +
  labs(fill = "mode_2", x = "", y = "") + 
  scale_fill_gradient2(low = "white", 
                       mid = "springgreen", 
                       high = "violet", 
                       midpoint = 40, 
                       limits = c(0, max(hm.mode2.df$V1, 
                                         na.rm = TRUE)
                       )
  ) +
  labs(title = "mode_2", subtitle = "")
p
ggsave(filename = "results/hm_mode2.pdf", dpi = 600, width = 240, height = 240, units = "mm")
#CRM

hm.CRM.df = getNo.CRM.df(results.gkb)
hm.CRM.df = symmetrizer_function(hm.CRM.df)
hm.CRM.df.clustering = clustering_analysis(hm.CRM.df)

p = ggplot(hm.CRM.df.clustering$ordered.df, aes(drug1, drug2, fill = V1)) +
  geom_raster() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) +
  theme(axis.text.y = element_text(size = 5)) +
  labs(fill = "CRM", x = "", y = "") + 
  scale_fill_gradient2(low = "white", 
                       mid = "seagreen1", 
                       high = "plum", 
                       midpoint = quantile(hm.CRM.df$V1, 
                                            0.4, 
                                            na.rm = TRUE),
                       #midpoint = 1, 
                       limits = c(0, max(hm.CRM.df$V1, 
                                         na.rm = TRUE)
                       )
  ) +
  labs(title = "CRM", subtitle = "")
p
ggsave(filename = "results/hm_CRM.pdf", dpi = 600, width = 240, height = 240, units = "mm")
