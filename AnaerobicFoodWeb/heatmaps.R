library(data.table)
library(ggplot2)
library(dplyr)
my_breaks <- c(0.01,0.1,1,10,100,1000)

# gapseq
a.gs <- readRDS("dat/ex.gapseq.RDS")
a.gs$recon <- "gapseq"
a.gs <- a.gs[time == 3]

# carveme
a.cm <- readRDS("dat/ex.carveme.RDS")
a.cm <- a.cm[time == 3]
a.cm.tmp <- copy(a.cm[grepl("lac__", sub)])
a.cm.tmp[, mflux := sum(mflux), by = spec.name]
a.cm.tmp <- a.cm.tmp[!grepl("lac__D", sub)]
a.cm <- rbind(a.cm[!grepl("lac__", sub)], a.cm.tmp)
a.cm[is.na(stat) & mflux < 0, stat := "uptake"]
a.cm[is.na(stat) & mflux > 0, stat := "production"]
a.cm$recon <- "CarveMe"

# modelseed
a.ms <- readRDS("dat/ex.modelseed.RDS")
a.ms$recon <- "ModelSEED"
a.ms <- a.ms[time == 5]

# combine
a <- rbind(a.gs, rbind(a.cm, a.ms))

# aestetics (names)
a[, spec.name := gsub("^ ","", spec.name)]
a[, spec.name := gsub(" chromosome| DNA| genome assembly| E_hallii.*$| Scfld.*$","",spec.name)]
a[, sub := gsub("-e0$","",sub)]

a <- a[abs(mflux)>0]
#a[stat == "uptake", mflux := mflux*100]
a[, mflux2 := abs(mflux)]

# repair modelseed sub names
a[recon == "ModelSEED" & sub == "methane", sub := "Methane"]
a[recon == "ModelSEED" & sub == "Hydrogen sulfide", sub := "H2S"]

# repaire carveme sub names
a[recon == "CarveMe" & sub == "EX_h2s(e)", sub := "H2S"]
a[recon == "CarveMe" & sub == "EX_ac(e)", sub := "Acetate"]
a[recon == "CarveMe" & sub == "EX_succ(e)", sub := "Succinate"]
a[recon == "CarveMe" & sub == "EX_for(e)", sub := "Formate"]
a[recon == "CarveMe" & sub == "EX_etoh(e)", sub := "Ethanol"]
a[recon == "CarveMe" & sub == "EX_h2(e)", sub := "H2"]
a[recon == "CarveMe" & sub == "EX_ch4(e)", sub := "Methane"]
a[recon == "CarveMe" & sub == "EX_lac__L(e)", sub := "L-Lactate"]
#a[recon == "CarveMe" & sub == "EX_h2s(e)", sub := "H2S"]

# plotting ( High: )
sub.incl <- "H2S|Acetate|Butyrate|Ethanol|Formate|H2|actate|Methane|Propio|Succ"
p <- ggplot(a[grepl(sub.incl, sub)], aes(sub, reorder(spec.name, desc(spec.name)))) + 
  #geom_tile(data = tmp[!grepl(sub.excl, V1)], aes(V1, V2, fill = N)) +
  geom_tile(aes(fill=mflux2), colour = "white", size = 1.5) + scale_fill_gradient(low = "white", high = "red", name = "Production\n [unit]", trans = "log",
                                                                                  breaks = my_breaks, labels = my_breaks) +
  theme_bw() +
  facet_grid(stat~recon) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                 panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
ggsave(filename = "plots/comm_heat_gg_red.pdf", plot = p, width = 13.25, height = 11)

p <- ggplot(a[grepl(sub.incl, sub)], aes(sub, reorder(spec.name, desc(spec.name)))) + 
  #geom_tile(data = tmp[!grepl(sub.excl, V1)], aes(V1, V2, fill = N)) +
  geom_tile(aes(fill=mflux2), colour = "white", size = 1.5) + scale_fill_gradient(low = "white", high = "blue", name = "Production\n [unit]", trans = "log",
                                                                                  breaks = my_breaks, labels = my_breaks) +
  theme_bw() +
  facet_grid(stat~recon) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                 panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
ggsave(filename = "plots/comm_heat_gg_blue.pdf", plot = p, width = 13.25, height = 11)


