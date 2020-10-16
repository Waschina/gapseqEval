library(ggplot2)
library(dplyr)
library(data.table)

#
# Comparison
# 

# (0) get table of test organisms
ft.orgs <- fread("organisms2.csv")

# (1) perform tests and merge data.tables
source("getFermProd_gapseq.R")
source("getFermProd_carveme.R")
source("getFermProd_modelseed.R")
dt <- rbindlist(list(gs.fermprod, cm.fermprod, ms.fermprod))

# (2) matching metabolite IDs
met.match <- data.table(matrix(c("acetic acid", "EX_ac(e)", "EX_cpd00029_e0",
                                 "formic acid", "EX_for(e)", "EX_cpd00047_e0",
                                 "ethanol","EX_etoh(e)","EX_cpd00363_e0",
                                 "butyric acid","EX_but(e)","EX_cpd00211_e0",
                                 "propionic acid","EX_ppa(e)","EX_cpd00141_e0",
                                 #"1,2-Propanediol","","EX_cpd00453_e0",
                                 #"(R)-1,2-Propanediol","","EX_cpd01861_e0",
                                 #"1,3-Propanediol-e0","","EX_cpd01618_e0",
                                 "n-butanol","EX_btoh(e)","EX_cpd03662_e0",
                                 #"propanal","EX_ppal(e)","EX_cpd00371_e0",
                                 "L-lactic acid","EX_lac__L(e)","EX_cpd00159_e0",
                                 "D-lactic acid","EX_lac__D(e)","EX_cpd00221_e0",
                                 "succinic acid","EX_succ(e)","EX_cpd00036_e0",
                                 "acetone","EX_acetone(e)","EX_cpd00178_e0",
                                 "H+", "EX_h(e)", "EX_cpd00067_e0",
                                 "methane", "EX_ch4(e)","EX_cpd01024_e0",
                                 "H2", "EX_h2(e)", "EX_cpd11640_e0"
), ncol = 3, byrow = T))
colnames(met.match) <- c("metabolite","bigg","seed")
dt <- merge(dt, met.match[,-"bigg"], by.x = "ex", by.y = "seed", all.x = T)
dt <- merge(dt, met.match[,-"seed"], by.x = "ex", by.y = "bigg", all.x = T, suffixes = c("",".y"))
dt[!is.na(metabolite.y), metabolite := metabolite.y]
dt[,metabolite.y := NULL]
dt <- dt[!is.na(metabolite)]

# (3) fill missing values with 0 and combine both lactate enantiomers to DL-lactic acid
#dt <- dcast(dt, .id + recon.method + gr.rate ~ metabolite, value.var = c("l","u","mtf.flux"))
dt[metabolite %in% c("L-lactic acid","D-lactic acid"), u := max(u), by = c(".id", "recon.method")]
dt[metabolite %in% c("L-lactic acid","D-lactic acid"), mtf.flux := max(mtf.flux), by = c(".id", "recon.method")]
dt[metabolite == "L-lactic acid", metabolite := "DL-lactic acid"]
dt <- dt[metabolite != "D-lactic acid"]
dt <- dt[metabolite != "H+"]
# fill missing 0 value entries
completeDT <- function(DT, cols, defs = NULL){
  mDT = do.call(CJ, c(DT[, ..cols], list(unique=TRUE)))
  res = DT[mDT, on=names(mDT)]
  if (length(defs)) 
    res[, names(defs) := Map(replace, .SD, lapply(.SD, is.na), defs), .SDcols=names(defs)]
  res[]
} 
dt <- completeDT(dt, cols = c(".id","metabolite","recon.method"), defs = c(l = 0, u = 0, mtf.flux = 0))

# (4) get human-readible organism names
dt <- merge(dt, ft.orgs, by.x = ".id", by.y = "id")
dt[, exp.measured := grepl(metabolite, fmets), by = c(".id","metabolite","recon.method")]
dt[l <= 0, l := NA]
dt[mtf.flux <= 1e-5, mtf.flux := NA]

#dt[mtf.flux/gr.rate > 600, mtf.flux := 0.042705]
#dt[u/gr.rate > 500, u := 0.042705]

# Remove Archaea from analysis since this is still experimental, though, it looks good with gapseq.
dt <- dt[organism != "Methanobrevibacter smithii ATCC 35061"]
dt <- dt[organism != "Methanosarcina barkeri str. Fusaro"]
dt <- dt[metabolite != "methane"]

dt[mtf.flux / gr.rate < 0.25, mtf.flux := NA]
dt[u / gr.rate < 0.25 | is.na(ex), u := NA]

# some aesthetics 
dt[, organism := gsub("\\[|\\]","",organism)]
dt$organism <- factor(dt$organism,levels = sort(unique(dt$organism), decreasing = T))
#levels(dt$organism) <- sort(unique(dt$organism), decreasing = T)

dt[metabolite == "acetic acid", metabolite := "Acetate"]
dt[metabolite == "butyric acid", metabolite := "Butyrate"]
dt[metabolite == "DL-lactic acid", metabolite := "Lactate"]
dt[metabolite == "ethanol", metabolite := "Ethanol"]
dt[metabolite == "formic acid", metabolite := "Formate"]
dt[metabolite == "propionic acid", metabolite := "Propionate"]
dt[metabolite == "succinic acid", metabolite := "Succinate"]

p <- ggplot(dt, aes(metabolite, organism)) +
  geom_tile(aes(fill = exp.measured), color = "white", size = 1.25) +
  scale_fill_manual(values = c("white","#2ca25f")) +
  geom_point(aes(size = u / gr.rate), color = "darkgrey", alpha = 1, pch = 16) +
  geom_point(aes(size = mtf.flux / gr.rate), color = "black") +
  #geom_point(aes(size = l / gr.rate), col = "white", alpha = "0.8") +
  scale_size(range = c(0, 8.5)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, color = "black"),
        axis.text.y=element_text(color = "black")) +
  facet_grid(.~recon.method)
p  

ggsave("plots/ferm.prod.comparison201015.eps", p, width = 14, height = 8.5)

# Barplots for TP, FP, TN, FN
table(dt$exp.measured, !is.na(dt$u), dt$recon.method)
table(dt$exp.measured, !is.na(dt$mtf.flux), dt$recon.method)
dt.lt <- as.data.table(table(dt$exp.measured, !is.na(dt$mtf.flux), dt$recon.method))
dt.lt[V1 == T & V2 == F, catg := "FN"]
dt.lt[V1 == F & V2 == T, catg := "FP"]
dt.lt[V1 == T & V2 == T, catg := "TP"]
dt.lt[V1 == F & V2 == F, catg := "TN"]

p <- ggplot(dt.lt, aes(fill=catg, y=N, x=V3)) + 
  geom_bar(position="stack", stat="identity", col = "black") + 
  coord_flip() + scale_fill_manual(values = c("#e34a33","#fdbb84","#2b8cbe","#a6bddb"))
p

ggsave("plots/ferm.prod.logicalTest201015.eps", p, width = 10, height = 2.5)

