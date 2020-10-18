# wrapper function for all GE Analysis (Gapseq, CarveMe, ModelSEED, Curated models)
source("calc_seq_overlap.R")
source("prepare_MS_gene_table.R")
source("../FermentationTest/sybil_toolkit.R")
source("prediction.stats.R")
source("predGeneEssentiality.GS.R")
source("predGeneEssentiality.CM.R")
source("predGeneEssentiality.MS.R")
library(sybilSBML)
library(data.table)
library(stringr)
library(ggplot2)
library(fmsb) # for spider web plot
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1

#~~~~~~~~~~~~~~~#
# Gapseq Models #
#~~~~~~~~~~~~~~~#
ecol.gs <- readRDS("models/gapseq/ecol.RDS")
ecol.gs <- constrain.mod(ecol.gs, mediadb = "media/media.tsv", media.id = "GS_MM_glc")
paer.gs <- readRDS("models/gapseq/paer.RDS")
paer.gs <- constrain.mod(paer.gs, mediadb = "media/media.tsv", media.id = "GS_MM_succ")
bsub.gs <- readRDS("models/gapseq/bsub.RDS")
bsub.gs <- constrain.mod(bsub.gs, mediadb = "media/media.tsv", media.id = "LB_marinos")
sone.gs <- readRDS("models/gapseq/sone.RDS")
sone.gs <- constrain.mod(sone.gs, mediadb = "media/media.tsv", media.id = "LB_marinos")
mgen.gs <- readRDS("models/gapseq/mgen.RDS")
mgen.gs <- constrain.mod(mgen.gs, mediadb = "media/media.tsv", media.id = "Complete2")

#~~~~~~~~~~~~~~~#
# Carveme Models#
#~~~~~~~~~~~~~~~#
bsub.cm <- readSBMLmod("models/carveme/bsub_protein.xml")
bsub.cm <- constrain.mod(bsub.cm, mediadb = "media/media.tsv", media.id = "LB_marinos")
ecol.cm <- readSBMLmod("models/carveme/ecol_protein.xml")
ecol.cm <- constrain.mod(ecol.cm, mediadb = "media/media.tsv", media.id = "GS_MM_glc")
paer.cm <- readSBMLmod("models/carveme/paer_protein.xml")
paer.cm <- constrain.mod(paer.cm, mediadb = "media/media.tsv", media.id = "GS_MM_succ")
mgen.cm <- readSBMLmod("models/carveme/mgen_protein.xml")
mgen.cm <- constrain.mod(mgen.cm, mediadb = "media/media.tsv", media.id = "Complete2")
sone.cm <- readSBMLmod("models/carveme/sone_protein.xml")
sone.cm <- constrain.mod(sone.cm, mediadb = "media/media.tsv", media.id = "LB_marinos")

## The following is a procedure to get the genomic position of the genes encoding the proteins that are the input for CarveMe Reconstruction.
## it needs to be performed only once and followed by the bash procedure saved in file "get_protein_gene_loci.sh"
#getGeneList <- function(mod) {
#  gl <- mod@allGenes
#  gl <- gsub("^G_","",gl)
#  gl <- gl[gl != "spontaneous"]
#  fwrite(data.table(X = gl), col.names = F, file = paste0("models/carveme/",deparse(substitute(mod)),"_genes.lst"))
#}

#~~~~~~~~~~~~~~~~~~#
# ModelSEED Models #
#~~~~~~~~~~~~~~~~~~#
bsub.ms <- readSBMLmod("models/modelseed/bsub.sbml")
bsub.ms <- constrain.mod(bsub.ms, mediadb = "media/media.tsv", media.id = "LB_marinos")
ecol.ms <- readSBMLmod("models/modelseed/ecol.sbml")
ecol.ms <- constrain.mod(ecol.ms, mediadb = "media/media.tsv", media.id = "GS_MM_glc")
paer.ms <- readSBMLmod("models/modelseed/paer.sbml")
paer.ms <- constrain.mod(paer.ms, mediadb = "media/media.tsv", media.id = "GS_MM_succ")
mgen.ms <- readSBMLmod("models/modelseed/mgen.sbml")
#mgen.ms <- constrain.mod(mgen.ms, mediadb = "media/media.tsv", media.id = "Complete2") # modelseed failed to gapfill with gapseq's complete medium. Thus we use the modelseeds complete medium.
sone.ms <- readSBMLmod("models/modelseed/sone.sbml")
sone.ms <- constrain.mod(sone.ms, mediadb = "media/media.tsv", media.id = "LB_marinos")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Perform Gene essentiality analysis #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# ecol
getest.ecol.gs <- predGeneEssentiality.GS(ecol.gs, "essentiality.data/loci-ecol.csv", "essentiality.data/gess-ecol.csv")
getest.ecol.cm <- predGeneEssentiality.CM(ecol.cm, "essentiality.data/loci-ecol.csv", "essentiality.data/gess-ecol.csv", "models/carveme/ecol_genes_pos.tsv")
getest.ecol.ms <- predGeneEssentiality.MS(ecol.ms, "essentiality.data/loci-ecol.csv", "essentiality.data/gess-ecol.csv", "models/modelseed/ecol_genes.gtbl")

# bsub
getest.bsub.gs <- predGeneEssentiality.GS(bsub.gs, "essentiality.data/loci-bsub.csv", "essentiality.data/gess-bsub.csv")
getest.bsub.cm <- predGeneEssentiality.CM(bsub.cm, "essentiality.data/loci-bsub.csv", "essentiality.data/gess-bsub.csv", "models/carveme/bsub_genes_pos.tsv")
getest.bsub.ms <- predGeneEssentiality.MS(bsub.ms, "essentiality.data/loci-bsub.csv", "essentiality.data/gess-bsub.csv", "models/modelseed/bsub_genes.gtbl")

# sone
getest.sone.gs <- predGeneEssentiality.GS(sone.gs, "essentiality.data/loci-sone.csv", "essentiality.data/gess-sone.csv")
getest.sone.cm <- predGeneEssentiality.CM(sone.cm, "essentiality.data/loci-sone.csv", "essentiality.data/gess-sone.csv", "models/carveme/sone_genes_pos.tsv")
getest.sone.ms <- predGeneEssentiality.MS(sone.ms, "essentiality.data/loci-sone.csv", "essentiality.data/gess-sone.csv", "models/modelseed/sone_genes.gtbl")

# paer
getest.paer.gs <- predGeneEssentiality.GS(paer.gs, "essentiality.data/loci-paer.csv", "essentiality.data/gess-paer.csv")
getest.paer.cm <- predGeneEssentiality.CM(paer.cm, "essentiality.data/loci-paer.csv", "essentiality.data/gess-paer.csv", "models/carveme/paer_genes_pos.tsv")
getest.paer.ms <- predGeneEssentiality.MS(paer.ms, "essentiality.data/loci-paer.csv", "essentiality.data/gess-paer.csv", "models/modelseed/paer_genes.gtbl")

# mgen
getest.mgen.gs <- predGeneEssentiality.GS(mgen.gs, "essentiality.data/loci-mgen.csv", "essentiality.data/gess-mgen.csv")
getest.mgen.cm <- predGeneEssentiality.CM(mgen.cm, "essentiality.data/loci-mgen.csv", "essentiality.data/gess-mgen.csv", "models/carveme/mgen_genes_pos.tsv")
getest.mgen.ms <- predGeneEssentiality.MS(mgen.ms, "essentiality.data/loci-mgen.csv", "essentiality.data/gess-mgen.csv", "models/modelseed/mgen_genes.gtbl")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plotting the results               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
gess_spider_web <- function(ge.gs, ge.cm, ge.ms) {
  dt <- merge(ge.gs$loci[,.(locus, gene2, start, end, GSmod.genes = mod.genes, ess.curated.model, ess.experimental, ess.gapseq)],
              ge.cm$loci[,.(locus, CMmod.genes = mod.genes, ess.carveme)],
              by = "locus")
  dt <- merge(dt,
              ge.ms$loci[,.(locus, CMmod.genes = mod.genes, ess.modelse)],
              by = "locus")
  
  tmp <- dt[!is.na(ess.experimental) & !is.na(ess.curated.model) & !is.na(ess.gapseq) & !is.na(ess.carveme) & !is.na(ess.modelse)]
  #tmp <- dt[!is.na(ess.experimental) & !is.na(ess.curated.model)]
  #tmp <- dt[!is.na(ess.experimental)]
  
  tmp.cur <- table(tmp$ess.experimental, tmp$ess.curated.model)
  tmp.gs  <- table(tmp$ess.experimental, tmp$ess.gapseq)
  tmp.cm  <- table(tmp$ess.experimental, tmp$ess.carveme)
  tmp.ms  <- table(tmp$ess.experimental, tmp$ess.modelse)
  print(tmp.gs)
  print(tmp.cm)
  print(tmp.ms)
  
  eval.cur <- prediction.stats(tmp.cur)
  eval.gs  <- prediction.stats(tmp.gs)
  eval.cm  <- prediction.stats(tmp.cm)
  eval.ms  <- prediction.stats(tmp.ms)
  
  eval.mat <- matrix(0, ncol = 5, nrow = 4)
  eval.mat[1,] <- eval.cur
  eval.mat[2,] <- eval.gs
  eval.mat[3,] <- eval.cm
  eval.mat[4,] <- eval.ms
  
  eval.lim <- matrix(c(rep(1,5),rep(0,5)), ncol = 5, byrow = T)
  
  eval.mat <- rbind(eval.lim, eval.mat)
  colnames(eval.mat) <- names(eval.cur)
  colnames(eval.mat)[1] <- "F1-score"
  rownames(eval.mat) <- c("max","min","Curated","GapSeq","CarveMe","ModelSEED")
  
  library(fmsb)
  print(as.data.frame(eval.mat))
  radarchart(as.data.frame(eval.mat), axistype=0,
             #custom the grid
             cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,0.25,4), cglwd=0.8, centerzero = T,
             plty = c(1,1,1), pcol = c("#000000", "#e41a1c", "#377eb8", "#FFB200"), plwd = 2, pty = 32,
             title = paste("Nr. of common loci:", nrow(tmp)))
}

pdf("plots/gene.ess.spiderwebs_co0_01.pdf", width = 5.25, height = 5.25)
# ecol
gess_spider_web(getest.ecol.gs, getest.ecol.cm, getest.ecol.ms)
# bsub
gess_spider_web(getest.bsub.gs, getest.bsub.cm, getest.bsub.ms)
# sone
gess_spider_web(getest.sone.gs, getest.sone.cm, getest.sone.ms)
# paer
gess_spider_web(getest.paer.gs, getest.paer.cm, getest.paer.ms)
# mgen
gess_spider_web(getest.mgen.gs, getest.mgen.cm, getest.mgen.ms)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# basic model stat comparison (nr of genes, reactions, metabolites) #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
source("gs_getGeneNr.R")
ecol <- readSBMLmod("models/curated/iML1515.xml")
bsub <- readSBMLmod("models/curated/iYO844.xml")
paer <- readSBMLmod("models/curated/iMO1056.xml")
sone <- readSBMLmod("models/curated/iSO783.xml")
mgen <- readSBMLmod("models/curated/iPS189_v4.xml")


dt.stat <- data.table(organism = rep(NA_character_, 60), recon.method = NA_character_, metric = NA_character_, nr = NA_integer_)
# GENES
dt.stat[1, `:=`(organism = "E. coli"      , recon.method = "curated", metric = "Genes", nr = length(ecol@allGenes))]
dt.stat[2, `:=`(organism = "B. subtilis"  , recon.method = "curated", metric = "Genes", nr = length(bsub@allGenes))]
dt.stat[3, `:=`(organism = "P. aeruginosa", recon.method = "curated", metric = "Genes", nr = length(paer@allGenes))]
dt.stat[4, `:=`(organism = "S. oneidensis", recon.method = "curated", metric = "Genes", nr = length(sone@allGenes))]
dt.stat[5, `:=`(organism = "M. genitalium", recon.method = "curated", metric = "Genes", nr = length(mgen@allGenes))]

dt.stat[6, `:=`(organism = "E. coli"      , recon.method = "CarveMe", metric = "Genes", nr = length(ecol.cm@allGenes))]
dt.stat[7, `:=`(organism = "B. subtilis"  , recon.method = "CarveMe", metric = "Genes", nr = length(bsub.cm@allGenes))]
dt.stat[8, `:=`(organism = "P. aeruginosa", recon.method = "CarveMe", metric = "Genes", nr = length(paer.cm@allGenes))]
dt.stat[9, `:=`(organism = "S. oneidensis", recon.method = "CarveMe", metric = "Genes", nr = length(sone.cm@allGenes))]
dt.stat[10,`:=`(organism = "M. genitalium", recon.method = "CarveMe", metric = "Genes", nr = length(mgen.cm@allGenes))]

dt.stat[11,`:=`(organism = "E. coli"      , recon.method = "gapseq", metric = "Genes", nr = gs_getGeneNr(ecol.gs))]
dt.stat[12,`:=`(organism = "B. subtilis"  , recon.method = "gapseq", metric = "Genes", nr = gs_getGeneNr(bsub.gs))]
dt.stat[13,`:=`(organism = "P. aeruginosa", recon.method = "gapseq", metric = "Genes", nr = gs_getGeneNr(paer.gs))]
dt.stat[14,`:=`(organism = "S. oneidensis", recon.method = "gapseq", metric = "Genes", nr = gs_getGeneNr(sone.gs))]
dt.stat[15,`:=`(organism = "M. genitalium", recon.method = "gapseq", metric = "Genes", nr = gs_getGeneNr(mgen.gs))]

dt.stat[16,`:=`(organism = "E. coli"      , recon.method = "modelSEED", metric = "Genes", nr = length(ecol.ms@allGenes))]
dt.stat[17,`:=`(organism = "B. subtilis"  , recon.method = "modelSEED", metric = "Genes", nr = length(bsub.ms@allGenes))]
dt.stat[18,`:=`(organism = "P. aeruginosa", recon.method = "modelSEED", metric = "Genes", nr = length(paer.ms@allGenes))]
dt.stat[19,`:=`(organism = "S. oneidensis", recon.method = "modelSEED", metric = "Genes", nr = length(sone.ms@allGenes))]
dt.stat[20,`:=`(organism = "M. genitalium", recon.method = "modelSEED", metric = "Genes", nr = length(mgen.ms@allGenes))]

# Reactions
dt.stat[21, `:=`(organism = "E. coli"      , recon.method = "curated", metric = "Reactions", nr = ecol@react_num)]
dt.stat[22, `:=`(organism = "B. subtilis"  , recon.method = "curated", metric = "Reactions", nr = bsub@react_num)]
dt.stat[23, `:=`(organism = "P. aeruginosa", recon.method = "curated", metric = "Reactions", nr = paer@react_num)]
dt.stat[24, `:=`(organism = "S. oneidensis", recon.method = "curated", metric = "Reactions", nr = sone@react_num)]
dt.stat[25, `:=`(organism = "M. genitalium", recon.method = "curated", metric = "Reactions", nr = mgen@react_num)]

dt.stat[26, `:=`(organism = "E. coli"      , recon.method = "CarveMe", metric = "Reactions", nr = ecol.cm@react_num)]
dt.stat[27, `:=`(organism = "B. subtilis"  , recon.method = "CarveMe", metric = "Reactions", nr = bsub.cm@react_num)]
dt.stat[28, `:=`(organism = "P. aeruginosa", recon.method = "CarveMe", metric = "Reactions", nr = paer.cm@react_num)]
dt.stat[29, `:=`(organism = "S. oneidensis", recon.method = "CarveMe", metric = "Reactions", nr = sone.cm@react_num)]
dt.stat[30, `:=`(organism = "M. genitalium", recon.method = "CarveMe", metric = "Reactions", nr = mgen.cm@react_num)]

dt.stat[31, `:=`(organism = "E. coli"      , recon.method = "gapseq", metric = "Reactions", nr = ecol.gs@react_num)]
dt.stat[32, `:=`(organism = "B. subtilis"  , recon.method = "gapseq", metric = "Reactions", nr = bsub.gs@react_num)]
dt.stat[33, `:=`(organism = "P. aeruginosa", recon.method = "gapseq", metric = "Reactions", nr = paer.gs@react_num)]
dt.stat[34, `:=`(organism = "S. oneidensis", recon.method = "gapseq", metric = "Reactions", nr = sone.gs@react_num)]
dt.stat[35, `:=`(organism = "M. genitalium", recon.method = "gapseq", metric = "Reactions", nr = mgen.gs@react_num)]

dt.stat[36, `:=`(organism = "E. coli"      , recon.method = "modelSEED", metric = "Reactions", nr = ecol.ms@react_num)]
dt.stat[37, `:=`(organism = "B. subtilis"  , recon.method = "modelSEED", metric = "Reactions", nr = bsub.ms@react_num)]
dt.stat[38, `:=`(organism = "P. aeruginosa", recon.method = "modelSEED", metric = "Reactions", nr = paer.ms@react_num)]
dt.stat[39, `:=`(organism = "S. oneidensis", recon.method = "modelSEED", metric = "Reactions", nr = sone.ms@react_num)]
dt.stat[40, `:=`(organism = "M. genitalium", recon.method = "modelSEED", metric = "Reactions", nr = mgen.ms@react_num)]

# Metabolites
dt.stat[41, `:=`(organism = "E. coli"      , recon.method = "curated", metric = "Metabolites", nr = ecol@met_num)]
dt.stat[42, `:=`(organism = "B. subtilis"  , recon.method = "curated", metric = "Metabolites", nr = bsub@met_num)]
dt.stat[43, `:=`(organism = "P. aeruginosa", recon.method = "curated", metric = "Metabolites", nr = paer@met_num)]
dt.stat[44, `:=`(organism = "S. oneidensis", recon.method = "curated", metric = "Metabolites", nr = sone@met_num)]
dt.stat[45, `:=`(organism = "M. genitalium", recon.method = "curated", metric = "Metabolites", nr = mgen@met_num)]

dt.stat[46, `:=`(organism = "E. coli"      , recon.method = "CarveMe", metric = "Metabolites", nr = ecol.cm@met_num)]
dt.stat[47, `:=`(organism = "B. subtilis"  , recon.method = "CarveMe", metric = "Metabolites", nr = bsub.cm@met_num)]
dt.stat[48, `:=`(organism = "P. aeruginosa", recon.method = "CarveMe", metric = "Metabolites", nr = paer.cm@met_num)]
dt.stat[49, `:=`(organism = "S. oneidensis", recon.method = "CarveMe", metric = "Metabolites", nr = sone.cm@met_num)]
dt.stat[50, `:=`(organism = "M. genitalium", recon.method = "CarveMe", metric = "Metabolites", nr = mgen.cm@met_num)]

dt.stat[51, `:=`(organism = "E. coli"      , recon.method = "gapseq", metric = "Metabolites", nr = ecol.gs@met_num)]
dt.stat[52, `:=`(organism = "B. subtilis"  , recon.method = "gapseq", metric = "Metabolites", nr = bsub.gs@met_num)]
dt.stat[53, `:=`(organism = "P. aeruginosa", recon.method = "gapseq", metric = "Metabolites", nr = paer.gs@met_num)]
dt.stat[54, `:=`(organism = "S. oneidensis", recon.method = "gapseq", metric = "Metabolites", nr = sone.gs@met_num)]
dt.stat[55, `:=`(organism = "M. genitalium", recon.method = "gapseq", metric = "Metabolites", nr = mgen.gs@met_num)]

dt.stat[56, `:=`(organism = "E. coli"      , recon.method = "modelSEED", metric = "Metabolites", nr = ecol.ms@met_num)]
dt.stat[57, `:=`(organism = "B. subtilis"  , recon.method = "modelSEED", metric = "Metabolites", nr = bsub.ms@met_num)]
dt.stat[58, `:=`(organism = "P. aeruginosa", recon.method = "modelSEED", metric = "Metabolites", nr = paer.ms@met_num)]
dt.stat[59, `:=`(organism = "S. oneidensis", recon.method = "modelSEED", metric = "Metabolites", nr = sone.ms@met_num)]
dt.stat[60, `:=`(organism = "M. genitalium", recon.method = "modelSEED", metric = "Metabolites", nr = mgen.ms@met_num)]

dt.stat[, recon.method := factor(recon.method, levels = c("curated","gapseq","CarveMe","modelSEED"))]
dt.stat[, organism := factor(organism, levels = c("E. coli","B. subtilis","P. aeruginosa","S. oneidensis","M. genitalium"))]
dt.stat[, metric := factor(metric, levels = c("Genes","Reactions","Metabolites"))]

p <- ggplot(dt.stat, aes(organism, nr, by = organism)) +
  geom_bar(col = "black", aes(fill=recon.method), stat = "identity",position = "dodge") + facet_grid(metric~., scales = "free_y") +
  scale_fill_manual(values = c("#000000", "#e41a1c", "#377eb8", "#FFB200")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "italic"))
p
ggsave("plots/model.stats.eps", p, width = 6, height = 3.5)

