library(sybilSBML)
#library(sybil.tools)
source("sybil_toolkit.R")
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1


prod.mets <- list()
#mtf.flux  <- list()
growth    <- list()


for(i in 1:nrow(ft.orgs)) {
  mod <- readRDS(paste0("models/gapseq/",ft.orgs[i,id],"_genomic.RDS"))
  mod <- constrain.mod(mod, mediadb = "media/media.tsv", media.id = ft.orgs[i,media])
  
  growth[[ft.orgs[i, id]]] <- optimizeProb(mod)@lp_obj
  
  prod.mets[[ft.orgs[i, id]]] <- get.produced.metabolites(mod)
  #[[ft.orgs[i, id]]]  <- optimizeProb_SW(mod, exclude.unused = F)
}

gs.fermprod <- rbindlist(prod.mets, idcol = T)

for(i in 1:length(growth)) {
  gs.fermprod[.id == names(growth[i]), gr.rate := growth[[i]]]
}
gs.fermprod[, recon.method := "gapseq"]

