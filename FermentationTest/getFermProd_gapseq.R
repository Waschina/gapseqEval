library(sybilSBML)
library(sybil.tools)
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1


prod.mets <- list()
mtf.flux  <- list()
growth    <- list()


for(i in 1:nrow(ft.orgs)) {
  mod <- readRDS(paste0("models/gapseq/",ft.orgs[i,id],"_genomic.RDS"))
  mod <- changeBounds(mod, react = "EX_cpd10516_e0", lb = -0.5)
  
  growth[[ft.orgs[i, id]]] <- optimizeProb(mod)@lp_obj
  
  prod.mets[[ft.orgs[i, id]]] <- get.produced.metabolites(mod)
  mtf.flux[[ft.orgs[i, id]]]  <- optimizeProb_SW(mod, exclude.unused = F)
}

gs.fermprod <- rbindlist(prod.mets, idcol = T)

for(i in 1:length(growth)) {
  gs.fermprod[.id == names(growth[i]), gr.rate := growth[[i]]]
}
gs.fermprod[, recon.method := "gapseq"]

