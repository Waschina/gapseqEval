library(sybilSBML)
library(sybil.tools)
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1


prod.mets <- list()
mtf.flux  <- list()
growth    <- list()

for(i in 1:nrow(ft.orgs)) {
  mod <- readSBMLmod(paste0("models/carveme/",ft.orgs[i,id],".xml"))
  mod <- constrain.mod(mod, mediadb = "media/media.tsv", media.id = "FT")
  mod <- changeBounds(mod, react = "EX_fe3(e)", lb = -0.5)
  
  growth[[ft.orgs[i, id]]] <- optimizeProb(mod)@lp_obj
  
  prod.mets[[ft.orgs[i, id]]] <- get.produced.metabolites(mod)
  mtf.flux[[ft.orgs[i, id]]]  <- optimizeProb_SW(mod, exclude.unused = F)
}

cm.fermprod <- rbindlist(prod.mets, idcol = T)

for(i in 1:length(growth)) {
  cm.fermprod[.id == names(growth[i]), gr.rate := growth[[i]]]
}

cm.fermprod[, recon.method := "carveme"]
