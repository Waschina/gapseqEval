library(sybilSBML)
#library(sybil.tools)
source("sybil_toolkit.R")
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1


prod.mets <- list()
#mtf.flux  <- list()
growth    <- list()

for(i in 1:nrow(ft.orgs)) {
  tmp.name <- sub("^(.*?_.*?)_.*", "\\1", ft.orgs[i,id])
  tmp.name <- gsub("\\.","_",tmp.name)
  mod <- readSBMLmod(paste0("models/modelseed/",tmp.name,".sbml"))
  mod <- constrain.mod(mod, mediadb = "media/media.tsv", media.id = ft.orgs[i,media])
  
  growth[[ft.orgs[i,id]]] <- optimizeProb(mod)@lp_obj
  
  prod.mets[[ft.orgs[i,id]]] <- get.produced.metabolites(mod)
  #mtf.flux[[ft.orgs[i,id]]]  <- optimizeProb_SW(mod, exclude.unused = F)
}

ms.fermprod <- rbindlist(prod.mets, idcol = T)

for(i in 1:length(growth)) {
  ms.fermprod[.id == names(growth[i]), gr.rate := growth[[i]]]
}

ms.fermprod[, recon.method := "modelseed"]