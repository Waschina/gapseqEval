# Function to contrain a model and which automatically detects biochemistry namespace
constrain.mod <- function(model, media.id, mediadb = NA, nasp = "auto") {
  require(data.table)
  mod <- model
  nasp = "auto"
  if(nasp == "auto")
    nasp <- get.mod.namespace(mod)
  
  if(is.na(mediadb)[1])
    mediadb <- "~/./media.tsv"
  #mediadb <- "/nuuk/Resources/diets/media.tsv"
  
  mediadb <- fread(mediadb)
  mediadb <- mediadb[medium == media.id]
  if(nrow(mediadb) == 0) {
    stop(paste("The medium with the ID",meida.id,"is not (yet) supported or part of the Medium-Database."))
  }
  
  # 1. block all inflows
  mod@lowbnd[grep("^EX_", mod@react_id, fixed = F)] <- 0
  
  # 2. contrain model according to medium
  if(nasp == "bigg")
    mediadb[,ex.rxns := paste0("EX_",compound,"(e)")]
  if(nasp == "seed")
    mediadb[,ex.rxns := paste0("EX_",modelseed,"_e0")]
  if(nasp == "vmh")
    mediadb[,ex.rxns := paste0("EX_",vmh.id,"(e)")]
  
  mediadb$mod.rxn.id <- match(mediadb$ex.rxns,mod@react_id)
  
  # Metabolited alread present
  media1 <- copy(mediadb[!is.na(mod.rxn.id)])
  mod@lowbnd[media1$mod.rxn.id] <- -media1$maxFlux
  
  return(mod)
}

# detect biochemistry namespace
get.mod.namespace <- function(mod) {
  nasp <- NA
  if(any(grepl("cpd00002",mod@met_id))) {
    nasp = "seed"
    cat("Namespace: SEED\n")
  }
  else if(any(grepl("atp",mod@met_id))) {
    if(any(grepl("^ala_L|^arg_L|^asp_L|^cys_L|^glu_L|^glc_D",mod@met_id))) {
      nasp = "vmh"
      cat("Namespace: VMH\n")
    }
    else if(any(grepl("^ala__L|^arg__L|^asp__L|^cys__L|^glu__L|^glc__D",mod@met_id))) {
      nasp = "bigg"
      cat("Namespace: BIGG\n")
    }
  }
  
  return(nasp)
}

# get metabolic by-products
get.produced.metabolites <- function(mod) {
  
  # get MTF solution
  sol.mtf <- optimizeProb(mod, algorithm = "mtf")
  dt.mtf  <- data.table(ex = mod@react_id,
                        mtf.flux = sol.mtf@fluxdist@fluxes[1:mod@react_num])
  dt.mtf.tmp <- copy(dt.mtf[grepl("^EX_", ex)])
  
  # this following two lines are there to prevent the case that a nutrient (e.g. L-Lactate)
  # from the environment is taken up, and thus enables the production of D-Lactate.
  dt.mtf.tmp[mtf.flux > 0, mtf.flux := 0]
  model.tmp <- changeBounds(mod, react = dt.mtf.tmp$ex, lb = dt.mtf.tmp$mtf.flux)
  
  
  # get FV solution
  sol.fv <- fluxVar(model.tmp, react = mod@react_id[grep("^EX_", mod@react_id)])
  
  dt <- data.table(ex = rep(mod@react_id[grep("^EX_", mod@react_id)],2),
                   rxn.name = rep(mod@react_name[grep("^EX_", mod@react_id)],2),
                   dir = c(rep("l",length(grep("^EX_", mod@react_id))),rep("u",length(grep("^EX_", mod@react_id)))),
                   fv = sol.fv@lp_obj)
  dt <- dcast(dt, ex + rxn.name ~ dir, value.var = "fv")[(u>1e-4 & l >= 0) | (u > 1)]
  
  
  
  dt <- merge(dt, dt.mtf, by = "ex")
  
  return(dt[order(-u)])
}

